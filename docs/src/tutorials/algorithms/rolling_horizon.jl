# Copyright (c) 2024 Diego Tejada and contributors                               #src
#                                                                                #src
# Permission is hereby granted, free of charge, to any person obtaining a copy   #src
# of this software and associated documentation files (the "Software"), to deal  #src
# in the Software without restriction, including without limitation the rights   #src
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell      #src
# copies of the Software, and to permit persons to whom the Software is          #src
# furnished to do so, subject to the following conditions:                       #src
#                                                                                #src
# The above copyright notice and this permission notice shall be included in all #src
# copies or substantial portions of the Software.                                #src
#                                                                                #src
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR     #src
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,       #src
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE    #src
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER         #src
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,  #src
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE  #src
# SOFTWARE.                                                                      #src

# # Rolling horizon problems
#
# **This tutorial was originally contributed by Diego Tejada.**
#
# The term "rolling horizon" refers to solving a time-dependent model
# repeatedly, where the planning interval is shifted forward in time during each
# solution step.
#
# With the features provided in [ParametricOptInterface.jl](@ref), setting up
# such a model is quite straightforward. This tutorial explains the necessary
# steps to implement a basic model with a rolling horizon in a simple generation
# expansion problem (GEP).

# ## Required packages
#
# This tutorial uses the following packages

using JuMP
import CSV
import DataFrames
import HiGHS
import ParametricOptInterface as POI
import Plots

# ## The optimization model
#
# The model is a simplified GEP problem in which we decide the new capacity in
# renewables for a power system with a given thermal and storage capacity.
#
# ### [Variables](@id rolling_horizon_variables)
#
# - Investment: $i \geq 0$
# - Renewable production: $r_t \geq 0$
# - Thermal production: $0 \leq p_t \leq \overline{P}$
# - Storage level: $0 \leq s_t \leq \overline{S}$
# - Storage charging: $0 \leq c_t \leq \overline{C}$
# - Storage discharging: $0 \leq d_t \leq \overline{D}$
#
# ### Parameters that will change each window
#
# - Demand at time $t$: $D_t$
# - Renewable availability at time $t$: $A_t$
# - Initial storage: $S_0$
#
# ### [Constraints](@id rolling_horizon_constraints)
#
# 1. **Balance Constraint:**
#
#    $p_t + r_t + d_t = D_t + c_t, \quad \forall t$
#
# 2. **Storage Dynamics for $t \geq 2$:**
#
#    $s_t = s_{t-1} + \eta^c \cdot c_t - \frac{d_t}{\eta^d}, \quad \forall t \in \{2, \ldots, \mathcal{T}\}$
#
# 3. **Initial Storage:**
#
#    $s_1 = S_0 +\eta^c \cdot c_1 - \frac{d_1}{\eta^d}$
#
# 4. **Maximum Renewable Availability:**
#
#    $r_t \leq A_t \cdot i, \quad \forall t$
#
# ### Objective Function
#
# The objective function to minimize is the total cost:
#
# $\min \left(I \cdot i + \sum_{t} O \cdot p_t\right)$

# In large-scale optimization problems, this model is solved in two steps:
#
# 1. *Investment decisions step*: This involves a simplified version of the
#    model, for example, without integer variables and representative periods.
# 2. *Operation decisions step*: After determining the values of the investments
#    from the previous step, this step involves solving an operational problem
#    to decide on production, storage levels, charging, and discharging.
#
# The second step is also computationally intensive. A common practice is to use
# the rolling horizon idea to solve multiple identical problems of a smaller
# size. These problems differ only in parameters such as demand, renewable
# profiles, or initial conditions.
#
# This example focuses on the second step, aiming to determine the operational
# variables for a given investment using a rolling horizon strategy.
#
# ## Parameter definition and input data

# There are two main parameters for a rolling horizon basic implementation: the
# optimization window and the move forward.

# **Optimization Window**: this value defines how many periods (for example,
# hours) we will optimize each time. For this example, we set the default value
# to 48 hours, meaning we will optimize two days each time.

optimization_window = 48

# **Move Forward**: this value defines how many periods (for example, hours) we
# will move forward to optimize the next optimization window. For this example,
# we set the default value in 24 hours, meaning we will move one day ahead each
# time.

move_forward = 24

# Note that the move forward parameter must be lower or equal to the
# optimization window parameter to work correctly.

@assert optimization_window >= move_forward

# Let's explore the input data in file [rolling_horizon.csv](rolling_horizon.csv).
# We have a total time horizon of a week (that is, 168 hours), an electricity
# demand, and a solar production profile.

filename = joinpath(@__DIR__, "rolling_horizon.csv")
time_series = CSV.read(filename, DataFrames.DataFrame);

# We define the solar investment (for example, 150 MW) to determine the solar
# production during the operation optimization step.

solar_investment = 150
time_series.solar_MW = solar_investment * time_series.solar_pu

# In addition, we can determine some basic information about the rolling
# horizon, such as the number of windows that we are going to optimize given the
# problem's time horizon.

total_time_length = size(time_series, 1)

# The total number of time windows we will solve for is:

ceil(Int, total_time_length / move_forward)

# Finally, we can see a plot representing the first two optimization windows and
# the move forward parameter to have a better idea of how the rolling horizon
# works.

x_series = 1:total_time_length
y_series = [time_series.demand_MW, time_series.solar_MW]
plot_1 = Plots.plot(x_series, y_series; label = ["demand" "solar"])
plot_2 = Plots.plot(x_series, y_series; label = false)
window = [0, optimization_window]
Plots.vspan!(plot_1, window; alpha = 0.25, label = false)
Plots.vspan!(plot_2, move_forward .+ window; alpha = 0.25, label = false)
text_1 = Plots.text("optimization\n  window 1", :top, :left, 8)
Plots.annotate!(plot_1, 18, time_series.solar_MW[12], text_1)
text_2 = Plots.text("optimization\n  window 2", :top, :left, 8)
Plots.annotate!(plot_2, 42, time_series.solar_MW[12], text_2)
Plots.plot(
    plot_1,
    plot_2;
    layout = (2, 1),
    linewidth = 3,
    xticks = 0:12:total_time_length,
    xlabel = "Hours",
    ylabel = "MW",
)

# ## JuMP model

# We have all the information we need to create a JuMP model to solve a single
# window of our rolling horizon problem.

model = Model(() -> POI.Optimizer(HiGHS.Optimizer()))
set_silent(model)
@variables(model, begin
    i == solar_investment
    0 <= r[1:optimization_window]
    0 <= p[1:optimization_window] <= 150
    0 <= s[1:optimization_window] <= 40
    0 <= c[1:optimization_window] <= 10
    0 <= d[1:optimization_window] <= 10
    ## Initialize empty parameters. These values will get updated layer
    D[t in 1:optimization_window] in Parameter(0)
    A[t in 1:optimization_window] in Parameter(0)
    So in Parameter(0)
end)
@objective(model, Min, 100 * i + 50 * sum(p))
@constraints(
    model,
    begin
        p .+ r .+ d .== D .+ c
        s[1] == So + 0.9 * c[1] - d[1] / 0.9
        [t in 2:optimization_window], s[t] == s[t-1] + 0.9 * c[t] - d[t] / 0.9
        r .<= A .* i
    end
)
model

# After the optimization, we can store the results in vectors. It's important to
# note that despite optimizing for 48 hours (the default value), we only store
# the values for the "move forward" parameter (for example, 24 hours or one day
# using the default value). This approach ensures that there is a buffer of
# additional periods or hours beyond the "move forward" parameter to prevent the
# storage from depleting entirely at the end of the specified hours.

objective_function_per_window = Float64[]
renewable_production = Float64[]
storage_level = Float64[0.0]  # Include an initial storage level

# Now we can iterate across the windows of our rolling horizon problem, and at
# each window, we:

# 1. update the parameters in the models
# 2. solve the model for that window
# 3. store the results for later analysis

for offset in 0:move_forward:total_time_length-1
    ## Step 1: update the parameter values over the optimization_window
    for t in 1:optimization_window
        ## This row computation just let's us "wrap around" the `time_series`
        ## DataFrame, so that the forecase for demand and solar PU in day 8 is
        ## the same as day 1. In real models, you might choose to do something
        ## different.
        row = 1 + mod(offset + t, size(time_series, 1))
        set_parameter_value(model[:D][t], time_series[row, :demand_MW])
        set_parameter_value(model[:A][t], time_series[row, :solar_pu])
    end
    set_parameter_value(model[:So], storage_level[end])
    ## Step 2: solve the model
    optimize!(model)
    ## Step 3: store the results of the move_forward values
    push!(objective_function_per_window, objective_value(model))
    for t in 1:move_forward
        push!(renewable_production, value(model[:r][t]))
        push!(storage_level, value(model[:s][t]))
    end
end

# We can explore the outputs in the following graphs:

Plots.plot(
    objective_function_per_window ./ 10^3;
    label = false,
    linewidth = 3,
    xlabel = "Window",
    ylabel = "[000'] \$",
)

#-

Plots.plot(
    [time_series.demand_MW, renewable_production, storage_level[2:end]];
    label = ["demand" "solar" "battery"],
    linewidth = 3,
    xlabel = "Hours",
    ylabel = "MW",
    xticks = 0:12:total_time_length,
)

# ## Final remark

# [ParametricOptInterface.jl](@ref) offers an easy way to update the parameters
# of an optimization problem that will be solved several times, as in the
# rolling horizon implementation. It has the benefit of avoiding rebuilding the
# model each time we want to solve it with new information in a new window.

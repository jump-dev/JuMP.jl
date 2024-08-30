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
# The purpose of this tutorial is to demonstrate how to use [ParametricOptInterface.jl](@ref)
# to solve a rolling horizon optimization problem.
#
# The term "rolling horizon" refers to solving a time-dependent model
# repeatedly, where the planning interval is shifted forward in time during each
# solution step.
#
# As a motivating example, this tutorial models the operations of a power system
# with solar generation and a battery.

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
# The model is a simplified model of a power system's operations with battery
# storage.
#
# We model the system of a set of time-steps $t \in 1,\ldots,T$, where each time
# step is a period of one hour.
#
# There are five types of decision variables in the model:
#
# - Renewable production: $r_t \geq 0$
# - Thermal production: $0 \leq p_t \leq \overline{P}$
# - Storage level: $0 \leq s_t \leq \overline{S}$
# - Storage charging: $0 \leq c_t \leq \overline{C}$
# - Storage discharging: $0 \leq d_t \leq \overline{D}$
#
# For the purpose of this tutorial, there are three parameters of interest:
#
# - Demand at time $t$: $D_t$
# - Renewable availability at time $t$: $A_t$
# - Initial storage: $S_0$
#
# The objective function to minimize is the total cost of thermal generation:
# $$\min \sum_{t} O \cdot p_t$$
#
# For the constraints, we must balance power generation and consumption in all
# time periods:
# $$p_t + r_t + d_t = D_t + c_t, \quad \forall t$$
#
# We need to account for the dynamics of the battery storage:
# $$s_t = s_{t-1} + \eta^c \cdot c_t - \frac{d_t}{\eta^d}, \quad \forall t$$
# with the boundary condition that $s_0 = S_0$.
#
# Finally, the level of renewable energy production is limited by the
# availability factor $A$ and the installed capacity $i$:
# $$r_t \leq A_t \cdot i, \quad \forall t$$

# Solving this problem with a large number of time steps is computationally
# challenging. A common practice is to use the rolling horizon idea to solve
# multiple identical problems of a smaller size. These problems differ only in
# parameters such as demand, renewable availability, and initial storage. By
# combining the solution of many smaller problems, we can recover a feasible
# solution to the full problem. However, because we don't optimize the full set
# of decisions in a single optimization problem, the recovered solution might be
# suboptimal.

# ## Parameter definition and input data

# There are two main parameters for a rolling horizon implementation: the
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

# We multiple the level of solar investment by the time series of availability
# to get actual MW generated.

time_series.solar_MW = solar_investment * time_series.solar_pu;

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
    0 <= r[1:optimization_window]
    0 <= p[1:optimization_window] <= 150
    0 <= s[1:optimization_window] <= 40
    0 <= c[1:optimization_window] <= 10
    0 <= d[1:optimization_window] <= 10
    ## Initialize empty parameters. These values will get updated layer
    D[t in 1:optimization_window] in Parameter(0)
    A[t in 1:optimization_window] in Parameter(0)
    S_0 in Parameter(0)
end)
@objective(model, Min, 50 * sum(p))
@constraints(
    model,
    begin
        p .+ r .+ d .== D .+ c
        s[1] == S_0 + 0.9 * c[1] - d[1] / 0.9
        [t in 2:optimization_window], s[t] == s[t-1] + 0.9 * c[t] - d[t] / 0.9
        r .<= A .* solar_investment
    end
)
model

# After the optimization, we can store the results in vectors. It's important to
# note that despite optimizing for 48 hours (the default value), we only store
# the values for the "move forward" parameter (for example, 24 hours or one day
# using the default value). This approach ensures that there is a buffer of
# additional periods or hours beyond the "move forward" parameter to prevent the
# storage from depleting entirely at the end of the specified hours.

renewable_production = Float64[]
storage_level = Float64[0.0]  # Include an initial storage level

# We'll also plot the solution at each of the time-steps to help visualize the
# solution to the rolling horizon problems.

function plot_solution(model, offset)
    plot = Plots.plot(;
        ylabel = "MW",
        xlims = (0, total_time_length),
        xticks = 0:12:total_time_length,
    )
    x = offset .+ (1:optimization_window)
    y = hcat(value.(model[:p]), value.(model[:r]), value.(model[:d]))
    if offset == 0
        Plots.areaplot!(x, y; label = ["thermal" "solar" "discharge"])
        Plots.areaplot!(x, -value.(model[:c]); label = "charge")
    else
        Plots.areaplot!(x, y; label = false)
        Plots.areaplot!(x, -value.(model[:c]); label = false)
    end
    return plot
end
plots = Any[]

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
    set_parameter_value(model[:S_0], storage_level[end])
    ## Step 2: solve the model
    optimize!(model)
    ## Step 3: store the results of the move_forward values
    for t in 1:move_forward
        push!(renewable_production, value(model[:r][t]))
        push!(storage_level, value(model[:s][t]))
    end
    push!(plots, plot_solution(model, offset))
end

# We can now plot the solution to the week-long problem:

Plots.plot(
    [time_series.demand_MW, renewable_production, storage_level[2:end]];
    label = ["demand" "solar" "battery"],
    linewidth = 3,
    xlabel = "Hours",
    ylabel = "MW",
    xticks = 0:12:total_time_length,
)

# and visualize each of the rolling horizon subplots:

Plots.plot(
    plots...;
    layout = (length(plots), 1),
    size = (600, 800),
    margin = 3Plots.mm,
)

# ## Final remark

# [ParametricOptInterface.jl](@ref) offers an easy way to update the parameters
# of an optimization problem that will be solved several times, as in the
# rolling horizon implementation. It has the benefit of avoiding rebuilding the
# model each time we want to solve it with new information in a new window.

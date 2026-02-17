# Copyright (c) 2019 Yury Dvorkin, Miles Lubin, and contributors                #src
#                                                                               #src
# This work is licensed under a Creative Commons Attribution-ShareAlike 4.0     #src
# International License. URL: http://creativecommons.org/licenses/by-sa/4.0.    #src

# # Power Systems

# **Originally Contributed by**: Yury Dvorkin and Miles Lubin

# This tutorial demonstrates how to formulate basic power systems engineering
# models in JuMP.

# We will consider basic "economic dispatch" and "unit commitment" models
# without taking into account transmission constraints.

# For this tutorial, we use the following packages:

using JuMP
import DataFrames
import HiGHS
import Plots
import StatsPlots

# ## Economic dispatch

# Economic dispatch (ED) is an optimization problem that minimizes the cost of
# supplying energy demand subject to operational constraints on power system
# assets. In its simplest modification, ED is an LP problem solved for an
# aggregated load and wind forecast and for a single infinitesimal moment.

# Mathematically, the ED problem can be written as follows:

# ```math
# \min \sum_{i \in I} c^g_{i} \cdot g_{i} + c^w \cdot w,
# ```

# where $c_{i}$ and $g_{i}$ are the incremental cost (\$/MWh) and power
# output (MW) of the $i^{th}$ generator, respectively, and $c^w$ and $w$ are the
# incremental cost (\$/MWh) and wind power injection (MW), respectively.

# Subject to the constraints:

# * Minimum ($g^{\min}$) and maximum ($g^{\max}$) limits on power outputs of
#   generators:
#   $g^{\min}_{i} \leq g_{i} \leq g^{\max}_{i}.$
# * Constraint on the wind power injection:
#   $0 \leq w \leq w^f,$
#   where $w$ and $w^f$ are the wind power injection and wind power forecast,
#   respectively.
# * Power balance constraint:
#   $\sum_{i \in I} g_{i} + w = d^f,$
#   where $d^f$ is the demand forecast.

# Further reading on ED models can be found in A. J. Wood, B. F. Wollenberg, and
# G. B. Sheblé, "Power Generation, Operation and Control", Wiley, 2013.

# Define some input data about the test system.

# We define some thermal generators:

function ThermalGenerator(
    min::Float64,
    max::Float64,
    fixed_cost::Float64,
    variable_cost::Float64,
)
    return (
        min = min,
        max = max,
        fixed_cost = fixed_cost,
        variable_cost = variable_cost,
    )
end

generators = [
    ThermalGenerator(0.0, 1000.0, 1000.0, 50.0),
    ThermalGenerator(300.0, 1000.0, 0.0, 100.0),
]

# A wind generator

WindGenerator(variable_cost::Float64) = (variable_cost = variable_cost,)

wind_generator = WindGenerator(50.0)

# And a scenario

function Scenario(demand::Float64, wind::Float64)
    return (demand = demand, wind = wind)
end

scenario = Scenario(1500.0, 200.0)

# Create a function `solve_ed`, which solves the economic dispatch problem for a
# given set of input parameters.

function solve_ed(generators::Vector, wind, scenario)
    ## Define the economic dispatch (ED) model
    ed = Model(HiGHS.Optimizer)
    set_silent(ed)
    ## Define decision variables
    ## power output of generators
    N = length(generators)
    @variable(ed, generators[i].min <= g[i = 1:N] <= generators[i].max)
    ## wind power injection
    @variable(ed, 0 <= w <= scenario.wind)
    ## Define the objective function
    @objective(
        ed,
        Min,
        sum(generators[i].variable_cost * g[i] for i in 1:N) +
        wind.variable_cost * w,
    )
    ## Define the power balance constraint
    @constraint(ed, sum(g[i] for i in 1:N) + w == scenario.demand)
    ## Solve statement
    optimize!(ed)
    ## return the optimal value of the objective function and its minimizers
    return (
        g = value.(g),
        w = value(w),
        wind_spill = scenario.wind - value(w),
        total_cost = objective_value(ed),
    )
end

# Solve the economic dispatch problem

solution = solve_ed(generators, wind_generator, scenario);

println("Dispatch of Generators: ", solution.g, " MW")
println("Dispatch of Wind: ", solution.w, " MW")
println("Wind spillage: ", solution.wind_spill, " MW")
println("Total cost: \$", solution.total_cost)

# ## Economic dispatch with adjustable incremental costs

# In the following exercise we adjust the incremental cost of generator G1 and
# observe its impact on the total cost.

function scale_generator_cost(g, scale)
    return ThermalGenerator(g.min, g.max, g.fixed_cost, scale * g.variable_cost)
end

start = time()
c_g_scale_df = DataFrames.DataFrame(;
    ## Scale factor
    scale = Float64[],
    ## Dispatch of Generator 1 [MW]
    dispatch_G1 = Float64[],
    ## Dispatch of Generator 2 [MW]
    dispatch_G2 = Float64[],
    ## Dispatch of Wind [MW]
    dispatch_wind = Float64[],
    ## Spillage of Wind [MW]
    spillage_wind = Float64[],
    ## Total cost [$]
    total_cost = Float64[],
)
for c_g1_scale in 0.5:0.1:3.0
    ## Update the incremental cost of the first generator at every iteration.
    new_generators = scale_generator_cost.(generators, [c_g1_scale, 1.0])
    ## Solve the ed problem with the updated incremental cost
    sol = solve_ed(new_generators, wind_generator, scenario)
    push!(
        c_g_scale_df,
        (c_g1_scale, sol.g[1], sol.g[2], sol.w, sol.wind_spill, sol.total_cost),
    )
end
print(string("elapsed time: ", time() - start, " seconds"))

#-

c_g_scale_df

# ## Modifying the JuMP model in-place

# Note that in the previous exercise we entirely rebuilt the optimization model
# at every iteration of the internal loop, which incurs an additional
# computational burden. This burden can be alleviated if instead of re-building
# the entire model, we modify a specific constraint(s) or the objective
# function, as it shown in the example below.

# Compare the computing time in case of the above and below models.

function solve_ed_inplace(
    generators::Vector,
    wind,
    scenario,
    scale::AbstractVector{Float64},
)
    obj_out = Float64[]
    w_out = Float64[]
    g1_out = Float64[]
    g2_out = Float64[]
    ## This function only works for two generators
    @assert length(generators) == 2
    ed = Model(HiGHS.Optimizer)
    set_silent(ed)
    N = length(generators)
    @variable(ed, generators[i].min <= g[i = 1:N] <= generators[i].max)
    @variable(ed, 0 <= w <= scenario.wind)
    @objective(
        ed,
        Min,
        sum(generators[i].variable_cost * g[i] for i in 1:N) +
        wind.variable_cost * w,
    )
    @constraint(ed, sum(g[i] for i in 1:N) + w == scenario.demand)
    for c_g1_scale in scale
        @objective(
            ed,
            Min,
            c_g1_scale * generators[1].variable_cost * g[1] +
            generators[2].variable_cost * g[2] +
            wind.variable_cost * w,
        )
        optimize!(ed)
        push!(obj_out, objective_value(ed))
        push!(w_out, value(w))
        push!(g1_out, value(g[1]))
        push!(g2_out, value(g[2]))
    end
    df = DataFrames.DataFrame(;
        scale = scale,
        dispatch_G1 = g1_out,
        dispatch_G2 = g2_out,
        dispatch_wind = w_out,
        spillage_wind = scenario.wind .- w_out,
        total_cost = obj_out,
    )
    return df
end

start = time()
inplace_df = solve_ed_inplace(generators, wind_generator, scenario, 0.5:0.1:3.0)
print(string("elapsed time: ", time() - start, " seconds"))

# Adjusting specific constraints or the objective function is faster than
# re-building the entire model.

inplace_df

# ## Inefficient usage of wind generators

# The economic dispatch problem does not perform commitment decisions and, thus,
# assumes that all generators must be dispatched at least at their minimum power
# output limit. This approach is not cost efficient and may lead to absurd
# decisions. For example, if $d = \sum_{i \in I} g^{\min}_{i}$, the wind power
# injection must be zero, i.e. all available wind generation is spilled, to meet
# the minimum power output constraints on generators.

# In the following example, we adjust the total demand and observed how it
# affects wind spillage.

demand_scale_df = DataFrames.DataFrame(;
    demand = Float64[],
    dispatch_G1 = Float64[],
    dispatch_G2 = Float64[],
    dispatch_wind = Float64[],
    spillage_wind = Float64[],
    total_cost = Float64[],
)

function scale_demand(scenario, scale)
    return Scenario(scale * scenario.demand, scenario.wind)
end

for demand_scale in 0.2:0.1:1.4
    new_scenario = scale_demand(scenario, demand_scale)
    sol = solve_ed(generators, wind_generator, new_scenario)
    push!(
        demand_scale_df,
        (
            new_scenario.demand,
            sol.g[1],
            sol.g[2],
            sol.w,
            sol.wind_spill,
            sol.total_cost,
        ),
    )
end

demand_scale_df

#-

dispatch_plot = StatsPlots.@df(
    demand_scale_df,
    Plots.plot(
        :demand,
        [:dispatch_G1, :dispatch_G2],
        labels = ["G1" "G2"],
        title = "Thermal Dispatch",
        legend = :bottomright,
        linewidth = 3,
        xlabel = "Demand",
        ylabel = "Dispatch [MW]",
    ),
)

wind_plot = StatsPlots.@df(
    demand_scale_df,
    Plots.plot(
        :demand,
        [:dispatch_wind, :spillage_wind],
        labels = ["Dispatch" "Spillage"],
        title = "Wind",
        legend = :bottomright,
        linewidth = 3,
        xlabel = "Demand [MW]",
        ylabel = "Energy [MW]",
    ),
)

Plots.plot(dispatch_plot, wind_plot)

# This particular drawback can be overcome by introducing binary decisions on
# the "on/off" status of generators. This model is called unit commitment and
# considered later in these notes.

# For further reading on the interplay between wind generation and the minimum
# power output constraints of generators, we refer interested readers to R.
# Baldick, "Wind and Energy Markets: A Case Study of Texas," IEEE Systems
# Journal, vol. 6, pp. 27-34, 2012.

# ## Unit commitment

# The Unit Commitment (UC) model can be obtained from ED model by introducing
# binary variable associated with each generator. This binary variable can
# attain two values: if it is "1", the generator is synchronized and, thus, can
# be dispatched, otherwise, i.e. if the binary variable is "0", that generator
# is not synchronized and its power output is set to 0.

# To obtain the mathematical formulation of the UC model, we will modify the
# constraints of the ED model as follows:

# ```math
# g^{\min}_{i} \cdot u_{t,i} \leq g_{i} \leq g^{\max}_{i} \cdot u_{t,i},
# ```
# where $u_{i} \in \{0,1\}.$ In this constraint, if $u_{i} = 0$, then
# $g_{i}  = 0$. On the other hand, if $u_{i} = 1$, then
# $g^{min}_{i} \leq g_{i} \leq g^{max}_{i}$.

# For further reading on the UC problem we refer interested readers to G.
# Morales-Espana, J. M. Latorre, and A. Ramos, "Tight and Compact MILP
# Formulation for the Thermal Unit Commitment Problem," IEEE Transactions on
# Power Systems, vol. 28, pp. 4897-4908, 2013.

# In the following example we convert the ED model explained above to the UC
# model.

function solve_uc(generators::Vector, wind, scenario)
    uc = Model(HiGHS.Optimizer)
    set_silent(uc)
    N = length(generators)
    @variable(uc, generators[i].min <= g[i = 1:N] <= generators[i].max)
    @variable(uc, 0 <= w <= scenario.wind)
    @constraint(uc, sum(g[i] for i in 1:N) + w == scenario.demand)
    ## !!! New: add binary on-off variables for each generator
    @variable(uc, u[i = 1:N], Bin)
    @constraint(uc, [i = 1:N], g[i] <= generators[i].max * u[i])
    @constraint(uc, [i = 1:N], g[i] >= generators[i].min * u[i])
    @objective(
        uc,
        Min,
        sum(generators[i].variable_cost * g[i] for i in 1:N) +
        wind.variable_cost * w +
        ## !!! new
        sum(generators[i].fixed_cost * u[i] for i in 1:N)
    )
    optimize!(uc)
    status = termination_status(uc)
    if status != OPTIMAL
        return (status = status,)
    end
    return (
        status = status,
        g = value.(g),
        w = value(w),
        wind_spill = scenario.wind - value(w),
        u = value.(u),
        total_cost = objective_value(uc),
    )
end

# Solve the economic dispatch problem
solution = solve_uc(generators, wind_generator, scenario)

println("Dispatch of Generators: ", solution.g, " MW")
println("Commitments of Generators: ", solution.u)
println("Dispatch of Wind: ", solution.w, " MW")
println("Wind spillage: ", solution.wind_spill, " MW")
println("Total cost: \$", solution.total_cost)

# ## Unit commitment as a function of demand

# After implementing the UC model, we can now assess the interplay between the
# minimum power output constraints on generators and wind generation.

uc_df = DataFrames.DataFrame(;
    demand = Float64[],
    commitment_G1 = Float64[],
    commitment_G2 = Float64[],
    dispatch_G1 = Float64[],
    dispatch_G2 = Float64[],
    dispatch_wind = Float64[],
    spillage_wind = Float64[],
    total_cost = Float64[],
)

for demand_scale in 0.2:0.1:1.4
    new_scenario = scale_demand(scenario, demand_scale)
    sol = solve_uc(generators, wind_generator, new_scenario)
    if sol.status == OPTIMAL
        push!(
            uc_df,
            (
                new_scenario.demand,
                sol.u[1],
                sol.u[2],
                sol.g[1],
                sol.g[2],
                sol.w,
                sol.wind_spill,
                sol.total_cost,
            ),
        )
    end
    println("Status: $(sol.status) for demand_scale = $(demand_scale)")
end

#-

uc_df

#-

commitment_plot = StatsPlots.@df(
    uc_df,
    Plots.plot(
        :demand,
        [:commitment_G1, :commitment_G2],
        labels = ["G1" "G2"],
        title = "Committment",
        legend = :bottomright,
        linewidth = 3,
        xlabel = "Demand [MW]",
        ylabel = "Commitment decision {0, 1}",
    ),
)

dispatch_plot = StatsPlots.@df(
    uc_df,
    Plots.plot(
        :demand,
        [:dispatch_G1, :dispatch_G2, :dispatch_wind],
        labels = ["G1" "G2" "Wind"],
        title = "Dispatch [MW]",
        legend = :bottomright,
        linewidth = 3,
        xlabel = "Demand",
        ylabel = "Dispatch [MW]",
    ),
)

Plots.plot(commitment_plot, dispatch_plot)

# ## Nonlinear economic dispatch

# As a final example, we modify our economic dispatch problem in two ways:
#
#  * The thermal cost function is user-defined
#  * The output of the wind is only the square-root of the dispatch

import Ipopt

"""
    thermal_cost_function(g)

A user-defined thermal cost function in pure-Julia! You can include
nonlinearities, and even things like control flow.

!!! warning
    It's still up to you to make sure that the function has a meaningful
    derivative.
"""
function thermal_cost_function(g)
    if g <= 500
        return g
    else
        return g + 1e-2 * (g - 500)^2
    end
end

function solve_nonlinear_ed(
    generators::Vector,
    wind,
    scenario;
    silent::Bool = false,
)
    model = Model(Ipopt.Optimizer)
    if silent
        set_silent(model)
    end
    register(model, :tcf, 1, thermal_cost_function; autodiff = true)
    N = length(generators)
    @variable(model, generators[i].min <= g[i = 1:N] <= generators[i].max)
    @variable(model, 0 <= w <= scenario.wind)
    @NLobjective(
        model,
        Min,
        sum(generators[i].variable_cost * tcf(g[i]) for i in 1:N) +
        wind.variable_cost * w,
    )
    @NLconstraint(model, sum(g[i] for i in 1:N) + sqrt(w) == scenario.demand)
    optimize!(model)
    return (
        g = value.(g),
        w = value(w),
        wind_spill = scenario.wind - value(w),
        total_cost = objective_value(model),
    )
end

solution = solve_nonlinear_ed(generators, wind_generator, scenario)

# Now let's see how the wind is dispatched as a function of the cost:

wind_cost = 0.0:1:100
wind_dispatch = Float64[]
for c in wind_cost
    sol = solve_nonlinear_ed(
        generators,
        WindGenerator(c),
        scenario;
        silent = true,
    )
    push!(wind_dispatch, sol.w)
end

Plots.plot(
    wind_cost,
    wind_dispatch;
    xlabel = "Cost",
    ylabel = "Dispatch [MW]",
    label = false,
)

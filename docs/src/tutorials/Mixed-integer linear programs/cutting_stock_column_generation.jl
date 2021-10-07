# Copyright 2017, Iain Dunning, Joey Huchette, Miles Lubin, and contributors    #src
# This Source Code Form is subject to the terms of the Mozilla Public License   #src
# v.2.0. If a copy of the MPL was not distributed with this file, You can       #src
# obtain one at https://mozilla.org/MPL/2.0/.                                   #src

# # Cutting stock

# This example solves the cutting stock problem (sometimes also called the
# cutting rod problem) using a column-generation technique. It is based on
# [https://doi.org/10.5281/zenodo.3329388](https://doi.org/10.5281/zenodo.3329388).

# Intuitively, this problem is about cutting large rolls of paper into smaller
# pieces. There is an exact demand of pieces to meet, and all rolls have the
# same size. The goal is to meet the demand while maximizing the profits (each
# paper roll has a fixed cost, each sold piece allows earning some money),
# which is roughly equivalent to using the smallest amount of rolls
# to cut (or, equivalently, to minimize the amount of paper waste).

# This function takes five parameters:

#   * `maxwidth`: the maximum width of a roll (or length of a rod)
#   * `widths`: an array of the requested widths
#   * `rollcost`: the cost of a complete roll
#   * `demand`: the demand, in number of pieces, for each width
#   * `prices`: the selling price for each width

# Mathematically, this problem might be formulated with two variables:

#   * `x[i, j] ∈ ℕ`: the number of times the width `i` is cut out of the roll `j`
#   * `y[j] ∈ 𝔹`: whether the roll `j` is used

# Several constraints are needed:

#   * the demand must be satisfied, for each width `i`:
#     ∑j x[i, j] = demand[i]
#   * the roll size cannot be exceed, for each roll `j` that is used:
#     ∑i x[i, j] width[i] ≤ maxwidth y[j]

# If you want to implement this naïve model, you will need an upper bound on the
# number of rolls to use: the simplest one is to consider that each required
# width is cut from its own roll, i.e. `j` varies from 1 to ∑i demand[i].

# This example prefers a more advanced technique to solve this problem:
# column generation.

# It considers a different set of variables: patterns of width to cut a roll.
# The decisions then become the number of times each pattern is used (i.e. the
# number of rolls that are cut following this pattern).

# The intelligence comes from the way these patterns are chosen: not all of them
# are considered, but only the "interesting" ones, within the master problem.

# A "pricing" problem is used to decide whether a new pattern should be
# generated or not (it is implemented in the function `solve_pricing`).
# "Interesting" means, for a pattern, that the optimal solution may use this
# cutting pattern.

# In more detail, the solving process is the following. First, a series of dumb
# patterns are generated (just one width per roll, repeated until the roll is
# completely cut). Then, the master problem is solved with these first patterns
# and its dual solution is passed on to the pricing problem. The latter decides
# if there is a new pattern to include in the formulation or not; if so,
# it returns it to the master problem. The master is solved again, the new dual
# variables are given to the pricing problem, until there is no more pattern to
# generate from the pricing problem: all "interesting" patterns have been
# generated, and the master can take its optimal decision.

# In the implementation, the variables deciding how many times a pattern is
# chosen are called `θ`.

# For more information on column-generation techniques applied on the cutting
# stock problem, you can see:

#   * [Integer programming column generation strategies for the cutting stock
#     problem and its variants](https://tel.archives-ouvertes.fr/tel-00011657/document)
#   * [Tackling the cutting stock problem](https://doi.org/10.5281/zenodo.3329388)

# This example uses the following packages:

using JuMP
import GLPK
import SparseArrays
import Test  #src

# The function `solve_pricing` implements the pricing problem for the function
# `example_cutting_stock`.
#
# It takes, as input, the dual solution from the master problem and the cutting
# stock instance.
#
# It outputs either a new cutting pattern, or `nothing` if no pattern could
# improve the current cost.

function solve_pricing(
    dual_demand_satisfaction,
    maxwidth,
    widths,
    rollcost,
    demand,
    prices,
)
    reduced_costs = dual_demand_satisfaction + prices
    n = length(reduced_costs)
    ## The actual pricing model.
    submodel = Model(GLPK.Optimizer)
    set_silent(submodel)
    @variable(submodel, xs[1:n] >= 0, Int)
    @constraint(submodel, sum(xs .* widths) <= maxwidth)
    @objective(submodel, Max, sum(xs .* reduced_costs))
    optimize!(submodel)
    new_pattern = round.(Int, value.(xs))
    net_cost =
        rollcost - sum(new_pattern .* (dual_demand_satisfaction .+ prices))
    ## If the net cost of this new pattern is nonnegative, no more patterns to add.
    return net_cost >= 0 ? nothing : new_pattern
end

function example_cutting_stock(; max_gen_cols::Int = 5_000)
    maxwidth = 100.0
    rollcost = 500.0
    prices = [
        167.0,
        197.0,
        281.0,
        212.0,
        225.0,
        111.0,
        93.0,
        129.0,
        108.0,
        106.0,
        55.0,
        85.0,
        66.0,
        44.0,
        47.0,
        15.0,
        24.0,
        13.0,
        16.0,
        14.0,
    ]
    widths = [
        75.0,
        75.0,
        75.0,
        75.0,
        75.0,
        53.8,
        53.0,
        51.0,
        50.2,
        32.2,
        30.8,
        29.8,
        20.1,
        16.2,
        14.5,
        11.0,
        8.6,
        8.2,
        6.6,
        5.1,
    ]
    demand = [
        38,
        44,
        30,
        41,
        36,
        33,
        36,
        41,
        35,
        37,
        44,
        49,
        37,
        36,
        42,
        33,
        47,
        35,
        49,
        42,
    ]
    nwidths = length(prices)
    n = length(widths)
    ncols = length(widths)
    ## Initial set of patterns (stored in a sparse matrix: a pattern won't
    ## include many different cuts).
    patterns = SparseArrays.spzeros(UInt16, n, ncols)
    for i in 1:n
        patterns[i, i] =
            min(floor(Int, maxwidth / widths[i]), round(Int, demand[i]))
    end
    ## Write the master problem with this "reduced" set of patterns.
    ## Not yet integer variables: otherwise, the dual values may make no sense
    ## (actually, GLPK will yell at you if you're trying to get duals for
    ## integer problems).
    m = Model(GLPK.Optimizer)
    set_silent(m)
    @variable(m, θ[1:ncols] >= 0)
    @objective(
        m,
        Min,
        sum(
            θ[p] * (rollcost - sum(patterns[j, p] * prices[j] for j in 1:n)) for
            p in 1:ncols
        )
    )
    @constraint(
        m,
        demand_satisfaction[j = 1:n],
        sum(patterns[j, p] * θ[p] for p in 1:ncols) >= demand[j]
    )
    ## First solve of the master problem.
    optimize!(m)
    if termination_status(m) != MOI.OPTIMAL
        warn("Master not optimal ($ncols patterns so far)")
    end
    ## Then, generate new patterns, based on the dual information.
    while ncols - n <= max_gen_cols ## Generate at most max_gen_cols columns.
        if !has_duals(m)
            break
        end
        new_pattern = solve_pricing(
            dual.(demand_satisfaction),
            maxwidth,
            widths,
            rollcost,
            demand,
            prices,
        )
        ## No new pattern to add to the formulation: done!
        if new_pattern === nothing
            break
        end
        ## Otherwise, add the new pattern to the master problem, recompute the
        ## duals, and go on waltzing one more time with the pricing problem.
        ncols += 1
        patterns = hcat(patterns, new_pattern)
        ## One new variable.
        push!(θ, @variable(m, base_name = "θ", lower_bound = 0))
        ## Update the objective function.
        set_objective_coefficient(
            m,
            θ[ncols],
            rollcost - sum(patterns[j, ncols] * prices[j] for j in 1:n),
        )
        ## Update the constraint number j if the new pattern impacts this production.
        for j in 1:n
            if new_pattern[j] > 0
                set_normalized_coefficient(
                    demand_satisfaction[j],
                    θ[ncols],
                    new_pattern[j],
                )
            end
        end
        ## Solve the new master problem to update the dual variables.
        optimize!(m)
        if termination_status(m) != MOI.OPTIMAL
            @warn("Master not optimal ($ncols patterns so far)")
        end
    end
    ## Finally, impose the master variables to be integer and resolve.
    ## To be exact, at each node in the branch-and-bound tree, we would need to
    ## restart the column generation process (just in case a new column would be
    ## interesting to add). This way, we only get an upper bound (a feasible
    ## solution).
    set_integer.(θ)
    optimize!(m)
    if termination_status(m) != MOI.OPTIMAL
        @warn("Final master not optimal ($ncols patterns)")
        return
    end
    Test.@test objective_value(m) ≈ 78374.0 atol = 1e-3  #src
    println("Final solution:")
    for i in 1:length(θ)
        if value(θ[i]) > 0.5
            println("$(round(Int, value(θ[i]))) units of pattern $(i)")
        end
    end
    return
end

example_cutting_stock()

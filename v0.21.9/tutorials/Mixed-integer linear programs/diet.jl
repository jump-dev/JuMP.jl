# Copyright 2017, Iain Dunning, Joey Huchette, Miles Lubin, and contributors    #src
# This Source Code Form is subject to the terms of the Mozilla Public License   #src
# v.2.0. If a copy of the MPL was not distributed with this file, You can       #src
# obtain one at https://mozilla.org/MPL/2.0/.                                   #src

# # The diet problem

# Solve the classic "diet problem", also known as the [Stigler diet](https://en.wikipedia.org/wiki/Stigler_diet).
# The code is based on [an example from Gurobi](https://www.gurobi.com/documentation/9.0/examples/diet_cpp_cpp.html).

using JuMP
import GLPK
import Test

function print_solution(is_optimal, foods, buy)
    println("RESULTS:")
    if is_optimal
        for food in foods
            println("  $(food) = $(value(buy[food]))")
        end
    else
        println("The solver did not find an optimal solution.")
    end
end

function example_diet(; verbose = true)
    ## Nutrition guidelines
    categories = ["calories", "protein", "fat", "sodium"]
    category_data = Containers.DenseAxisArray([
        1800 2200;
        91   Inf;
        0    65;
        0    1779
        ], categories, ["min", "max"]
    )
    Test.@test category_data["protein", "min"] == 91.0
    Test.@test category_data["sodium", "max"] == 1779.0
    ## Foods
    foods = [
        "hamburger", "chicken", "hot dog", "fries", "macaroni", "pizza",
        "salad", "milk", "ice cream",
    ]
    cost = Containers.DenseAxisArray(
        [2.49, 2.89, 1.50, 1.89, 2.09, 1.99, 2.49, 0.89, 1.59],
        foods
    )
    food_data = Containers.DenseAxisArray(
        [
            410 24 26 730;
            420 32 10 1190;
            560 20 32 1800;
            380  4 19 270;
            320 12 10 930;
            320 15 12 820;
            320 31 12 1230;
            100  8 2.5 125;
            330  8 10 180
        ], foods, categories
    )
    Test.@test food_data["hamburger", "calories"] == 410.0
    Test.@test food_data["milk", "fat"] == 2.5
    ## Build model
    model = Model(GLPK.Optimizer)
    @variables(model, begin
        ## Variables for nutrition info
        category_data[c, "min"] <= nutrition[c = categories] <= category_data[c, "max"]
        ## Variables for which foods to buy
        buy[foods] >= 0
    end)
    ## Objective - minimize cost
    @objective(model, Min, sum(cost[f] * buy[f] for f in foods))
    ## Nutrition constraints
    @constraint(model, [c in categories],
        sum(food_data[f, c] * buy[f] for f in foods) == nutrition[c]
    )
    ## Solve
    if verbose
        println("Solving original problem...")
    end
    optimize!(model)
    term_status = termination_status(model)
    is_optimal = term_status == MOI.OPTIMAL
    Test.@test primal_status(model) == MOI.FEASIBLE_POINT
    Test.@test objective_value(model) ≈ 11.8288 atol = 1e-4
    if verbose
        print_solution(is_optimal, foods, buy)
    end
    ## Limit dairy (note that the problem will become infeasible).
    @constraint(model, buy["milk"] + buy["ice cream"] <= 6)
    if verbose
        println("Solving dairy-limited problem...")
    end
    optimize!(model)
    Test.@test termination_status(model) == MOI.INFEASIBLE
    Test.@test primal_status(model) == MOI.NO_SOLUTION
    if verbose
        print_solution(false, foods, buy)
    end
    return
end

example_diet()

#  Copyright 2017, Iain Dunning, Joey Huchette, Miles Lubin, and contributors
#  This Source Code Form is subject to the terms of the Mozilla Public
#  License, v. 2.0. If a copy of the MPL was not distributed with this
#  file, You can obtain one at http://mozilla.org/MPL/2.0/.
#############################################################################
# JuMP
# An algebraic modeling language for Julia
# See http://github.com/JuliaOpt/JuMP.jl
#############################################################################

using JuMP, GLPK, Test
const MOI = JuMP.MathOptInterface

"""
    example_diet()

Solve the classic "diet problem". Based on
http://www.gurobi.com/documentation/5.6/example-tour/diet_cpp_cpp
"""
function example_diet(; verbose = true)
    function print_solution(is_optimal, foods, buy)
        println("RESULTS:")
        if is_optimal
            for food in foods
                println("  $(food) = $(JuMP.value(buy[food]))")
            end
        else
            println("The solver did not find an optimal solution.")
        end
    end

    # Nutrition guidelines
    categories = ["calories", "protein", "fat", "sodium"]
    category_data = JuMP.Containers.DenseAxisArray([
        1800 2200;
        91   Inf;
        0    65;
        0    1779
        ], categories, ["min", "max"]
    )
    @test category_data["protein", "min"] == 91.0
    @test category_data["sodium", "max"] == 1779.0

    # Foods
    foods = ["hamburger", "chicken", "hot dog", "fries", "macaroni", "pizza",
             "salad", "milk", "ice cream"]
    cost = JuMP.Containers.DenseAxisArray(
        [2.49, 2.89, 1.50, 1.89, 2.09, 1.99, 2.49, 0.89, 1.59],
        foods
    )
    food_data = JuMP.Containers.DenseAxisArray(
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
    @test food_data["hamburger", "calories"] == 410.0
    @test food_data["milk", "fat"] == 2.5

    # Build model
    model = Model(GLPK.Optimizer)

    @variables(model, begin
        # Variables for nutrition info
        category_data[c, "min"] <= nutrition[c = categories] <= category_data[c, "max"]
        # Variables for which foods to buy
        buy[foods] >= 0
    end)

    # Objective - minimize cost
    @objective(model, Min, sum(cost[f] * buy[f] for f in foods))

    # Nutrition constraints
    @constraint(model, [c in categories],
        sum(food_data[f, c] * buy[f] for f in foods) == nutrition[c]
    )

    # Solve
    verbose && println("Solving original problem...")
    JuMP.optimize!(model)
    term_status = JuMP.termination_status(model)
    is_optimal = term_status == MOI.OPTIMAL
    @test JuMP.primal_status(model) == MOI.FEASIBLE_POINT
    @test JuMP.objective_value(model) ≈ 11.8288 atol = 1e-4
    verbose && print_solution(is_optimal, foods, buy)

    # Limit dairy (note that the problem will become infeasible).
    @constraint(model, buy["milk"] + buy["ice cream"] <= 6)
    verbose && println("Solving dairy-limited problem...")
    JuMP.optimize!(model)
    @test JuMP.termination_status(model) == MOI.INFEASIBLE
    @test JuMP.primal_status(model) == MOI.NO_SOLUTION
    verbose && print_solution(false, foods, buy)
end

example_diet(verbose = false)

#  Copyright 2017, Iain Dunning, Joey Huchette, Miles Lubin, and contributors
#  This Source Code Form is subject to the terms of the Mozilla Public
#  License, v. 2.0. If a copy of the MPL was not distributed with this
#  file, You can obtain one at http://mozilla.org/MPL/2.0/.
#############################################################################
# JuMP
# An algebraic modeling langauge for Julia
# See http://github.com/JuliaOpt/JuMP.jl
#############################################################################
# diet.jl
#
# Solve the classic "diet problem".
# Based on
#  http://www.gurobi.com/documentation/5.6/example-tour/diet_cpp_cpp
#############################################################################

using JuMP, Clp

solver = ClpSolver()

function PrintSolution(status, foods, buy)
    println("RESULTS:")
    if status == :Optimal
        for i = 1:length(foods)
            println("  $(foods[i]) = $(getvalue(buy[i]))")
        end
    else
        println("  No solution")
    end
    println("")
end

function SolveDiet()

    # Nutrition guidelines
    numCategories = 4
    categories = ["calories", "protein", "fat", "sodium"]
    minNutrition = [1800, 91, 0, 0]
    maxNutrition = [2200, Inf, 65, 1779]

    # Foods
    numFoods = 9
    foods = ["hamburger", "chicken", "hot dog", "fries",
                     "macaroni", "pizza", "salad", "milk", "ice cream"]
    cost = [2.49, 2.89, 1.50, 1.89, 2.09, 1.99, 2.49, 0.89, 1.59]
    nutritionValues = [410 24 26 730;
                       420 32 10 1190;
                       560 20 32 1800;
                       380  4 19 270;
                       320 12 10 930;
                       320 15 12 820;
                       320 31 12 1230;
                       100  8 2.5 125;
                       330  8 10 180]

    # Build model
    m = Model(solver=solver)

    # Variables for nutrition info
    @variable(m, minNutrition[i] <= nutrition[i=1:numCategories] <= maxNutrition[i])
    # Variables for which foods to buy
    @variable(m, buy[i=1:numFoods] >= 0)

    # Objective - minimize cost
    @objective(m, Min, dot(cost, buy))

    # Nutrition constraints
    for j = 1:numCategories
        @constraint(m, sum(nutritionValues[i,j]*buy[i] for i=1:numFoods) == nutrition[j])
    end

    # Solve
    println("Solving original problem...")
    status = solve(m)
    PrintSolution(status, foods, buy)

    # Limit dairy
    @constraint(m, buy[8] + buy[9] <= 6)
    println("Solving dairy-limited problem...")
    status = solve(m)
    PrintSolution(status, foods, buy)

end

SolveDiet()


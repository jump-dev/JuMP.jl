# Copyright 2017, Iain Dunning, Joey Huchette, Miles Lubin, and contributors    #src
# This Source Code Form is subject to the terms of the Mozilla Public License   #src
# v.2.0. If a copy of the MPL was not distributed with this file, You can       #src
# obtain one at https://mozilla.org/MPL/2.0/.                                   #src

# # The knapsack problem

# Formulate and solve a simple knapsack problem:
#
#     max sum(p_j x_j)
#      st sum(w_j x_j) <= C
#         x binary

using JuMP
import GLPK
import Test

function example_knapsack(; verbose = true)
    profit = [5, 3, 2, 7, 4]
    weight = [2, 8, 4, 2, 5]
    capacity = 10
    model = Model(GLPK.Optimizer)
    @variable(model, x[1:5], Bin)
    ## Objective: maximize profit
    @objective(model, Max, profit' * x)
    ## Constraint: can carry all
    @constraint(model, weight' * x <= capacity)
    ## Solve problem using MIP solver
    optimize!(model)
    if verbose
        println("Objective is: ", objective_value(model))
        println("Solution is:")
        for i in 1:5
            print("x[$i] = ", value(x[i]))
            println(", p[$i]/w[$i] = ", profit[i] / weight[i])
        end
    end
    Test.@test termination_status(model) == MOI.OPTIMAL
    Test.@test primal_status(model) == MOI.FEASIBLE_POINT
    Test.@test objective_value(model) == 16.0
    return
end

example_knapsack()

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
    example_knapsack(; verbose = true)

Formulate and solve a simple knapsack problem:
    max sum(p_j x_j)
     st sum(w_j x_j) <= C
        x binary
"""
function example_knapsack(; verbose = true)
    profit = [5, 3, 2, 7, 4]
    weight = [2, 8, 4, 2, 5]
    capacity = 10
    model = Model(GLPK.Optimizer)
    @variable(model, x[1:5], Bin)
    # Objective: maximize profit
    @objective(model, Max, profit' * x)
    # Constraint: can carry all
    @constraint(model, weight' * x <= capacity)
    # Solve problem using MIP solver
    JuMP.optimize!(model)
    if verbose
        println("Objective is: ", JuMP.objective_value(model))
        println("Solution is:")
        for i in 1:5
            print("x[$i] = ", JuMP.value(x[i]))
            println(", p[$i]/w[$i] = ", profit[i] / weight[i])
        end
    end
    @test JuMP.termination_status(model) == MOI.OPTIMAL
    @test JuMP.primal_status(model) == MOI.FEASIBLE_POINT
    @test JuMP.objective_value(model) == 16.0
end

example_knapsack(verbose = false)

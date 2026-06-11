# Copyright 2017, Iain Dunning, Joey Huchette, Miles Lubin, and contributors    #src
# This Source Code Form is subject to the terms of the Mozilla Public License   #src
# v.2.0. If a copy of the MPL was not distributed with this file, You can       #src
# obtain one at https://mozilla.org/MPL/2.0/.                                   #src

# # Quadratically constrained programs

# A simple quadratically constrained program based on an [example from Gurobi](https://www.gurobi.com/documentation/9.0/examples/qcp_c_c.html).

using JuMP
import Ipopt
import Test

function example_qcp(; verbose = true)
    model = Model(Ipopt.Optimizer)
    set_silent(model)
    @variable(model, x)
    @variable(model, y >= 0)
    @variable(model, z >= 0)
    @objective(model, Max, x)
    @constraint(model, x + y + z == 1)
    @constraint(model, x * x + y * y - z * z <= 0)
    @constraint(model, x * x - y * z <= 0)
    optimize!(model)
    if verbose
        print(model)
        println("Objective value: ", objective_value(model))
        println("x = ", value(x))
        println("y = ", value(y))
    end
    Test.@test termination_status(model) == LOCALLY_SOLVED
    Test.@test primal_status(model) == FEASIBLE_POINT
    Test.@test objective_value(model) ≈ 0.32699 atol = 1e-5
    Test.@test value(x) ≈ 0.32699 atol = 1e-5
    Test.@test value(y) ≈ 0.25707 atol = 1e-5
    return
end

example_qcp()

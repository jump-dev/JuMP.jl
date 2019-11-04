#  Copyright 2017, Iain Dunning, Joey Huchette, Miles Lubin, and contributors
#  This Source Code Form is subject to the terms of the Mozilla Public
#  License, v. 2.0. If a copy of the MPL was not distributed with this
#  file, You can obtain one at http://mozilla.org/MPL/2.0/.
#############################################################################
# JuMP
# An algebraic modeling language for Julia
# See http://github.com/JuliaOpt/JuMP.jl
#############################################################################

using JuMP, Ipopt, Test
const MOI = JuMP.MathOptInterface

"""
    example_qcp(; verbose = true)

A simple quadratically constrained program based on
http://www.gurobi.com/documentation/5.5/example-tour/node25
"""
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
    JuMP.optimize!(model)
    if verbose
        print(model)
        println("Objective value: ", JuMP.objective_value(model))
        println("x = ", JuMP.value(x))
        println("y = ", JuMP.value(y))
    end
    @test JuMP.termination_status(model) == MOI.LOCALLY_SOLVED
    @test JuMP.primal_status(model) == MOI.FEASIBLE_POINT
    @test JuMP.objective_value(model) ≈ 0.32699 atol = 1e-5
    @test JuMP.value(x) ≈ 0.32699 atol = 1e-5
    @test JuMP.value(y) ≈ 0.25707 atol = 1e-5
end

example_qcp(verbose = false)

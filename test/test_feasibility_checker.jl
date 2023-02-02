#  Copyright 2017, Iain Dunning, Joey Huchette, Miles Lubin, and contributors
#  This Source Code Form is subject to the terms of the Mozilla Public
#  License, v. 2.0. If a copy of the MPL was not distributed with this
#  file, You can obtain one at https://mozilla.org/MPL/2.0/.

module TestFeasibilityChecker

using JuMP
using Test

function test_no_solution()
    model = Model()
    @variable(model, x, Bin)
    @test_throws ErrorException primal_feasibility_report(model)
end

function test_primal_solution()
    model = Model(() -> MOIU.MockOptimizer(MOIU.Model{Float64}()))
    @variable(model, x, Bin)
    optimize!(model)
    mock = unsafe_backend(model)
    MOI.set(mock, MOI.TerminationStatus(), MOI.OPTIMAL)
    MOI.set(mock, MOI.PrimalStatus(), MOI.FEASIBLE_POINT)
    MOI.set(mock, MOI.VariablePrimal(), optimizer_index(x), 1.0)
    report = primal_feasibility_report(model)
    @test isempty(report)
end

function test_primal_solution_func()
    model = Model(() -> MOIU.MockOptimizer(MOIU.Model{Float64}()))
    @variable(model, x, Bin)
    optimize!(model)
    mock = unsafe_backend(model)
    MOI.set(mock, MOI.TerminationStatus(), MOI.OPTIMAL)
    MOI.set(mock, MOI.PrimalStatus(), MOI.FEASIBLE_POINT)
    MOI.set(mock, MOI.VariablePrimal(), optimizer_index(x), 1.0)
    report = primal_feasibility_report(model) do xi
        return value(xi)
    end
    @test isempty(report)
    return
end

function test_feasible()
    model = Model()
    @variable(model, x, Bin)
    @variable(model, 0 <= y <= 2, Int)
    @variable(model, z == 0.5)
    @constraint(model, x + y + z >= 0.5)
    report =
        primal_feasibility_report(model, Dict(x => 0.0, y => 0.0, z => 0.5))
    @test isempty(report)
end

function test_missing()
    model = Model()
    @variable(model, x, Bin)
    @variable(model, 0 <= y <= 2, Int)
    @variable(model, z == 0.5)
    @constraint(model, x + y + z >= 0.5)
    report =
        primal_feasibility_report(model, Dict(z => 0.0); skip_missing = true)
    @test report[FixRef(z)] == 0.5
    @test length(report) == 1
end

function test_missing_error()
    model = Model()
    @variable(model, x, Bin)
    @variable(model, 0 <= y <= 2, Int)
    point = Dict(x => 0.1)
    err = ErrorException(
        "point does not contain a value for variable $(y). Provide a value, " *
        "or pass `skip_missing = true`.",
    )
    @test_throws err primal_feasibility_report(model, point)
end

function test_bounds()
    model = Model()
    @variable(model, x, Bin)
    @variable(model, 0 <= y <= 2, Int)
    @variable(model, z == 0.5)
    point = Dict(x => 0.1, y => 2.1, z => 0.0)
    report = primal_feasibility_report(model, point)
    @test report[BinaryRef(x)] ≈ 0.1
    @test report[UpperBoundRef(y)] ≈ 0.1
    @test report[IntegerRef(y)] ≈ 0.1
    @test report[FixRef(z)] ≈ 0.5
    @test length(report) == 4
end

function test_scalar_affine()
    model = Model()
    @variable(model, x)
    @constraint(model, c1, x <= 0.5)
    @constraint(model, c2, x >= 1.25)
    @constraint(model, c3, x == 1.1)
    @constraint(model, c4, 0 <= x <= 0.5)
    report = primal_feasibility_report(model, Dict(x => 1.0))
    @test report[c1] ≈ 0.5
    @test report[c2] ≈ 0.25
    @test report[c3] ≈ 0.1
    @test report[c4] ≈ 0.5
    @test length(report) == 4
end

function test_scalar_affine_func()
    model = Model()
    @variable(model, x)
    @constraint(model, c1, x <= 0.5)
    @constraint(model, c2, x >= 1.25)
    @constraint(model, c3, x == 1.1)
    @constraint(model, c4, 0 <= x <= 0.5)
    report = primal_feasibility_report(model) do _
        return 1.0
    end
    @test report[c1] ≈ 0.5
    @test report[c2] ≈ 0.25
    @test report[c3] ≈ 0.1
    @test report[c4] ≈ 0.5
    @test length(report) == 4
    return
end

function test_scalar_quadratic()
    model = Model()
    @variable(model, x)
    @constraint(model, c1, x^2 + x <= 0.5)
    @constraint(model, c2, x^2 - x >= 1.25)
    @constraint(model, c3, x^2 + x == 1.1)
    @constraint(model, c4, 0 <= x^2 + x <= 0.5)
    report = primal_feasibility_report(model, Dict(x => 1.0))
    @test report[c1] ≈ 1.5
    @test report[c2] ≈ 1.25
    @test report[c3] ≈ 0.9
    @test report[c4] ≈ 1.5
    @test length(report) == 4
end

function test_vector()
    model = Model()
    @variable(model, x[1:3])
    @constraint(model, c1, x in MOI.Nonnegatives(3))
    @constraint(model, c2, x in MOI.Nonpositives(3))
    @constraint(model, c3, x in MOI.Reals(3))
    @constraint(model, c4, x in MOI.Zeros(3))
    point = Dict(x[1] => 1.0, x[2] => -1.0, x[3] => 0.0)
    report = primal_feasibility_report(model, point)
    @test report[c1] ≈ 1.0
    @test report[c2] ≈ 1.0
    @test !haskey(report, c3)
    @test report[c4] ≈ sqrt(2)
    @test length(report) == 3
end

function test_vector_affine()
    model = Model()
    @variable(model, x[1:3])
    @constraint(model, c1, 2 * x in MOI.Nonnegatives(3))
    @constraint(model, c2, 2 * x in MOI.Nonpositives(3))
    @constraint(model, c3, 2 * x in MOI.Reals(3))
    @constraint(model, c4, 2 * x in MOI.Zeros(3))
    point = Dict(x[1] => 1.0, x[2] => -1.0, x[3] => 0.0)
    report = primal_feasibility_report(model, point)
    @test report[c1] ≈ 2.0
    @test report[c2] ≈ 2.0
    @test !haskey(report, c3)
    @test report[c4] ≈ sqrt(8)
    @test length(report) == 3
end

function test_nonlinear()
    model = Model()
    @variable(model, x)
    @NLconstraint(model, c1, sin(x) <= 0.0)
    @NLconstraint(model, c2, sin(x) <= 1.0)
    @NLconstraint(model, c3, exp(x) >= 2.0)
    @NLconstraint(model, c4, x + x^2 - x^3 == 0.5)
    report = primal_feasibility_report(model, Dict(x => 0.5))
    @test report[c1] ≈ sin(0.5)
    @test !haskey(report, c2)
    @test report[c3] ≈ 2 - exp(0.5)
    @test report[c4] ≈ 0.125
end

function test_nonlinear_missing()
    model = Model()
    @variable(model, x)
    @NLconstraint(model, c1, sin(x) <= 0.0)
    @test_throws(
        ErrorException(
            "`skip_missing = true` is not allowed when nonlinear constraints " *
            "are present.",
        ),
        primal_feasibility_report(model, Dict(x => 0.5); skip_missing = true)
    )
end

end

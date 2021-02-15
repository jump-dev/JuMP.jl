#  Copyright 2017, Iain Dunning, Joey Huchette, Miles Lubin, and contributors
#  This Source Code Form is subject to the terms of the Mozilla Public
#  License, v. 2.0. If a copy of the MPL was not distributed with this
#  file, You can obtain one at https://mozilla.org/MPL/2.0/.

module TestFeasibilityChecker

using JuMP
using Test

function test_unsupported()
    @test_throws(
        ErrorException,
        JuMP._distance_to_set([1.0, 1.0], MOI.Complements(1)),
    )
end

function test_lessthan()
    @test JuMP._distance_to_set(1.0, MOI.LessThan(2.0)) ≈ 0.0
    @test JuMP._distance_to_set(1.0, MOI.LessThan(0.5)) ≈ 0.5
end

function test_greaterthan()
    @test JuMP._distance_to_set(1.0, MOI.GreaterThan(2.0)) ≈ 1.0
    @test JuMP._distance_to_set(1.0, MOI.GreaterThan(0.5)) ≈ 0.0
end

function test_equalto()
    @test JuMP._distance_to_set(1.0, MOI.EqualTo(2.0)) ≈ 1.0
    @test JuMP._distance_to_set(1.0, MOI.EqualTo(0.5)) ≈ 0.5
end

function test_interval()
    @test JuMP._distance_to_set(1.0, MOI.Interval(1.0, 2.0)) ≈ 0.0
    @test JuMP._distance_to_set(0.5, MOI.Interval(1.0, 2.0)) ≈ 0.5
    @test JuMP._distance_to_set(2.75, MOI.Interval(1.0, 2.0)) ≈ 0.75
end

function test_zeroone()
    @test JuMP._distance_to_set(0.6, MOI.ZeroOne()) ≈ 0.4
    @test JuMP._distance_to_set(-0.01, MOI.ZeroOne()) ≈ 0.01
    @test JuMP._distance_to_set(1.01, MOI.ZeroOne()) ≈ 0.01
end

function test_integer()
    @test JuMP._distance_to_set(0.6, MOI.Integer()) ≈ 0.4
    @test JuMP._distance_to_set(3.1, MOI.Integer()) ≈ 0.1
    @test JuMP._distance_to_set(-0.01, MOI.Integer()) ≈ 0.01
    @test JuMP._distance_to_set(1.01, MOI.Integer()) ≈ 0.01
end

function test_semicontinuous()
    s = MOI.Semicontinuous(2.0, 4.0)
    @test JuMP._distance_to_set(-2.0, s) ≈ 2.0
    @test JuMP._distance_to_set(0.5, s) ≈ 0.5
    @test JuMP._distance_to_set(1.9, s) ≈ 0.1
    @test JuMP._distance_to_set(2.1, s) ≈ 0.0
    @test JuMP._distance_to_set(4.1, s) ≈ 0.1
end

function test_semiintger()
    s = MOI.Semiinteger(1.9, 4.0)
    @test JuMP._distance_to_set(-2.0, s) ≈ 2.0
    @test JuMP._distance_to_set(0.5, s) ≈ 0.5
    @test JuMP._distance_to_set(1.9, s) ≈ 0.1
    @test JuMP._distance_to_set(2.1, s) ≈ 0.1
    @test JuMP._distance_to_set(4.1, s) ≈ 0.1
end

function test_nonnegatives()
    @test JuMP._distance_to_set([-1.0, 1.0], MOI.Nonnegatives(2)) ≈ 1.0
end

function test_nonpositives()
    @test JuMP._distance_to_set([-1.0, 1.0], MOI.Nonpositives(2)) ≈ 1.0
end

function test_reals()
    @test JuMP._distance_to_set([-1.0, 1.0], MOI.Reals(2)) ≈ 0.0
end

function test_zeros()
    @test JuMP._distance_to_set([-1.0, 1.0], MOI.Zeros(2)) ≈ sqrt(2)
end

function test_feasible()
    model = Model()
    @variable(model, x, Bin)
    @variable(model, 0 <= y <= 2, Int)
    @variable(model, z == 0.5)
    @constraint(model, x + y + z >= 0.5)
    @test primal_feasibility_report(model, Dict(z => 0.5)) === nothing
end

function test_bounds()
    model = Model()
    @variable(model, x, Bin)
    @variable(model, 0 <= y <= 2, Int)
    @variable(model, z == 0.5)
    report = primal_feasibility_report(model, Dict(x => 0.1, y => 2.1))
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

function test_vector_affine()
    model = Model()
    @variable(model, x[1:3])
    @constraint(model, c1, x in MOI.Nonnegatives(3))
    @constraint(model, c2, x in MOI.Nonpositives(3))
    @constraint(model, c3, x in MOI.Reals(3))
    @constraint(model, c4, x in MOI.Zeros(3))
    report = primal_feasibility_report(model, Dict(x[1] => 1.0, x[2] => -1.0))
    @test report[c1] ≈ 1.0
    @test report[c2] ≈ 1.0
    @test !haskey(report, c3)
    @test report[c4] ≈ sqrt(2)
    @test length(report) == 3
end

function runtests()
    for name in names(@__MODULE__; all = true)
        if !startswith("$(name)", "test_")
            continue
        end
        @testset "$(name)" begin
            getfield(@__MODULE__, name)()
        end
    end
end

end

TestFeasibilityChecker.runtests()

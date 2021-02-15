#  Copyright 2017, Iain Dunning, Joey Huchette, Miles Lubin, and contributors
#  This Source Code Form is subject to the terms of the Mozilla Public
#  License, v. 2.0. If a copy of the MPL was not distributed with this
#  file, You can obtain one at https://mozilla.org/MPL/2.0/.

module TestFeasibilityChecker

using JuMP
using Test

function test_distance_to_set()
    @test JuMP._distance_to_set(1.0, MOI.LessThan(2.0)) ≈ 0.0
    @test JuMP._distance_to_set(1.0, MOI.LessThan(0.5)) ≈ 0.5
    @test JuMP._distance_to_set(1.0, MOI.GreaterThan(2.0)) ≈ 1.0
    @test JuMP._distance_to_set(1.0, MOI.GreaterThan(0.5)) ≈ 0.0
    @test JuMP._distance_to_set(1.0, MOI.EqualTo(2.0)) ≈ 1.0
    @test JuMP._distance_to_set(1.0, MOI.EqualTo(0.5)) ≈ 0.5
    @test JuMP._distance_to_set(1.0, MOI.Interval(1.0, 2.0)) ≈ 0.0
    @test JuMP._distance_to_set(0.5, MOI.Interval(1.0, 2.0)) ≈ 0.5
    @test JuMP._distance_to_set(2.75, MOI.Interval(1.0, 2.0)) ≈ 0.75
    @test JuMP._distance_to_set(0.6, MOI.ZeroOne()) ≈ 0.4
    @test JuMP._distance_to_set(-0.01, MOI.ZeroOne()) ≈ 0.01
    @test JuMP._distance_to_set(1.01, MOI.ZeroOne()) ≈ 0.01
    @test JuMP._distance_to_set(0.6, MOI.Integer()) ≈ 0.4
    @test JuMP._distance_to_set(3.1, MOI.Integer()) ≈ 0.1
    @test JuMP._distance_to_set(-0.01, MOI.Integer()) ≈ 0.01
    @test JuMP._distance_to_set(1.01, MOI.Integer()) ≈ 0.01
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

function test_affine()
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

#  Copyright 2017, Iain Dunning, Joey Huchette, Miles Lubin, and contributors
#  This Source Code Form is subject to the terms of the Mozilla Public
#  License, v. 2.0. If a copy of the MPL was not distributed with this
#  file, You can obtain one at https://mozilla.org/MPL/2.0/.

using JuMP
using Test

struct DummyCallbackData end

@testset "LazyConstraint" begin
    mock = MOI.Utilities.MockOptimizer(
        MOI.Utilities.UniversalFallback(MOI.Utilities.Model{Float64}()),
    )
    model = Model(() -> mock)
    @variable(model, 0 <= x <= 2.5, Int)
    # We attach optimizer here, because when submitting the lazy constraint
    # in reality, we will be in `optimize!`, so the optimizer will be
    # attached.
    MOIU.attach_optimizer(model)
    con = @build_constraint(x <= 2)
    MOI.submit(model, MOI.LazyConstraint(DummyCallbackData()), con)
    @test length(mock.submitted) == 1
    c = mock.submitted[MOI.LazyConstraint(DummyCallbackData())]
    @test length(c) == 1
    @test c[1][1] isa MOI.ScalarAffineFunction{Float64}
    @test c[1][1].constant == 0.0
    @test length(c[1][1].terms) == 1
    @test c[1][1].terms[1].coefficient == 1.0
    @test c[1][2] == MOI.LessThan(2.0)
end

@testset "UserCut" begin
    mock = MOI.Utilities.MockOptimizer(
        MOI.Utilities.UniversalFallback(MOI.Utilities.Model{Float64}()),
    )
    model = direct_model(mock)
    @variable(model, 0 <= x <= 2.5, Int)
    con = @build_constraint(x <= 2)
    MOI.submit(model, MOI.UserCut(DummyCallbackData()), con)
    @test length(mock.submitted) == 1
    c = mock.submitted[MOI.UserCut(DummyCallbackData())]
    @test length(c) == 1
    @test c[1][1] ≈ moi_function(1.0 * x)
    @test c[1][2] == MOI.LessThan(2.0)
end

@testset "HeuristicSolution" begin
    mock = MOI.Utilities.MockOptimizer(
        MOI.Utilities.UniversalFallback(MOI.Utilities.Model{Float64}()),
    )
    model = direct_model(mock)
    @variable(model, 0 <= x <= 2.5, Int)
    con = @build_constraint(x <= 2)
    MOI.submit(model, MOI.HeuristicSolution(DummyCallbackData()), [x], [0.0])
    @test length(mock.submitted) == 1
    c = mock.submitted[MOI.HeuristicSolution(DummyCallbackData())]
    @test length(c) == 1
    @test c[1][1] == [index(x)]
    @test c[1][2] == [0.0]
    # Non-Float64 case
    MOI.submit(model, MOI.HeuristicSolution(DummyCallbackData()), [x], [1])
    @test length(mock.submitted) == 1
    c = mock.submitted[MOI.HeuristicSolution(DummyCallbackData())]
    @test length(c) == 2
    @test c[2][1] == [index(x)]
    @test c[2][2] == [1.0]
end

@testset "callback_value" begin
    mock = MOI.Utilities.MockOptimizer(
        MOI.Utilities.UniversalFallback(MOI.Utilities.Model{Float64}()),
    )
    cb = DummyCallbackData()
    model = direct_model(mock)
    @variable(model, 0 <= x <= 2.5, Int)
    MOIU.set_mock_optimize!(
        mock,
        mock -> begin
            MOI.set(mock, MOI.TerminationStatus(), MOI.OPTIMAL)
            MOI.set(mock, MOI.CallbackVariablePrimal(cb), index(x), 1)
        end,
    )
    optimize!(model)
    @test callback_value(cb, x) == 1
    expr = x + 1
    @test callback_value(cb, expr) == 2
    quad_expr = expr^2
    @test callback_value(cb, quad_expr) == 4
end

@testset "callback_node_status" begin
    mock = MOI.Utilities.MockOptimizer(
        MOI.Utilities.UniversalFallback(MOI.Utilities.Model{Float64}()),
    )
    cb = DummyCallbackData()
    model = direct_model(mock)
    @variable(model, 0 <= x <= 2.5, Int)
    MOIU.set_mock_optimize!(
        mock,
        mock -> begin
            MOI.set(mock, MOI.TerminationStatus(), MOI.OPTIMAL)
            MOI.set(
                mock,
                MOI.CallbackNodeStatus(cb),
                MOI.CALLBACK_NODE_STATUS_INTEGER,
            )
        end,
    )
    optimize!(model)
    @test callback_node_status(cb, model) == MOI.CALLBACK_NODE_STATUS_INTEGER
end

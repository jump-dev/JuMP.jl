#  Copyright 2017, Iain Dunning, Joey Huchette, Miles Lubin, and contributors
#  This Source Code Form is subject to the terms of the Mozilla Public
#  License, v. 2.0. If a copy of the MPL was not distributed with this
#  file, You can obtain one at http://mozilla.org/MPL/2.0/.
#############################################################################
# JuMP
# An algebraic modeling langauge for Julia
# See http://github.com/JuliaOpt/JuMP.jl
#############################################################################
# test/model.jl
#############################################################################

using JuMP

using Compat
using Compat.Test

using MathOptInterface
const MOI = MathOptInterface
const MOIT = MOI.Test
const MOIU = MOI.Utilities

# Simple LP model not supporting Interval
@MOIU.model SimpleLPModel () (EqualTo, GreaterThan, LessThan) () () (SingleVariable,) (ScalarAffineFunction,) () ()

struct Optimizer
    a::Int
    b::Int
end
function opt_build(a::Int; b::Int = 1)
    return Optimizer(a, b)
end

function test_model()
    @testset "optimize_hook" begin
        m = Model()
        @test m.optimize_hook === nothing
        called = false
        function hook(m)
            called = true
        end
        JuMP.set_optimize_hook(m, hook)
        @test !called
        optimize!(m)
        @test called
    end

    @testset "UniversalFallback" begin
        m = Model()
        MOI.set!(m, MOIT.UnknownModelAttribute(), 1)
        @test MOI.get(m, MOIT.UnknownModelAttribute()) == 1
    end

    @testset "Bridges" begin
        @testset "Automatic bridging" begin
            # optimizer not supporting Interval
            model = Model(with_optimizer(MOIU.MockOptimizer,
                                         SimpleLPModel{Float64}()))
            @variable model x
            cref = @constraint model 0 <= x + 1 <= 1
            @test cref isa JuMP.ConstraintRef{JuMP.Model,MOI.ConstraintIndex{MOI.ScalarAffineFunction{Float64},MOI.Interval{Float64}}}
            JuMP.optimize!(model)
        end
        @testset "Automatic bridging disabled with `bridge_constraints` keyword" begin
            model = Model(with_optimizer(MOIU.MockOptimizer,
                                         SimpleLPModel{Float64}()),
                          bridge_constraints=false)
            @test model.moi_backend isa MOIU.CachingOptimizer
            @test model.moi_backend === JuMP.caching_optimizer(model)
            @variable model x
            @test_throws ErrorException @constraint model 0 <= x + 1 <= 1
        end
        @testset "No bridge automatically added in Direct mode" begin
            optimizer = MOIU.MockOptimizer(SimpleLPModel{Float64}())
            model = JuMP.direct_model(optimizer)
            @variable model x
            @test_throws ErrorException @constraint model 0 <= x + 1 <= 1
        end
    end

    @testset "Factories" begin
        factory = with_optimizer(Optimizer, 1, 2)
        @test factory.constructor == Optimizer
        @test factory.args == (1, 2)
        optimizer = factory()
        @test optimizer isa Optimizer
        @test optimizer.a == 1
        @test optimizer.b == 2
        @test_throws ErrorException factory = with_optimizer(opt_build, 1, 2)
        factory = with_optimizer(opt_build, 1, b = 2)
        @test factory.constructor == opt_build
        @test factory.args == (1,)
        optimizer = factory()
        @test optimizer isa Optimizer
        @test optimizer.a == 1
        @test optimizer.b == 2
    end
end

@testset "Model" begin
    test_model()
end

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

using Test

using MathOptInterface
const MOI = MathOptInterface
const MOIT = MOI.Test
const MOIU = MOI.Utilities

# Simple LP model not supporting Interval
@MOIU.model(SimpleLPModel,
            (),
            (MOI.EqualTo, MOI.GreaterThan, MOI.LessThan),
            (),
            (),
            (MOI.SingleVariable,),
            (MOI.ScalarAffineFunction,),
            (),
            ())

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
        MOI.set(m, MOIT.UnknownModelAttribute(), 1)
        @test MOI.get(m, MOIT.UnknownModelAttribute()) == 1
    end

    @testset "Bridges" begin
        @testset "Automatic bridging" begin
            # optimizer not supporting Interval
            model = Model(with_optimizer(MOIU.MockOptimizer,
                                         SimpleLPModel{Float64}()))
            @test JuMP.bridge_constraints(model)
            @test JuMP.backend(model) isa MOIU.CachingOptimizer
            @test JuMP.backend(model).optimizer isa MOI.Bridges.LazyBridgeOptimizer
            @test JuMP.backend(model).optimizer.model isa MOIU.MockOptimizer
            @variable model x
            cref = @constraint model 0 <= x + 1 <= 1
            @test cref isa JuMP.ConstraintRef{JuMP.Model,MOI.ConstraintIndex{MOI.ScalarAffineFunction{Float64},MOI.Interval{Float64}}}
            JuMP.optimize!(model)
        end
        @testset "Automatic bridging with cache for bridged model" begin
            # optimizer not supporting Interval and not supporting `default_copy_to`
            model = Model(with_optimizer(MOIU.MockOptimizer,
                                         SimpleLPModel{Float64}(),
                                         needs_allocate_load=true))
            @test JuMP.bridge_constraints(model)
            @test JuMP.backend(model) isa MOIU.CachingOptimizer
            @test JuMP.backend(model).optimizer isa MOI.Bridges.LazyBridgeOptimizer
            @test JuMP.backend(model).optimizer.model isa MOIU.CachingOptimizer
            @test JuMP.backend(model).optimizer.model.optimizer isa MOIU.MockOptimizer
            @variable model x
            cref = @constraint model 0 <= x + 1 <= 1
            @test cref isa JuMP.ConstraintRef{JuMP.Model,MOI.ConstraintIndex{MOI.ScalarAffineFunction{Float64},MOI.Interval{Float64}}}
            JuMP.optimize!(model)
        end
        @testset "Automatic bridging disabled with `bridge_constraints` keyword" begin
            model = Model(with_optimizer(MOIU.MockOptimizer,
                                         SimpleLPModel{Float64}()),
                          bridge_constraints=false)
            @test !JuMP.bridge_constraints(model)
            @test JuMP.backend(model) isa MOIU.CachingOptimizer
            @test !(JuMP.backend(model).optimizer isa MOI.Bridges.LazyBridgeOptimizer)
            @variable model x
            err = ErrorException("Constraints of type MathOptInterface.ScalarAffineFunction{Float64}-in-MathOptInterface.Interval{Float64} are not supported by the solver, try using `bridge_constraints=true` in the `JuMP.Model` constructor if you believe the constraint can be reformulated to constraints supported by the solver.")
            @test_throws err @constraint model 0 <= x + 1 <= 1
        end
        @testset "No bridge automatically added in Direct mode" begin
            optimizer = MOIU.MockOptimizer(SimpleLPModel{Float64}())
            model = JuMP.direct_model(optimizer)
            @test !JuMP.bridge_constraints(model)
            @variable model x
            err = ErrorException("Constraints of type MathOptInterface.ScalarAffineFunction{Float64}-in-MathOptInterface.Interval{Float64} are not supported by the solver.")
            @test_throws err @constraint model 0 <= x + 1 <= 1
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

    @testset "solver_name" begin
        @testset "Not attached" begin
            model = Model()
            @test JuMP.solver_name(model) == "No optimizer attached."
        end

        @testset "Mock" begin
            model = Model(with_optimizer(MOIU.MockOptimizer,
                                         SimpleLPModel{Float64}()))
            @test JuMP.solver_name(model) == "Mock"
        end
    end
end

@testset "Model" begin
    test_model()
end

struct DummyExtensionData
    model::JuMP.Model
end
function JuMP.copy_extension_data(data::DummyExtensionData,
                                  new_model::JuMP.AbstractModel,
                                  model::JuMP.AbstractModel)
    @test data.model === model
    return DummyExtensionData(new_model)
end
function dummy_optimizer_hook(::JuMP.AbstractModel) end

@testset "Model copy" begin
    for copy_model in (true, true)
        @testset "Using $(copy_model ? "JuMP.copy_model" : "Base.copy")" begin
            for caching_mode in (MOIU.Automatic, MOIU.Manual)
                @testset "In $caching_mode mode" begin
                    model = Model(caching_mode = caching_mode)
                    model.optimize_hook = dummy_optimizer_hook
                    data = DummyExtensionData(model)
                    model.ext[:dummy] = data
                    @variable(model, x ≥ 0, Bin)
                    @variable(model, y ≤ 1, Int)
                    @variable(model, z == 0)
                    @constraint(model, cref, x + y == 1)

                    if copy_model
                        new_model, reference_map = JuMP.copy_model(model)
                    else
                        new_model = copy(model)
                        reference_map = Dict{Union{JuMP.VariableRef,
                                                   JuMP.ConstraintRef},
                                             Union{JuMP.VariableRef,
                                                   JuMP.ConstraintRef}}()
                        reference_map[x] = new_model[:x]
                        reference_map[y] = new_model[:y]
                        reference_map[z] = new_model[:z]
                        reference_map[cref] = new_model[:cref]
                    end
                    @test MOIU.mode(JuMP.backend(new_model)) == caching_mode
                    @test new_model.optimize_hook === dummy_optimizer_hook
                    @test new_model.ext[:dummy].model === new_model
                    x_new = reference_map[x]
                    @test JuMP.owner_model(x_new) === new_model
                    @test JuMP.name(x_new) == "x"
                    y_new = reference_map[y]
                    @test JuMP.owner_model(y_new) === new_model
                    @test JuMP.name(y_new) == "y"
                    z_new = reference_map[z]
                    @test JuMP.owner_model(z_new) === new_model
                    @test JuMP.name(z_new) == "z"
                    if copy_model
                        @test JuMP.LowerBoundRef(x_new) == reference_map[JuMP.LowerBoundRef(x)]
                        @test JuMP.BinaryRef(x_new) == reference_map[JuMP.BinaryRef(x)]
                        @test JuMP.UpperBoundRef(y_new) == reference_map[JuMP.UpperBoundRef(y)]
                        @test JuMP.IntegerRef(y_new) == reference_map[JuMP.IntegerRef(y)]
                        @test JuMP.FixRef(z_new) == reference_map[JuMP.FixRef(z)]
                    end
                    cref_new = reference_map[cref]
                    @test cref_new.model === new_model
                    @test JuMP.name(cref_new) == "cref"
                end
            end
        end
    end
    @testset "In Direct mode" begin
        mock = MOIU.MockOptimizer(JuMP.JuMPMOIModel{Float64}())
        model = JuMP.direct_model(mock)
        @test_throws ErrorException JuMP.copy(model)
    end
end

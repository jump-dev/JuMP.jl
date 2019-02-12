#  Copyright 2017, Iain Dunning, Joey Huchette, Miles Lubin, and contributors
#  This Source Code Form is subject to the terms of the Mozilla Public
#  License, v. 2.0. If a copy of the MPL was not distributed with this
#  file, You can obtain one at http://mozilla.org/MPL/2.0/.
#############################################################################
# JuMP
# An algebraic modeling language for Julia
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

# Custom set Nonnegative with bridge NonnegativeBridge
include("nonnegative_bridge.jl")

function test_model()
    @testset "Result attributes" begin
        err = JuMP.OptimizeNotCalled()
        model = Model()
        @variable(model, x)
        c = @constraint(model, x ≤ 0)
        @objective(model, Max, x)
        @test_throws err JuMP.objective_value(model)
        @test_throws err JuMP.objective_bound(model)
        @test_throws err JuMP.value(x)
        @test_throws err JuMP.value(c)
        @test_throws err JuMP.dual(c)
    end

    @testset "Test variable/model 'hygiene'" begin
        model_x = Model()
        @variable(model_x, x)
        model_y = Model()
        @variable(model_y, y)
        err = JuMP.VariableNotOwned(y)
        @testset "Variable" begin
            @testset "constraint" begin
                @test_throws err @constraint(model_x, y in MOI.EqualTo(1.0))
                @test_throws err @constraint(model_x, [x, y] in MOI.Zeros(2))
            end
            @testset "objective" begin
                @test_throws err @objective(model_x, Min, y)
            end
        end
        @testset "Linear" begin
            @testset "constraint" begin
                @test_throws err @constraint(model_x, x + y == 1)
                @test_throws err begin
                    @constraint(model_x, [x, x + y] in MOI.Zeros(2))
                end
                @test_throws err begin
                    @constraint(model_x, [x + y, x] in MOI.Zeros(2))
                end
            end
            @testset "objective" begin
                @test_throws err @objective(model_x, Min, x + y)
            end
        end
        @testset "Quadratic" begin
            @testset "constraint" begin
                @test_throws err @constraint(model_x, x * y >= 1)
                @test_throws err begin
                    @constraint(model_x, [x, x * y] in MOI.Zeros(2))
                end
                @test_throws err begin
                    @constraint(model_x, [x * y, x] in MOI.Zeros(2))
                end
                @test_throws err @constraint(model_x, x * x + x + y <= 1)
                @test_throws err begin
                    @constraint(model_x, [x, x * x + x + y] in MOI.Zeros(2))
                end
                @test_throws err begin
                    @constraint(model_x, [x * x + x + y, x] in MOI.Zeros(2))
                end
            end
            @testset "objective" begin
                @test_throws err @objective(model_x, Min, x * y)
                @test_throws err @objective(model_x, Min, x * x + x + y)
            end
        end
        @testset "Attribute" begin
            cy = @constraint(model_y, y in MOI.EqualTo(1.0))
            cerr = JuMP.ConstraintNotOwned(cy)
            @testset "get" begin
                @test_throws err begin
                    MOI.get(model_x, MOI.VariablePrimalStart(), y)
                end
                @test_throws cerr begin
                    MOI.get(model_x, MOI.ConstraintPrimalStart(), cy)
                end
            end
            @testset "set" begin
                @test_throws err begin
                    MOI.set(model_x, MOI.VariablePrimalStart(), y, 1.0)
                end
                @test_throws cerr begin
                    MOI.set(model_x, MOI.ConstraintPrimalStart(), cy, 1.0)
                end
            end
        end
    end

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

        @testset "Add bridge" begin
            function mock()
                mock = MOIU.MockOptimizer(JuMP._MOIModel{Float64}(),
                                          eval_variable_constraint_dual=false)
                optimize!(mock) = MOIU.mock_optimize!(mock, [1.0],
                        (MOI.SingleVariable, MOI.GreaterThan{Float64}) => [2.0])
                MOIU.set_mock_optimize!(mock, optimize!)
                return mock
            end
            factory = with_optimizer(mock)
            @testset "before loading the constraint to the optimizer" begin
                @testset "with_optimizer at Model" begin
                    model = Model(factory)
                    @variable(model, x)
                    JuMP.add_bridge(model, NonnegativeBridge)
                    c = @constraint(model, x in Nonnegative())
                    JuMP.optimize!(model)
                    @test 1.0 == @inferred JuMP.value(x)
                    @test 1.0 == @inferred JuMP.value(c)
                    @test 2.0 == @inferred JuMP.dual(c)
                end
                @testset "with_optimizer at optimize!" begin
                    model = Model()
                    @variable(model, x)
                    c = @constraint(model, x in Nonnegative())
                    JuMP.add_bridge(model, NonnegativeBridge)
                    JuMP.optimize!(model, factory)
                    @test 1.0 == @inferred JuMP.value(x)
                    @test 1.0 == @inferred JuMP.value(c)
                    @test 2.0 == @inferred JuMP.dual(c)
                end
            end
            @testset "after loading the constraint to the optimizer" begin
                @testset "with_optimizer at Model" begin
                    err = ErrorException(string("Constraints of type ",
                    "MathOptInterface.SingleVariable-in-Nonnegative are not ",
                    "supported by the solver and there are no bridges that ",
                    "can reformulate it into supported constraints."))
                    model = Model(factory)
                    @variable(model, x)
                    @test_throws err @constraint(model, x in Nonnegative())
                    JuMP.add_bridge(model, NonnegativeBridge)
                    c = @constraint(model, x in Nonnegative())
                    JuMP.optimize!(model)
                    @test 1.0 == @inferred JuMP.value(x)
                    @test 1.0 == @inferred JuMP.value(c)
                    @test 2.0 == @inferred JuMP.dual(c)
                end
                @testset "with_optimizer at optimize!" begin
                    err = MOI.UnsupportedConstraint{MOI.SingleVariable,
                                                    Nonnegative}()
                    model = Model()
                    @variable(model, x)
                    c = @constraint(model, x in Nonnegative())
                    @test_throws err JuMP.optimize!(model, factory)
                    JuMP.add_bridge(model, NonnegativeBridge)
                    JuMP.optimize!(model)
                    @test 1.0 == @inferred JuMP.value(x)
                    @test 1.0 == @inferred JuMP.value(c)
                    @test 2.0 == @inferred JuMP.dual(c)
                end
            end
            @testset "automatically with BridgeableConstraint" begin
                @testset "with_optimizer at Model" begin
                    model = Model(factory)
                    @variable(model, x)
                    constraint = ScalarConstraint(x, Nonnegative())
                    bc = BridgeableConstraint(constraint, NonnegativeBridge)
                    c = add_constraint(model, bc)
                    JuMP.optimize!(model)
                    @test 1.0 == @inferred JuMP.value(x)
                    @test 1.0 == @inferred JuMP.value(c)
                    @test 2.0 == @inferred JuMP.dual(c)
                end
                @testset "with_optimizer at optimize!" begin
                    model = Model()
                    @variable(model, x)
                    constraint = ScalarConstraint(x, Nonnegative())
                    bc = BridgeableConstraint(constraint, NonnegativeBridge)
                    c = add_constraint(model, bc)
                    JuMP.optimize!(model, factory)
                    @test 1.0 == @inferred JuMP.value(x)
                    @test 1.0 == @inferred JuMP.value(c)
                    @test 2.0 == @inferred JuMP.dual(c)
                end
            end
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
            @test "No optimizer attached." == @inferred JuMP.solver_name(model)
        end

        @testset "Mock" begin
            model = Model(with_optimizer(MOIU.MockOptimizer,
                                         SimpleLPModel{Float64}()))
            @test "Mock" == @inferred JuMP.solver_name(model)
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
    for copy_model in (true, false)
        @testset "Using $(copy_model ? "JuMP.copy_model" : "Base.copy")" begin
            for caching_mode in (MOIU.AUTOMATIC, MOIU.MANUAL)
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
                    @test caching_mode == @inferred MOIU.mode(JuMP.backend(new_model))
                    @test new_model.optimize_hook === dummy_optimizer_hook
                    @test new_model.ext[:dummy].model === new_model
                    x_new = reference_map[x]
                    @test JuMP.owner_model(x_new) === new_model
                    @test "x" == @inferred JuMP.name(x_new)
                    y_new = reference_map[y]
                    @test JuMP.owner_model(y_new) === new_model
                    @test "y" == @inferred JuMP.name(y_new)
                    z_new = reference_map[z]
                    @test JuMP.owner_model(z_new) === new_model
                    @test "z" == @inferred JuMP.name(z_new)
                    if copy_model
                        @test reference_map[JuMP.LowerBoundRef(x)] == @inferred JuMP.LowerBoundRef(x_new)
                        @test reference_map[JuMP.BinaryRef(x)] == @inferred JuMP.BinaryRef(x_new)
                        @test reference_map[JuMP.UpperBoundRef(y)] == @inferred JuMP.UpperBoundRef(y_new)
                        @test reference_map[JuMP.IntegerRef(y)] == @inferred JuMP.IntegerRef(y_new)
                        @test reference_map[JuMP.FixRef(z)] == @inferred JuMP.FixRef(z_new)
                    end
                    cref_new = reference_map[cref]
                    @test cref_new.model === new_model
                    @test "cref" == @inferred JuMP.name(cref_new)
                end
            end
        end
    end
    @testset "In Direct mode" begin
        mock = MOIU.MockOptimizer(JuMP._MOIModel{Float64}())
        model = JuMP.direct_model(mock)
        @test_throws ErrorException JuMP.copy(model)
    end
end

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

# Simple LP model not supporting Interval
MOIU.@model(
    SimpleLPModel,
    (), (MOI.EqualTo, MOI.GreaterThan, MOI.LessThan), (), (),
    (), (MOI.ScalarAffineFunction,), (), ()
)

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
    @testset "NoOptimizer" begin
        err = NoOptimizer()
        model = Model()
        @variable(model, x)
        @test_throws err optimizer_index(x)
        cref = @constraint(model, x == 1)
        @test_throws err JuMP.optimizer_index(cref)
        @test_throws err JuMP.optimize!(model)
        @test_throws err JuMP.value(x)
        @test_throws err JuMP.value(cref)
        @test_throws err JuMP.dual(cref)
    end

    @testset "Result attributes" begin
        err = JuMP.OptimizeNotCalled()
        model = Model(() -> MOIU.MockOptimizer(SimpleLPModel{Float64}()))
        @variable(model, x)
        c = @constraint(model, x ≤ 0)
        @objective(model, Max, x)
        @test_throws err JuMP.objective_value(model)
        @test_throws err JuMP.dual_objective_value(model)
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

        m = Model()
        err = ErrorException("Unrecognized keyword arguments: unexpected_arg")
        @test_throws err optimize!(m, unexpected_arg=1)
        JuMP.set_optimize_hook(m, (m ; my_new_arg=nothing) -> my_new_arg)
        @test optimize!(m) === nothing
        @test optimize!(m, my_new_arg = 1) == 1
    end

    @testset "UniversalFallback" begin
        m = Model()
        MOI.set(m, MOI.Test.UnknownModelAttribute(), 1)
        @test MOI.get(m, MOI.Test.UnknownModelAttribute()) == 1
    end

    @testset "Bridges" begin
        @testset "Automatic bridging" begin
            # optimizer not supporting Interval
            model = Model(() -> MOIU.MockOptimizer(SimpleLPModel{Float64}()))
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
            model = Model(() -> MOIU.MockOptimizer(SimpleLPModel{Float64}(),
                                                   needs_allocate_load=true))
            @test JuMP.bridge_constraints(model)
            @test JuMP.backend(model) isa MOIU.CachingOptimizer
            @test JuMP.backend(model).optimizer isa MOI.Bridges.LazyBridgeOptimizer
            @test JuMP.backend(model).optimizer.model isa MOIU.CachingOptimizer
            @test JuMP.backend(model).optimizer.model.optimizer isa MOIU.MockOptimizer
            @variable model x
            err = ErrorException(
                "There is no `optimizer_index` as the optimizer is not " *
                "synchronized with the cached model. Call " *
                "`MOIU.attach_optimizer(model)` to synchronize it.")
            @test_throws err optimizer_index(x)
            cref = @constraint model 0 <= x + 1 <= 1
            @test cref isa JuMP.ConstraintRef{JuMP.Model,MOI.ConstraintIndex{MOI.ScalarAffineFunction{Float64},MOI.Interval{Float64}}}
            @test_throws err optimizer_index(cref)
            JuMP.optimize!(model)
            err = ErrorException(
                "There is no `optimizer_index` for $(typeof(index(cref))) " *
                "constraints because they are bridged.")
            @test_throws err optimizer_index(cref)
        end
        @testset "Automatic bridging disabled with `bridge_constraints` keyword" begin
            model = Model(() -> MOIU.MockOptimizer(SimpleLPModel{Float64}()),
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
            function mock_factory()
                mock = MOIU.MockOptimizer(MOIU.Model{Float64}(),
                                          eval_variable_constraint_dual=false)
                optimize!(mock) = MOIU.mock_optimize!(mock, [1.0],
                        (MOI.SingleVariable, MOI.GreaterThan{Float64}) => [2.0])
                MOIU.set_mock_optimize!(mock, optimize!)
                return mock
            end
            @testset "before loading the constraint to the optimizer" begin
                @testset "optimizer set at Model" begin
                    model = Model(mock_factory)
                    @variable(model, x)
                    JuMP.add_bridge(model, NonnegativeBridge)
                    c = @constraint(model, x in Nonnegative())
                    JuMP.optimize!(model)
                    @test 1.0 == @inferred JuMP.value(x)
                    @test 1.0 == @inferred JuMP.value(c)
                    @test 2.0 == @inferred JuMP.dual(c)
                end
                @testset "optimizer set with set_optimizer" begin
                    model = Model()
                    @variable(model, x)
                    c = @constraint(model, x in Nonnegative())
                    JuMP.add_bridge(model, NonnegativeBridge)
                    set_optimizer(model, mock_factory)
                    JuMP.optimize!(model)
                    @test 1.0 == @inferred JuMP.value(x)
                    @test 1.0 == @inferred JuMP.value(c)
                    @test 2.0 == @inferred JuMP.dual(c)
                end
            end
            @testset "after loading the constraint to the optimizer" begin
                @testset "optimizer set at Model" begin
                    err = ErrorException(string("Constraints of type ",
                    "MathOptInterface.SingleVariable-in-Nonnegative are not ",
                    "supported by the solver and there are no bridges that ",
                    "can reformulate it into supported constraints."))
                    model = Model(mock_factory)
                    @variable(model, x)
                    @test_throws err @constraint(model, x in Nonnegative())
                    JuMP.add_bridge(model, NonnegativeBridge)
                    c = @constraint(model, x in Nonnegative())
                    JuMP.optimize!(model)
                    @test 1.0 == @inferred JuMP.value(x)
                    @test 1.0 == @inferred JuMP.value(c)
                    @test 2.0 == @inferred JuMP.dual(c)
                end
                @testset "optimizer set with set_optimizer" begin
                    err = MOI.UnsupportedConstraint{MOI.SingleVariable,
                                                    Nonnegative}()
                    model = Model()
                    @variable(model, x)
                    c = @constraint(model, x in Nonnegative())
                    set_optimizer(model, mock_factory)
                    @test_throws err JuMP.optimize!(model)
                    JuMP.add_bridge(model, NonnegativeBridge)
                    JuMP.optimize!(model)
                    @test 1.0 == @inferred JuMP.value(x)
                    @test 1.0 == @inferred JuMP.value(c)
                    @test 2.0 == @inferred JuMP.dual(c)
                end
            end
            @testset "automatically with BridgeableConstraint" begin
                @testset "optimizer set at Model" begin
                    model = Model(mock_factory)
                    @variable(model, x)
                    constraint = ScalarConstraint(x, Nonnegative())
                    bc = BridgeableConstraint(constraint, NonnegativeBridge)
                    c = add_constraint(model, bc)
                    JuMP.optimize!(model)
                    @test 1.0 == @inferred JuMP.value(x)
                    @test 1.0 == @inferred JuMP.value(c)
                    @test 2.0 == @inferred JuMP.dual(c)
                end
                @testset "optimizer set with set_optimizer" begin
                    model = Model()
                    @variable(model, x)
                    constraint = ScalarConstraint(x, Nonnegative())
                    bc = BridgeableConstraint(constraint, NonnegativeBridge)
                    c = add_constraint(model, bc)
                    set_optimizer(model, mock_factory)
                    JuMP.optimize!(model)
                    @test 1.0 == @inferred JuMP.value(x)
                    @test 1.0 == @inferred JuMP.value(c)
                    @test 2.0 == @inferred JuMP.dual(c)
                end
            end
        end
    end

    @testset "solve_time" begin
        @testset "NoOptimizer()" begin
            err = NoOptimizer()
            model = Model()
            @test_throws err solve_time(model)
        end

        @testset "OptimizeNotCalled()" begin
            err = OptimizeNotCalled()
            model = Model(() -> MOIU.MockOptimizer(SimpleLPModel{Float64}()))
            @test_throws err solve_time(model)
        end

        @testset "Solved model" begin
            # TODO
        end
    end

    @testset "solver_name" begin
        @testset "Not attached" begin
            model = Model()
            @test "No optimizer attached." == @inferred JuMP.solver_name(model)
        end

        @testset "Mock" begin
            model = Model(() -> MOIU.MockOptimizer(SimpleLPModel{Float64}()))
            @test "Mock" == @inferred JuMP.solver_name(model)
        end
    end
    @testset "set_silent and unset_silent" begin
        mock = MOIU.UniversalFallback(MOIU.Model{Float64}())
        model = Model(() -> MOIU.MockOptimizer(mock))
        @test JuMP.set_silent(model)
        @test MOI.get(backend(model), MOI.Silent())
        @test MOI.get(model, MOI.Silent())
        @test !JuMP.unset_silent(model)
        @test !MOI.get(backend(model), MOI.Silent())
        @test !MOI.get(model, MOI.Silent())
    end

    @testset "set_parameter" begin
        mock = MOIU.UniversalFallback(MOIU.Model{Float64}())
        model = Model(() -> MOIU.MockOptimizer(mock))
        @test JuMP.set_parameter(model, "aaa", "bbb") == "bbb"
        @test MOI.get(backend(model), MOI.RawParameter("aaa")) == "bbb"
        @test MOI.get(model, MOI.RawParameter("aaa")) == "bbb"
    end

    @testset "set_parameters" begin
        mock = MOIU.UniversalFallback(MOIU.Model{Float64}())
        model = Model(() -> MOIU.MockOptimizer(mock))
        JuMP.set_parameters(model, "aaa" => "bbb", "abc" => 10)
        @test MOI.get(model, MOI.RawParameter("aaa")) == "bbb"
        @test MOI.get(model, MOI.RawParameter("abc")) == 10
    end

    @testset "set and retrieve time limit" begin
        mock = MOIU.UniversalFallback(MOIU.Model{Float64}())
        model = Model(() -> MOIU.MockOptimizer(mock))
        JuMP.set_time_limit_sec(model, 12.0)
        @test JuMP.time_limit_sec(model) == 12.0
        JuMP.set_time_limit_sec(model, nothing)
        @test JuMP.time_limit_sec(model) === nothing
        JuMP.set_time_limit_sec(model, 12.0)
        JuMP.unset_time_limit_sec(model)
        @test JuMP.time_limit_sec(model) === nothing
    end

    @testset "set_optimizer error cases" begin
        model = Model()
        @test_throws(ErrorException(JuMP._set_optimizer_not_callable_message),
                     set_optimizer(model,
                                   MOIU.MockOptimizer(MOIU.Model{Float64}())))
        err = ErrorException("The provided optimizer_factory returned an " *
            "object of type Int64. Expected a " *
            "MathOptInterface.AbstractOptimizer.")
        @test_throws err set_optimizer(model, () -> Int64(10))
        # TODO: A factory that returns a non-empty optimizer.
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
        mock = MOIU.MockOptimizer(MOIU.Model{Float64}())
        model = JuMP.direct_model(mock)
        @test_throws ErrorException JuMP.copy(model)
    end
end

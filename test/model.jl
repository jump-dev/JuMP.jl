#  Copyright 2017, Iain Dunning, Joey Huchette, Miles Lubin, and contributors
#  This Source Code Form is subject to the terms of the Mozilla Public
#  License, v. 2.0. If a copy of the MPL was not distributed with this
#  file, You can obtain one at http://mozilla.org/MPL/2.0/.
#############################################################################
# JuMP
# An algebraic modeling language for Julia
# See http://github.com/jump-dev/JuMP.jl
#############################################################################
# test/model.jl
#############################################################################

module TestModel

using JuMP
using Test

# Simple LP model not supporting Interval
MOIU.@model(
    SimpleLPModel,
    (), (MOI.EqualTo, MOI.GreaterThan, MOI.LessThan), (), (),
    (), (MOI.ScalarAffineFunction,), (), ()
)

include("nonnegative_bridge.jl")

function _fill_small_test_model!(model)
    # The model does not need to make sense, just use many different features.
    @variable(model, a[1:5] >= 0, Int)
    @variable(model, b[6:10], Bin)
    @variable(model, c[1:3] == 0)
    @variable(model, 10 <= d[1:3] <= 20)
    @constraint(model, con1, sum(a) + sum(b) <= 5)
    @constraint(model, con2, sum(b) >= 3)
    @constraint(model, con3, sum(d[1:2]) >= 5)
    @constraint(model, con4, sum(d) <= (sum(c) + 10))
    @objective(model, Max, sum(a) - sum(b) + sum(d))
    return model
end

function test_NoOptimizer()
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

function test_result_attributes()
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

function test_result_attributes_after_empty()
    err = JuMP.OptimizeNotCalled()
    model = Model(() -> MOIU.MockOptimizer(SimpleLPModel{Float64}()))
    @variable(model, x)
    c = @constraint(model, x ≤ 0)
    @objective(model, Max, x)
    optimize!(model)
    empty!(model)
    @test_throws err JuMP.objective_value(model)
    @test_throws err JuMP.dual_objective_value(model)
    @test_throws err JuMP.objective_bound(model)
    @test_throws err JuMP.value(x)
    @test_throws err JuMP.value(c)
    @test_throws err JuMP.dual(c)
end

function test_empty!()
    model = Model()
    backend_type = typeof(backend(model))
    model.optimize_hook === nothing
    hook(m) = nothing
    JuMP.set_optimize_hook(model, hook)
    @test model.optimize_hook === hook
    @test _fill_small_test_model!(model) === model
    @test_throws ErrorException _fill_small_test_model!(model)
    @test empty!(model) === model
    @test model.optimize_hook === hook # empty! does not touch the hook
    @test isa(backend(model), backend_type)
    @test _fill_small_test_model!(model) === model
end

function test_variable_model_hygiene()
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

function test_optimize_hook()
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

function test_universal_fallback()
    m = Model()
    MOI.set(m, MOI.Test.UnknownModelAttribute(), 1)
    @test MOI.get(m, MOI.Test.UnknownModelAttribute()) == 1
end

function test_automatic_bridging()
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

function test_automatic_bridging_with_cache()
    # optimizer not supporting Interval and not supporting `default_copy_to`
    model = Model(() -> MOIU.MockOptimizer(SimpleLPModel{Float64}(),
                                            needs_allocate_load=true))
    @test JuMP.bridge_constraints(model)
    @test JuMP.backend(model) isa MOIU.CachingOptimizer
    @test JuMP.backend(model).optimizer isa MOI.Bridges.LazyBridgeOptimizer
    @test JuMP.backend(model).optimizer.model isa MOIU.CachingOptimizer
    @test JuMP.backend(model).optimizer.model.optimizer isa MOIU.MockOptimizer
    @variable(model, x)
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

function test_automatic_bridging_with_bridge_constraints_false()
    model = Model(
        () -> MOIU.MockOptimizer(SimpleLPModel{Float64}()),
        bridge_constraints = false
    )
    @test !JuMP.bridge_constraints(model)
    @test JuMP.backend(model) isa MOIU.CachingOptimizer
    @test !(JuMP.backend(model).optimizer isa MOI.Bridges.LazyBridgeOptimizer)
    @variable model x
    err = ErrorException("Constraints of type MathOptInterface.ScalarAffineFunction{Float64}-in-MathOptInterface.Interval{Float64} are not supported by the solver, try using `bridge_constraints=true` in the `JuMP.Model` constructor if you believe the constraint can be reformulated to constraints supported by the solver.")
    @test_throws err @constraint model 0 <= x + 1 <= 1
end

function test_no_bridge_direct_model()
    optimizer = MOIU.MockOptimizer(SimpleLPModel{Float64}())
    model = JuMP.direct_model(optimizer)
    @test !JuMP.bridge_constraints(model)
    @variable model x
    err = ErrorException("Constraints of type MathOptInterface.ScalarAffineFunction{Float64}-in-MathOptInterface.Interval{Float64} are not supported by the solver.")
    @test_throws err @constraint model 0 <= x + 1 <= 1
end

function _mock_factory()
    mock = MOIU.MockOptimizer(
        MOIU.Model{Float64}(), eval_variable_constraint_dual = false
    )
    optimize!(mock) = MOIU.mock_optimize!(
        mock, [1.0], (MOI.SingleVariable, MOI.GreaterThan{Float64}) => [2.0]
    )
    MOIU.set_mock_optimize!(mock, optimize!)
    return mock
end

function test_addbridge_at_model()
    model = Model(_mock_factory)
    @variable(model, x)
    JuMP.add_bridge(model, NonnegativeBridge)
    c = @constraint(model, x in Nonnegative())
    JuMP.optimize!(model)
    @test 1.0 == @inferred JuMP.value(x)
    @test 1.0 == @inferred JuMP.value(c)
    @test 2.0 == @inferred JuMP.dual(c)
end

function test_addbridge_at_set_optimizer()
    model = Model()
    @variable(model, x)
    c = @constraint(model, x in Nonnegative())
    JuMP.add_bridge(model, NonnegativeBridge)
    set_optimizer(model, _mock_factory)
    JuMP.optimize!(model)
    @test 1.0 == @inferred JuMP.value(x)
    @test 1.0 == @inferred JuMP.value(c)
    @test 2.0 == @inferred JuMP.dual(c)
end

function test_bridging_at_model()
    err = ErrorException("""
Constrained variables in `Nonnegative` are not supported and cannot be bridged into supported constrained variables and constraints. See details below:
[1] constrained variables in `Nonnegative` are not supported because no added bridge supports bridging it.
  Cannot add free variables and then constrain them because:
  (1) `MOI.SingleVariable`-in-`Nonnegative` constraints are not supported
(1) `MOI.SingleVariable`-in-`Nonnegative` constraints are not supported because:
  Cannot use `MOIB.Constraint.ScalarFunctionizeBridge{Float64,Nonnegative}` because:
  (2) `MOI.ScalarAffineFunction{Float64}`-in-`Nonnegative` constraints are not supported
(2) `MOI.ScalarAffineFunction{Float64}`-in-`Nonnegative` constraints are not supported because:
  Cannot use `MOIB.Constraint.ScalarSlackBridge{Float64,MOI.ScalarAffineFunction{Float64},Nonnegative}` because:
  [1] constrained variables in `Nonnegative` are not supported
""")

    model = Model(_mock_factory)
    @variable(model, x)
    @test_throws err try
        @constraint(model, x in Nonnegative())
    catch e
        error(replace(e.msg, "Main.TestModel." => ""))
    end
    JuMP.add_bridge(model, NonnegativeBridge)
    c = @constraint(model, x in Nonnegative())
    JuMP.optimize!(model)
    @test 1.0 == @inferred JuMP.value(x)
    @test 1.0 == @inferred JuMP.value(c)
    @test 2.0 == @inferred JuMP.dual(c)
end

function test_bridging_at_set_optimizer()
    err = MOI.UnsupportedConstraint{MOI.SingleVariable, Nonnegative}()
    model = Model()
    @variable(model, x)
    c = @constraint(model, x in Nonnegative())
    set_optimizer(model, _mock_factory)
    @test_throws err JuMP.optimize!(model)
    JuMP.add_bridge(model, NonnegativeBridge)
    JuMP.optimize!(model)
    @test 1.0 == @inferred JuMP.value(x)
    @test 1.0 == @inferred JuMP.value(c)
    @test 2.0 == @inferred JuMP.dual(c)
end

function test_bridgeableconstraint_at_model()
    model = Model(_mock_factory)
    @variable(model, x)
    constraint = ScalarConstraint(x, Nonnegative())
    bc = BridgeableConstraint(constraint, NonnegativeBridge)
    c = add_constraint(model, bc)
    JuMP.optimize!(model)
    @test 1.0 == @inferred JuMP.value(x)
    @test 1.0 == @inferred JuMP.value(c)
    @test 2.0 == @inferred JuMP.dual(c)
end

function test_bridgeableconstraint_at_set_optimizer()
    model = Model()
    @variable(model, x)
    constraint = ScalarConstraint(x, Nonnegative())
    bc = BridgeableConstraint(constraint, NonnegativeBridge)
    c = add_constraint(model, bc)
    set_optimizer(model, _mock_factory)
    JuMP.optimize!(model)
    @test 1.0 == @inferred JuMP.value(x)
    @test 1.0 == @inferred JuMP.value(c)
    @test 2.0 == @inferred JuMP.dual(c)
end

function test_solve_time_no_optimizer()
    err = NoOptimizer()
    model = Model()
    @test_throws err solve_time(model)
end

function test_solve_time_OptimizeNotCalled()
    err = OptimizeNotCalled()
    model = Model(() -> MOIU.MockOptimizer(SimpleLPModel{Float64}()))
    @test_throws err solve_time(model)
end

function test_solve_time()
    # TODO(odow): Add a test for getting the solve time of a solved model.
    @test_broken true == false
end

function test_solver_name_not_attached()
    model = Model()
    @test "No optimizer attached." == @inferred JuMP.solver_name(model)
end

function test_solver_name_mock()
    model = Model(() -> MOIU.MockOptimizer(SimpleLPModel{Float64}()))
    @test "Mock" == @inferred JuMP.solver_name(model)
end

function test_set_silent()
    mock = MOIU.UniversalFallback(MOIU.Model{Float64}())
    model = Model(() -> MOIU.MockOptimizer(mock))
    @test JuMP.set_silent(model)
    @test MOI.get(backend(model), MOI.Silent())
    @test MOI.get(model, MOI.Silent())
    @test !JuMP.unset_silent(model)
    @test !MOI.get(backend(model), MOI.Silent())
    @test !MOI.get(model, MOI.Silent())
end

function test_set_optimizer_attribute()
    mock = MOIU.UniversalFallback(MOIU.Model{Float64}())
    model = Model(() -> MOIU.MockOptimizer(mock))
    @test JuMP.set_optimizer_attribute(model, "aaa", "bbb") == "bbb"
    @test MOI.get(backend(model), MOI.RawParameter("aaa")) == "bbb"
    @test MOI.get(model, MOI.RawParameter("aaa")) == "bbb"
end

function test_set_optimizer_attributes()
    mock = MOIU.UniversalFallback(MOIU.Model{Float64}())
    model = Model(() -> MOIU.MockOptimizer(mock))
    JuMP.set_optimizer_attributes(model, "aaa" => "bbb", "abc" => 10)
    @test MOI.get(model, MOI.RawParameter("aaa")) == "bbb"
    @test MOI.get(model, MOI.RawParameter("abc")) == 10
end

function test_get_optimizer_attribute()
    mock = MOIU.UniversalFallback(MOIU.Model{Float64}())
    model = Model(() -> MOIU.MockOptimizer(mock))
    @test JuMP.set_optimizer_attribute(model, "aaa", "bbb") == "bbb"
    @test JuMP.get_optimizer_attribute(model, "aaa") == "bbb"
end

function test_time_limit_sec()
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

struct DummyExtensionData
    model::JuMP.Model
end

function JuMP.copy_extension_data(
    data::DummyExtensionData,
    new_model::JuMP.AbstractModel,
    model::JuMP.AbstractModel
)
    @test data.model === model
    return DummyExtensionData(new_model)
end

function dummy_optimizer_hook(::JuMP.AbstractModel) end

function test_copy_model()
    for caching_mode in (MOIU.AUTOMATIC, MOIU.MANUAL)
        model = Model(caching_mode = caching_mode)
        model.optimize_hook = dummy_optimizer_hook
        model.ext[:dummy] = DummyExtensionData(model)
        @variable(model, x ≥ 0, Bin)
        @variable(model, y ≤ 1, Int)
        @variable(model, z == 0)
        @constraint(model, cref, x + y == 1)
        new_model, reference_map = JuMP.copy_model(model)
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
        @test reference_map[JuMP.LowerBoundRef(x)] == @inferred JuMP.LowerBoundRef(x_new)
        @test reference_map[JuMP.BinaryRef(x)] == @inferred JuMP.BinaryRef(x_new)
        @test reference_map[JuMP.UpperBoundRef(y)] == @inferred JuMP.UpperBoundRef(y_new)
        @test reference_map[JuMP.IntegerRef(y)] == @inferred JuMP.IntegerRef(y_new)
        @test reference_map[JuMP.FixRef(z)] == @inferred JuMP.FixRef(z_new)
        cref_new = reference_map[cref]
        @test cref_new.model === new_model
        @test "cref" == @inferred JuMP.name(cref_new)
    end
end

function test_base_copy()
    caching_mode = MOIU.MANUAL
    model = Model(caching_mode = caching_mode)
    @variable(model, x)
    new_model = Base.copy(model)
    @test caching_mode == @inferred MOIU.mode(JuMP.backend(new_model))
    @test JuMP.owner_model(new_model[:x]) === new_model
    @test JuMP.owner_model(new_model[:x]) != model
    @test "x" == @inferred JuMP.name(new_model[:x])
end

function test_copy_in_direct_mode()
    mock = MOIU.MockOptimizer(MOIU.Model{Float64}())
    model = JuMP.direct_model(mock)
    @test_throws ErrorException JuMP.copy(model)
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

TestModel.runtests()

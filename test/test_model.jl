#  Copyright 2017, Iain Dunning, Joey Huchette, Miles Lubin, and contributors
#  This Source Code Form is subject to the terms of the Mozilla Public
#  License, v. 2.0. If a copy of the MPL was not distributed with this
#  file, You can obtain one at https://mozilla.org/MPL/2.0/.
#############################################################################
# JuMP
# An algebraic modeling language for Julia
# See https://github.com/jump-dev/JuMP.jl
#############################################################################
# test/model.jl
#############################################################################

module TestModels

using JuMP
using Test

# Simple LP model not supporting Interval
MOIU.@model(
    SimpleLPModel,
    (),
    (MOI.EqualTo, MOI.GreaterThan, MOI.LessThan),
    (),
    (),
    (),
    (MOI.ScalarAffineFunction,),
    (),
    ()
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

function _test_result_attributes(; test_empty = false)
    err = OptimizeNotCalled()
    model = Model(() -> MOIU.MockOptimizer(SimpleLPModel{Float64}()))
    @variable(model, x)
    c = @constraint(model, x ≤ 0)
    @objective(model, Max, x)
    if test_empty
        optimize!(model)
        @test !isempty(model)
        empty!(model)
        @test isempty(model)
    end
    @test_throws err objective_value(model)
    @test_throws err dual_objective_value(model)
    @test_throws err objective_bound(model)
    @test_throws err value(x)
    @test_throws err value(c)
    @test_throws err dual(c)
end

function test_result_attributes()
    return _test_result_attributes()
end

function test_result_attributes_after_empty()
    return _test_result_attributes(; test_empty = true)
end

function fill_small_test_model!(model)
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

function test_nooptimizer()
    err = NoOptimizer()
    model = Model()
    @variable(model, x)
    @test_throws err optimizer_index(x)
    cref = @constraint(model, x == 1)
    @test_throws err optimizer_index(cref)
    @test_throws err optimize!(model)
    @test_throws OptimizeNotCalled() value(x)
    @test_throws OptimizeNotCalled() value(cref)
    @test_throws OptimizeNotCalled() dual(cref)
    return
end

function test_empty!_model()
    model = Model()
    backend_type = typeof(backend(model))
    model.optimize_hook === nothing
    hook(m) = nothing
    set_optimize_hook(model, hook)
    @test model.optimize_hook === hook
    @test fill_small_test_model!(model) === model
    @test_throws ErrorException fill_small_test_model!(model)
    @test empty!(model) === model
    @test model.optimize_hook === hook # empty! does not touch the hook
    @test isa(backend(model), backend_type)
    @test fill_small_test_model!(model) === model
end

function test_innermost_error()
    model = Model()
    @test_throws ErrorException unsafe_backend(model)
end

function test_innermost()
    model = Model(() -> MOIU.MockOptimizer(MOIU.Model{Float64}()))
    MOI.Utilities.attach_optimizer(model)
    @test unsafe_backend(model) isa MOIU.MockOptimizer{MOIU.Model{Float64}}
end

function test_hygiene_variable()
    model_x = Model()
    @variable(model_x, x)
    model_y = Model()
    @variable(model_y, y)
    err = VariableNotOwned(y)
    @test_throws err @constraint(model_x, y in MOI.EqualTo(1.0))
    @test_throws err @constraint(model_x, [x, y] in MOI.Zeros(2))
    @test_throws err @objective(model_x, Min, y)
end

function test_hygiene_linear()
    model_x = Model()
    @variable(model_x, x)
    model_y = Model()
    @variable(model_y, y)
    err = VariableNotOwned(y)

    @test_throws err @constraint(model_x, x + y == 1)
    @test_throws err @constraint(model_x, [x, x + y] in MOI.Zeros(2))
    @test_throws err @constraint(model_x, [x + y, x] in MOI.Zeros(2))
    @test_throws err @objective(model_x, Min, x + y)
end

function test_hygiene_quadratic()
    model_x = Model()
    @variable(model_x, x)
    model_y = Model()
    @variable(model_y, y)
    err = VariableNotOwned(y)

    @test_throws err @constraint(model_x, x * y >= 1)
    @test_throws err @constraint(model_x, [x, x * y] in MOI.Zeros(2))
    @test_throws err @constraint(model_x, [x * y, x] in MOI.Zeros(2))
    @test_throws err @constraint(model_x, x * x + x + y <= 1)
    @test_throws err @constraint(model_x, [x, x * x + x + y] in MOI.Zeros(2))
    @test_throws err @constraint(model_x, [x * x + x + y, x] in MOI.Zeros(2))
    @test_throws err @objective(model_x, Min, x * y)
    @test_throws err @objective(model_x, Min, x * x + x + y)
end

function test_hygiene_attribute()
    model_x = Model()
    @variable(model_x, x)
    model_y = Model()
    @variable(model_y, y)
    err = VariableNotOwned(y)

    cy = @constraint(model_y, y in MOI.EqualTo(1.0))
    cerr = ConstraintNotOwned(cy)

    @test_throws err MOI.get(model_x, MOI.VariablePrimalStart(), y)
    @test_throws cerr MOI.get(model_x, MOI.ConstraintPrimalStart(), cy)
    @test_throws err MOI.set(model_x, MOI.VariablePrimalStart(), y, 1.0)
    @test_throws cerr MOI.set(model_x, MOI.ConstraintPrimalStart(), cy, 1.0)
end

function test_optimize_hook()
    m = Model()
    @test m.optimize_hook === nothing
    called = false
    function hook(m)
        return called = true
    end
    set_optimize_hook(m, hook)
    @test !called
    optimize!(m)
    @test called

    m = Model()
    err = ErrorException("Unrecognized keyword arguments: unexpected_arg")
    @test_throws err optimize!(m, unexpected_arg = 1)
    set_optimize_hook(m, (m; my_new_arg = nothing) -> my_new_arg)
    @test optimize!(m) === nothing
    @test optimize!(m; my_new_arg = 1) == 1
end

function test_universal_fallback()
    m = Model()
    MOI.set(m, MOI.Test.UnknownModelAttribute(), 1)
    @test MOI.get(m, MOI.Test.UnknownModelAttribute()) == 1
end

function test_bridges_automatic()
    # optimizer not supporting Interval
    model = Model(() -> MOIU.MockOptimizer(SimpleLPModel{Float64}()))
    @test bridge_constraints(model)
    @test backend(model) isa MOIU.CachingOptimizer
    @test backend(model).optimizer isa MOI.Bridges.LazyBridgeOptimizer
    @test backend(model).optimizer.model isa MOIU.MockOptimizer
    @variable model x
    cref = @constraint model 0 <= x + 1 <= 1
    @test cref isa ConstraintRef{
        Model,
        MOI.ConstraintIndex{
            MOI.ScalarAffineFunction{Float64},
            MOI.Interval{Float64},
        },
    }
    return optimize!(model)
end

function test_bridges_automatic_disabled()
    # Automatic bridging disabled with `add_bridges` keyword
    model = Model(
        () -> MOIU.MockOptimizer(SimpleLPModel{Float64}());
        add_bridges = false,
    )
    @test !bridge_constraints(model)
    @test backend(model) isa MOIU.CachingOptimizer
    @test !(backend(model).optimizer isa MOI.Bridges.LazyBridgeOptimizer)
    @variable model x
    F = MOI.ScalarAffineFunction{Float64}
    S = MOI.Interval{Float64}
    err = ErrorException(
        "Constraints of type $(F)-in-$(S) are not supported by the " *
        "solver.\n\nIf you expected the solver to support your problem, " *
        "you may have an error in your formulation. Otherwise, consider " *
        "using a different solver.\n\nThe list of available solvers, " *
        "along with the problem types they support, is available at " *
        "https://jump.dev/JuMP.jl/stable/installation/#Supported-solvers.",
    )
    @test_throws err @constraint model 0 <= x + 1 <= 1
end

function test_bridges_direct()
    # No bridge automatically added in Direct mode
    optimizer = MOIU.MockOptimizer(SimpleLPModel{Float64}())
    model = direct_model(optimizer)
    @test !bridge_constraints(model)
    @variable model x
    F = MOI.ScalarAffineFunction{Float64}
    S = MOI.Interval{Float64}
    err = ErrorException(
        "Constraints of type $(F)-in-$(S) are not supported by the " *
        "solver.\n\nIf you expected the solver to support your problem, " *
        "you may have an error in your formulation. Otherwise, consider " *
        "using a different solver.\n\nThe list of available solvers, " *
        "along with the problem types they support, is available at " *
        "https://jump.dev/JuMP.jl/stable/installation/#Supported-solvers.",
    )
    @test_throws err @constraint model 0 <= x + 1 <= 1
end

function mock_factory()
    mock = MOIU.MockOptimizer(
        MOIU.Model{Float64}();
        eval_variable_constraint_dual = false,
    )
    function optimize!(mock)
        return MOIU.mock_optimize!(
            mock,
            [1.0],
            (MOI.VariableIndex, MOI.GreaterThan{Float64}) => [2.0],
        )
    end
    MOIU.set_mock_optimize!(mock, optimize!)
    return mock
end

function test_bridges_add_before_con_model_optimizer()
    model = Model(mock_factory)
    @variable(model, x)
    add_bridge(model, NonnegativeBridge)
    c = @constraint(model, x in Nonnegative())
    optimize!(model)
    @test 1.0 == @inferred value(x)
    @test 1.0 == @inferred value(c)
    @test 2.0 == @inferred dual(c)
end

function test_bridges_add_before_con_set_optimizer()
    model = Model()
    @variable(model, x)
    c = @constraint(model, x in Nonnegative())
    add_bridge(model, NonnegativeBridge)
    set_optimizer(model, mock_factory)
    optimize!(model)
    @test 1.0 == @inferred value(x)
    @test 1.0 == @inferred value(c)
    @test 2.0 == @inferred dual(c)
end

function test_bridges_remove_bridge()
    model = Model(mock_factory)
    add_bridge(model, NonnegativeBridge)
    @variable(model, x)
    c = @constraint(model, x in Nonnegative())
    optimize!(model)
    @test 1.0 == @inferred value(x)
    @test 1.0 == @inferred value(c)
    @test 2.0 == @inferred dual(c)
    remove_bridge(model, NonnegativeBridge)
    @test_throws MOI.UnsupportedConstraint optimize!(model)
    return
end

function test_bridges_print_active_bridge()
    model = Model(mock_factory)
    add_bridge(model, NonnegativeBridge)
    @variable(model, x)
    c = @constraint(model, x in Nonnegative())
    optimize!(model)
    print_active_bridges(model)
    contents = sprint(print_active_bridges, model)
    @test occursin("Nonnegative", contents)
    @test occursin("NonnegativeBridge", contents)
    @test occursin("GreaterThan", contents)
    return
end

function test_bridges_print_active_bridge_objective()
    model = Model(mock_factory)
    add_bridge(model, NonnegativeBridge)
    contents = sprint(print_active_bridges, model, VariableRef)
    @test occursin(" * Supported", contents)
    @test occursin("objective", contents)
    @test occursin("VariableIndex", contents)
    return
end

function test_bridges_print_active_bridge_constraint()
    model = Model(mock_factory)
    add_bridge(model, NonnegativeBridge)
    contents = sprint(print_active_bridges, model, VariableRef, Nonnegative)
    @test occursin("Nonnegative", contents)
    @test occursin("NonnegativeBridge", contents)
    @test occursin("GreaterThan", contents)
    return
end

function test_bridges_print_active_bridge_variable()
    model = Model(mock_factory)
    contents = sprint(print_active_bridges, model, MOI.ZeroOne)
    @test occursin(" * Supported", contents)
    @test occursin("ZeroOne", contents)
    return
end

function test_bridges_add_after_con_model_optimizer()
    model = Model(mock_factory)
    @variable(model, x)
    flag = true
    try
        @constraint(model, x in Nonnegative())
        flag = false
    catch err
        @test err isa ErrorException
        # Rather than test a particular bridging error, just check that the
        # bridge explanation has been called. The sequence of errors could vary
        # between MOI versions.
        @test occursin("Nonnegative", err.msg)
        @test occursin("are not supported", err.msg)
    end
    @test flag
    add_bridge(model, NonnegativeBridge)
    c = @constraint(model, x in Nonnegative())
    optimize!(model)
    @test 1.0 == @inferred value(x)
    @test 1.0 == @inferred value(c)
    @test 2.0 == @inferred dual(c)
end

function test_bridges_add_after_con_set_optimizer()
    err = MOI.UnsupportedConstraint{MOI.VariableIndex,Nonnegative}()
    model = Model()
    @variable(model, x)
    c = @constraint(model, x in Nonnegative())
    set_optimizer(model, mock_factory)
    @test_throws err optimize!(model)
    add_bridge(model, NonnegativeBridge)
    optimize!(model)
    @test 1.0 == @inferred value(x)
    @test 1.0 == @inferred value(c)
    @test 2.0 == @inferred dual(c)
end

function test_bridges_add_bridgeable_con_model_optimizer()
    model = Model(mock_factory)
    @variable(model, x)
    constraint = ScalarConstraint(x, Nonnegative())
    bc = BridgeableConstraint(constraint, NonnegativeBridge)
    c = add_constraint(model, bc)
    optimize!(model)
    @test 1.0 == @inferred value(x)
    @test 1.0 == @inferred value(c)
    @test 2.0 == @inferred dual(c)
end

function test_bridges_add_bridgeable_con_set_optimizer()
    model = Model()
    @variable(model, x)
    constraint = ScalarConstraint(x, Nonnegative())
    bc = BridgeableConstraint(constraint, NonnegativeBridge)
    c = add_constraint(model, bc)
    set_optimizer(model, mock_factory)
    optimize!(model)
    @test 1.0 == @inferred value(x)
    @test 1.0 == @inferred value(c)
    @test 2.0 == @inferred dual(c)
end

function test_bridge_graph_false()
    model = Model(mock_factory; add_bridges = false)
    @variable(model, x)
    @test_throws(
        ErrorException(
            "Cannot use bridge if `add_bridges` was set to `false` in " *
            "the `Model` constructor.",
        ),
        add_bridge(model, NonnegativeBridge)
    )
    @test_throws(
        ErrorException(
            "Cannot use bridge if `add_bridges` was set to " *
            "`false` in the `Model` constructor.",
        ),
        print_bridge_graph(model)
    )
    optimize!(model)
    @test 1.0 == @inferred value(x)
end

function test_bridge_graph_true()
    model = Model(mock_factory)
    @variable(model, x)
    add_bridge(model, NonnegativeBridge)
    @test sprint(print_bridge_graph, model) ==
          "Bridge graph with 0 variable nodes, 0 constraint nodes and 0 objective nodes.\n"
    c = @constraint(model, x in Nonnegative())
    @test sprint(print_bridge_graph, model) == replace(
        "Bridge graph with 1 variable nodes, 2 constraint nodes and 0 objective nodes.\n" *
        " [1] constrained variables in `$(Nonnegative)` are supported (distance 2) by adding free variables and then constrain them, see (1).\n" *
        " (1) `MOI.VariableIndex`-in-`$(Nonnegative)` constraints are bridged (distance 1) by $(NonnegativeBridge{Float64,MOI.VariableIndex}).\n" *
        " (2) `MOI.ScalarAffineFunction{Float64}`-in-`$(Nonnegative)` constraints are bridged (distance 1) by $(NonnegativeBridge{Float64,MOI.ScalarAffineFunction{Float64}}).\n",
        "MathOptInterface." => "MOI.",
    )
    optimize!(model)
    @test 1.0 == @inferred value(x)
    @test 1.0 == @inferred value(c)
    @test 2.0 == @inferred dual(c)
end

function test_solve_time()
    err = NoOptimizer()
    model = Model()
    @test_throws err solve_time(model)

    err = OptimizeNotCalled()
    model = Model(() -> MOIU.MockOptimizer(SimpleLPModel{Float64}()))
    @test_throws err solve_time(model)

    # TODO: Solved model
end

function test_solver_name()
    model = Model()
    @test "No optimizer attached." == @inferred solver_name(model)

    model = Model(() -> MOIU.MockOptimizer(SimpleLPModel{Float64}()))
    @test "Mock" == @inferred solver_name(model)
end

function test_set_silent()
    mock = MOIU.UniversalFallback(MOIU.Model{Float64}())
    model = Model(() -> MOIU.MockOptimizer(mock))
    set_silent(model)
    @test MOI.get(backend(model), MOI.Silent())
    @test MOI.get(model, MOI.Silent())
    @test get_attribute(model, MOI.Silent())
    unset_silent(model)
    @test !MOI.get(backend(model), MOI.Silent())
    @test !MOI.get(model, MOI.Silent())
    @test !get_attribute(model, MOI.Silent())
    return
end

function test_set_optimizer_attribute()
    mock = MOIU.UniversalFallback(MOIU.Model{Float64}())
    model = Model(() -> MOIU.MockOptimizer(mock))
    @test set_optimizer_attribute(model, "aaa", "bbb") === nothing
    @test MOI.get(backend(model), MOI.RawOptimizerAttribute("aaa")) == "bbb"
    @test MOI.get(model, MOI.RawOptimizerAttribute("aaa")) == "bbb"
end

function test_set_optimizer_attributes()
    mock = MOIU.UniversalFallback(MOIU.Model{Float64}())
    model = Model(() -> MOIU.MockOptimizer(mock))
    set_optimizer_attributes(model, "aaa" => "bbb", "abc" => 10)
    @test MOI.get(model, MOI.RawOptimizerAttribute("aaa")) == "bbb"
    @test MOI.get(model, MOI.RawOptimizerAttribute("abc")) == 10
end

function test_get_optimizer_attribute()
    mock = MOIU.UniversalFallback(MOIU.Model{Float64}())
    model = Model(() -> MOIU.MockOptimizer(mock))
    @test set_optimizer_attribute(model, "aaa", "bbb") === nothing
    @test get_optimizer_attribute(model, "aaa") == "bbb"
end

function test_optimizer_attribute_abstract_string()
    mock = MOIU.UniversalFallback(MOIU.Model{Float64}())
    model = Model(() -> MOIU.MockOptimizer(mock))
    attribute = "abc"
    abstract_str = @views attribute[1:end]
    set_optimizer_attribute(model, abstract_str, "x")
    @test get_optimizer_attribute(model, abstract_str) == "x"
    @test get_optimizer_attribute(model, attribute) == "x"
    return
end

function test_set_retrieve_time_limit()
    mock = MOIU.UniversalFallback(MOIU.Model{Float64}())
    model = Model(() -> MOIU.MockOptimizer(mock))
    set_time_limit_sec(model, 12.0)
    @test time_limit_sec(model) == 12.0
    set_time_limit_sec(model, nothing)
    @test time_limit_sec(model) === nothing
    set_time_limit_sec(model, 12.0)
    unset_time_limit_sec(model)
    @test time_limit_sec(model) === nothing
    # Test setting something other than `Float64`.
    set_time_limit_sec(model, 0x0a)
    @test time_limit_sec(model) == 10.0
    return
end

struct DummyExtensionData
    model::GenericModel
end
function JuMP.copy_extension_data(
    data::DummyExtensionData,
    new_model::AbstractModel,
    model::AbstractModel,
)
    @test data.model === model
    return DummyExtensionData(new_model)
end
function dummy_optimizer_hook(::AbstractModel) end

function copy_model_style_mode(use_copy_model, filter_mode)
    model = Model()
    model.optimize_hook = dummy_optimizer_hook
    data = DummyExtensionData(model)
    model.ext[:dummy] = data
    @variable(model, w[i = 1:2] ≥ 0)
    @variable(model, x ≥ 0, Bin)
    @variable(model, y ≤ 1, Int)
    @variable(model, z == 0)
    @constraint(model, cref, x + y == 1)
    @constraint(model, cref2[i = 1:2], w[i] + z == 1)

    if use_copy_model
        if filter_mode
            filter_constraints = (cr) -> cr != cref
            new_model, reference_map =
                copy_model(model; filter_constraints = filter_constraints)
        else
            new_model, reference_map = copy_model(model)
        end
    else
        new_model = copy(model)
        reference_map = Dict{
            Union{VariableRef,ConstraintRef},
            Union{VariableRef,ConstraintRef},
        }()
        reference_map[w[1]] = new_model[:w][1]
        reference_map[w[2]] = new_model[:w][2]
        reference_map[x] = new_model[:x]
        reference_map[y] = new_model[:y]
        reference_map[z] = new_model[:z]
        reference_map[cref] = new_model[:cref]
        reference_map[cref2[1]] = new_model[:cref2][1]
        reference_map[cref2[2]] = new_model[:cref2][2]
    end
    @test new_model.optimize_hook === dummy_optimizer_hook
    @test new_model.ext[:dummy].model === new_model
    x_new = reference_map[x]
    @test owner_model(x_new) === new_model
    @test "x" == @inferred name(x_new)
    y_new = reference_map[y]
    @test owner_model(y_new) === new_model
    @test "y" == @inferred name(y_new)
    z_new = reference_map[z]
    @test owner_model(z_new) === new_model
    @test "z" == @inferred name(z_new)
    if use_copy_model
        @test reference_map[LowerBoundRef(x)] == @inferred LowerBoundRef(x_new)
        @test reference_map[BinaryRef(x)] == @inferred BinaryRef(x_new)
        @test reference_map[UpperBoundRef(y)] == @inferred UpperBoundRef(y_new)
        @test reference_map[IntegerRef(y)] == @inferred IntegerRef(y_new)
        @test reference_map[FixRef(z)] == @inferred FixRef(z_new)
    end

    cref2_1_new = reference_map[cref2[1]]
    @test cref2_1_new.model === new_model
    @test "cref2[1]" == @inferred name(cref2_1_new)
    cref2_2_new = reference_map[cref2[2]]
    @test cref2_2_new.model === new_model
    @test "cref2[2]" == @inferred name(cref2_2_new)

    if filter_mode
        @test_throws KeyError object_dictionary(new_model)[name(cref)]
    else
        cref_new = reference_map[cref]
        @test cref_new.model === new_model
        @test "cref" == @inferred name(cref_new)
    end
end

function test_copy_model_jump_auto()
    return copy_model_style_mode(true, false)
end

function test_compute_conflict()
    err = NoOptimizer()
    model = Model()
    @test_throws err compute_conflict!(model)
end

function test_copy_model_base_auto()
    return copy_model_style_mode(false, false)
end

function test_copy_direct_mode()
    mock = MOIU.MockOptimizer(MOIU.Model{Float64}())
    model = direct_model(mock)
    @test_throws ErrorException copy(model)
end

function test_direct_mode_using_OptimizerWithAttributes()
    function fake_optimizer()
        return MOIU.MockOptimizer(
            MOI.Utilities.UniversalFallback(MOI.Utilities.Model{Float64}()),
        )
    end

    optimizer = optimizer_with_attributes(fake_optimizer, "a" => 1, "b" => 2)
    model = direct_model(optimizer)
    @test model.moi_backend isa MOIU.MockOptimizer
    @test MOI.get(model.moi_backend, MOI.RawOptimizerAttribute("a")) == 1
    @test MOI.get(model.moi_backend, MOI.RawOptimizerAttribute("b")) == 2
end

function test_copy_expr_aff()
    model = Model()
    @variable(model, x)
    @expression(model, ex, 2 * x + 1)
    new_model, ref_map = copy_model(model)
    @test ref_map[ex] == new_model[:ex]
    @test new_model[:ex] == 2 * new_model[:x] + 1
end

function test_copy_expr_quad()
    model = Model()
    @variable(model, x)
    @expression(model, ex, 2 * x^2 + x + 1)
    new_model, ref_map = copy_model(model)
    @test ref_map[ex] == new_model[:ex]
    @test new_model[:ex] == 2 * new_model[:x]^2 + new_model[:x] + 1
end

function test_copy_dense()
    model = Model()
    @variable(model, x[[:a, :b]])
    new_model, ref_map = copy_model(model)
    @test ref_map[x] == new_model[:x]
end

function test_copy_sparse()
    model = Model()
    @variable(model, x[i = 1:4; isodd(i)])
    new_model, ref_map = copy_model(model)
    @test ref_map[x] == new_model[:x]
end

function test_copy_object_fail()
    model = Model()
    @variable(model, x)
    model[:d] = Dict(x => 2)
    new_model, ref_map = copy_model(model)
    @test !haskey(new_model, :d)
end

function test_copy_extension()
    model = Model()
    @variable(model, x)
    model.ext[:my_content] = "x"
    new_model, _ = copy_model(model)
    @test ismissing(new_model.ext[:my_content])
    return
end

function test_haskey()
    model = Model()
    @variable(model, p[i = 1:10] >= 0)
    @test haskey(model, :p)
    @test !haskey(model, :i)
end

function test_copy_refmap_expr()
    model = Model()
    @variable(model, x)
    @expression(model, expr[i = 1:2, j = 1:2], [i, j] * x)
    new_model, ref_map = copy_model(model)
    @test ref_map[expr] == new_model[:expr]
    for i in 1:2, j in 1:2
        @test new_model[:expr][i, j] == [i, j] * new_model[:x]
    end
end

function test_copy_dict_expr()
    model = Model()
    @variable(model, x)
    @expression(model, dictExpr[i = 1:2, j = 1:2], Dict(i => x, j => 2x))
    new_model, ref_map = copy_model(model)
    @test !haskey(new_model, :dictExpr)
end

function test_copy_filter()
    return copy_model_style_mode(true, true)
end

function test_copy_filter_array()
    model = Model()
    @variable(model, x[i = 1:2], container = Array)
    @constraint(model, cref[i = 1:2], x[i] == 1, container = Array)
    @test num_constraints(
        model,
        GenericAffExpr{Float64,VariableRef},
        MOI.EqualTo{Float64},
    ) == 2

    filter_constraints = (cr) -> cr != cref[1]
    new_model, reference_map =
        copy_model(model; filter_constraints = filter_constraints)
    @test num_constraints(
        new_model,
        GenericAffExpr{Float64,VariableRef},
        MOI.EqualTo{Float64},
    ) == 1

    x1_new = reference_map[x[1]]
    @test owner_model(x1_new) === new_model
    @test "x[1]" == @inferred name(x1_new)

    cref_2_new = reference_map[cref[2]]
    @test cref_2_new.model === new_model
    @test "cref[2]" == @inferred name(cref_2_new)
end

function test_copy_filter_denseaxisarray()
    model = Model()
    @variable(model, x[i = 1:2], container = DenseAxisArray)
    @constraint(model, cref[i = 1:2], x[i] == 1, container = DenseAxisArray)
    @test num_constraints(
        model,
        GenericAffExpr{Float64,VariableRef},
        MOI.EqualTo{Float64},
    ) == 2

    filter_constraints = (cr) -> cr != cref[1]
    new_model, reference_map =
        copy_model(model; filter_constraints = filter_constraints)
    @test num_constraints(
        new_model,
        GenericAffExpr{Float64,VariableRef},
        MOI.EqualTo{Float64},
    ) == 1

    x1_new = reference_map[x[1]]
    @test owner_model(x1_new) === new_model
    @test "x[1]" == @inferred name(x1_new)

    cref_2_new = reference_map[cref[2]]
    @test cref_2_new.model === new_model
    @test "cref[2]" == @inferred name(cref_2_new)
end

function test_copy_filter_sparseaxisarray()
    model = Model()
    @variable(model, x[i = 1:2], container = SparseAxisArray)
    @constraint(model, cref[i = 1:2], x[i] == 1, container = SparseAxisArray)
    @test num_constraints(
        model,
        GenericAffExpr{Float64,VariableRef},
        MOI.EqualTo{Float64},
    ) == 2

    filter_constraints = (cr) -> cr != cref[1]
    new_model, reference_map =
        copy_model(model; filter_constraints = filter_constraints)
    @test num_constraints(
        new_model,
        GenericAffExpr{Float64,VariableRef},
        MOI.EqualTo{Float64},
    ) == 1

    x1_new = reference_map[x[1]]
    @test owner_model(x1_new) === new_model
    @test "x[1]" == @inferred name(x1_new)

    cref_2_new = reference_map[cref[2]]
    @test cref_2_new.model === new_model
    @test "cref[2]" == @inferred name(cref_2_new)
end

function test_copy_conflict()
    model = Model()
    @variable(model, x[i = 1:2], container = SparseAxisArray)
    @constraint(model, cref[i = 1:2], x[i] == 1, container = SparseAxisArray)
    @test num_constraints(
        model,
        GenericAffExpr{Float64,VariableRef},
        MOI.EqualTo{Float64},
    ) == 2

    set_optimizer(
        model,
        () -> MOIU.MockOptimizer(
            MOIU.Model{Float64}();
            eval_objective_value = false,
        ),
    )
    optimize!(model)

    mockoptimizer = backend(model).optimizer.model
    MOI.set(mockoptimizer, MOI.TerminationStatus(), INFEASIBLE)
    MOI.set(mockoptimizer, MOI.ConflictStatus(), MOI.CONFLICT_FOUND)
    MOI.set(
        mockoptimizer,
        MOI.ConstraintConflictStatus(),
        optimizer_index(cref[1]),
        MOI.IN_CONFLICT,
    )
    MOI.set(
        mockoptimizer,
        MOI.ConstraintConflictStatus(),
        optimizer_index(cref[2]),
        MOI.NOT_IN_CONFLICT,
    )

    new_model, reference_map = copy_conflict(model)
    @test num_constraints(
        new_model,
        GenericAffExpr{Float64,VariableRef},
        MOI.EqualTo{Float64},
    ) == 1

    x1_new = reference_map[x[1]]
    @test owner_model(x1_new) === new_model
    @test "x[1]" == @inferred name(x1_new)

    cref_1_new = reference_map[cref[1]]
    @test cref_1_new.model === new_model
    @test "cref[1]" == @inferred name(cref_1_new)
end

function test_nlp_data_error()
    model = Model()
    err = ErrorException(
        "The internal field `.nlp_data` was removed from `Model` in JuMP " *
        "v.1.2.0. If you encountered this message without going " *
        "`model.nlp_data`, it means you are using a package that is " *
        "incompatible with your installed version of JuMP. As a " *
        "temporary fix, install a compatible version with " *
        "`import Pkg; Pkg.pkg\"add JuMP@1.1\"`, then restart Julia for " *
        "the changes to take effect. In addition, you should open a " *
        "GitHub issue for the package you are using so that the issue " *
        "can be fixed for future users.",
    )
    @test_throws(err, model.nlp_data)
    return
end

function test_reset_optimizer()
    direct = direct_model(
        MOI.Utilities.MockOptimizer(MOI.Utilities.Model{Float64}()),
    )
    @test_throws(
        ErrorException(
            "The `$(MOI.Utilities.reset_optimizer)` function is not " *
            "supported in DIRECT mode.",
        ),
        MOI.Utilities.reset_optimizer(direct),
    )
    inner = MOI.Utilities.Model{Float64}()
    model = Model(() -> MOI.Utilities.MockOptimizer(inner))
    @variable(model, x >= 0)
    @objective(model, Min, x)
    @test MOI.is_empty(inner)
    MOI.Utilities.attach_optimizer(model)
    @test !MOI.is_empty(inner)
    MOI.Utilities.reset_optimizer(model)
    @test MOI.is_empty(inner)
    return
end

function test_drop_optimizer()
    direct = direct_model(
        MOI.Utilities.MockOptimizer(MOI.Utilities.Model{Float64}()),
    )
    @test_throws(
        ErrorException(
            "The `$(MOI.Utilities.drop_optimizer)` function is not supported " *
            "in DIRECT mode.",
        ),
        MOI.Utilities.drop_optimizer(direct),
    )
    model = Model() do
        return MOI.Utilities.MockOptimizer(MOI.Utilities.Model{Float64}())
    end
    @test solver_name(model) == "Mock"
    MOI.Utilities.drop_optimizer(model)
    @test solver_name(model) == "No optimizer attached."
    return
end

function test_optimizer_attribute_get_set()
    opt = optimizer_with_attributes() do
        return MOI.Utilities.MockOptimizer(MOI.Utilities.Model{Float64}())
    end
    @test get_optimizer_attribute(opt, MOI.Silent()) === nothing
    set_optimizer_attribute(opt, MOI.Silent(), true)
    @test get_optimizer_attribute(opt, MOI.Silent()) == true
    attr_a = MOI.RawOptimizerAttribute("a")
    @test get_optimizer_attribute(opt, attr_a) === nothing
    @test get_optimizer_attribute(opt, "a") === nothing
    set_optimizer_attribute(opt, attr_a, 1.0)
    @test get_optimizer_attribute(opt, attr_a) == 1.0
    @test get_optimizer_attribute(opt, "a") == 1.0
    set_optimizer_attribute(opt, "a", 2.0)
    @test get_optimizer_attribute(opt, attr_a) == 2.0
    @test get_optimizer_attribute(opt, "a") == 2.0
    return
end

function test_optimize_not_called_warning()
    model = Model() do
        return MOI.Utilities.MockOptimizer(
            MOI.Utilities.UniversalFallback(MOI.Utilities.Model{Float64}()),
        )
    end
    @variable(model, x >= 0)
    @constraint(model, c, 2x <= 1)
    optimize!(model)
    set_start_value(x, 0.0)
    @test_logs (:warn,) (@test_throws OptimizeNotCalled objective_value(model))
    @test_logs (:warn,) (@test_throws OptimizeNotCalled value(x))
    @test_logs (:warn,) (@test_throws OptimizeNotCalled value(c))
    return
end

function test_get_set_attribute_model()
    model = Model()
    @test get_attribute(model, MOI.Name()) == ""
    set_attribute(model, MOI.Name(), "Test")
    @test get_attribute(model, MOI.Name()) == "Test"
    set_attributes(model, MOI.Name() => "Test2")
    @test get_attribute(model, MOI.Name()) == "Test2"
    return
end

function test_get_set_attribute_optimizer()
    model = Model() do
        return MOI.Utilities.MockOptimizer(
            MOI.Utilities.UniversalFallback(MOI.Utilities.Model{Float64}()),
        )
    end
    set_attribute(model, MOI.Silent(), true)
    @test get_attribute(model, MOI.Silent()) == true
    set_attributes(model, MOI.Silent() => false)
    @test get_attribute(model, MOI.Silent()) == false
    return
end

function test_get_set_attribute_optimizer_with_attributes()
    optimizer = optimizer_with_attributes() do
        return MOI.Utilities.MockOptimizer(
            MOI.Utilities.UniversalFallback(MOI.Utilities.Model{Float64}()),
        )
    end
    set_attribute(optimizer, MOI.Silent(), true)
    @test get_attribute(optimizer, MOI.Silent()) == true
    set_attributes(optimizer, MOI.Silent() => false)
    @test get_attribute(optimizer, MOI.Silent()) == false
    return
end

function test_get_set_attribute_variable()
    model = Model()
    @variable(model, x)
    @test get_attribute(x, MOI.VariableName()) == "x"
    set_attribute(x, MOI.VariableName(), "y")
    @test get_attribute(x, MOI.VariableName()) == "y"
    set_attributes(x, MOI.VariableName() => "x")
    @test get_attribute(x, MOI.VariableName()) == "x"
    return
end

function test_get_set_attribute_constraint()
    model = Model()
    @variable(model, x)
    @constraint(model, c, 2 * x <= 1)
    @test get_attribute(c, MOI.ConstraintName()) == "c"
    set_attribute(c, MOI.ConstraintName(), "y")
    @test get_attribute(c, MOI.ConstraintName()) == "y"
    set_attributes(c, MOI.ConstraintName() => "c")
    @test get_attribute(c, MOI.ConstraintName()) == "c"
    return
end

function test_set_start_values()
    model = Model()
    @variable(model, x >= 0, Int)
    @constraint(model, c, 2 * x <= 1)
    @NLconstraint(model, nl, x^2 <= 1)
    @test_throws OptimizeNotCalled set_start_values(model)
    set_start_values(
        model;
        variable_primal_start = x -> Float64(index(x).value) + 0.1,
        constraint_primal_start = c -> Float64(index(c).value) + 0.2,
        constraint_dual_start = c -> Float64(index(c).value) + 0.3,
        nonlinear_dual_start = m -> fill(0.4, num_nonlinear_constraints(m)),
    )
    @test start_value(x) == index(x).value + 0.1
    @test start_value(c) == index(c).value + 0.2
    @test dual_start_value(c) == index(c).value + 0.3
    @test nonlinear_dual_start_value(model) == [0.4]
    # Only dual start
    model = Model()
    @variable(model, x >= 0, Int)
    @constraint(model, c, 2 * x <= 1)
    @NLconstraint(model, nl, x^2 <= 1)
    set_start_values(
        model;
        variable_primal_start = nothing,
        constraint_primal_start = nothing,
        constraint_dual_start = c -> Float64(index(c).value) + 0.5,
        nonlinear_dual_start = m -> fill(0.6, num_nonlinear_constraints(m)),
    )
    @test start_value(x) === nothing
    @test start_value(c) === nothing
    @test dual_start_value(c) == index(c).value + 0.5
    @test nonlinear_dual_start_value(model) == [0.6]
    # Only primal start
    model = Model()
    @variable(model, x >= 0, Int)
    @constraint(model, c, 2 * x <= 1)
    @NLconstraint(model, nl, x^2 <= 1)
    set_start_values(
        model;
        variable_primal_start = x -> Float64(index(x).value) + 0.7,
        constraint_primal_start = c -> Float64(index(c).value) + 0.8,
        constraint_dual_start = nothing,
        nonlinear_dual_start = nothing,
    )
    @test start_value(x) === index(x).value + 0.7
    @test start_value(c) === index(c).value + 0.8
    @test dual_start_value(c) === nothing
    @test nonlinear_dual_start_value(model) === nothing
    return
end

function test_model_quad_to_soc_start_values()
    inner = MOI.Utilities.MockOptimizer(
        MOI.Bridges.Constraint.QuadtoSOC{Float64}(
            MOI.Utilities.UniversalFallback(MOI.Utilities.Model{Float64}()),
        ),
    )
    model = direct_model(inner)
    @variable(model, x)
    @constraint(model, c, x^2 <= 0)
    @test_throws ErrorException set_start_value(c, 0.0)
    set_start_values(
        model;
        variable_primal_start = x -> 1.0,
        constraint_primal_start = x -> 1.0,
        constraint_dual_start = x -> 1.0,
    )
    @test start_value(c) == 1.0
    return
end

function test_keyword_getindex()
    err = JuMP._get_index_keyword_indexing_error()
    model = Model()
    @variable(model, x[i = 1:2])
    # Vector{VariableRef}
    @test_throws err x[i = 1]
    @test_throws err x[i = 2]
    @test_throws err x[i = :]
    @test_throws err x[i = 1:2]
    @test_throws BoundsError x[]
    # Matrix{VariableRef}
    @variable(model, y[i = 1:2, j = 1:2])
    @test_throws err y[i = 1]
    @test_throws err y[i = 2, j = 2]
    @test_throws BoundsError y[]
    # Vector{AffExpr}
    @expression(model, ex[i = 1:2], x[i] + i)
    @test_throws err ex[i = 1]
    @test_throws err ex[i = 2, j = 2]
    @test_throws BoundsError ex[]
    # Vector{ConstraintRef}
    @constraint(model, c[i = 1:2], x[i] + i <= 1)
    @test_throws err c[i = 1]
    @test_throws err c[i = 1:2]
    @test_throws err c[i = :]
    @test_throws BoundsError c[]
    return
end

function test_getindex_no_arg()
    model = Model()
    @variable(model, x)
    y = [x]
    @test y[] === x
    ex = x + 1
    y = reshape([ex], 1, 1, 1)
    @test y[] === ex
    c = @constraint(model, x >= 1)
    y = [c]
    @test y[] === c
    return
end

end  # module TestModels

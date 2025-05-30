#  Copyright 2017, Iain Dunning, Joey Huchette, Miles Lubin, and contributors
#  This Source Code Form is subject to the terms of the Mozilla Public
#  License, v. 2.0. If a copy of the MPL was not distributed with this
#  file, You can obtain one at https://mozilla.org/MPL/2.0/.

# The tests here check JuMP's model generation and communication with solvers.
# Model generation is checked by comparing the internal model with a serialized
# test model (in MOI.Utilities's lightweight text format).
# Communication with solvers is tested by using a mock solver with solution data
# that we feed to it. Prior to using this testing approach, we would test JuMP
# by calling real solvers, which was flakey and slow.

# Note: No attempt is made to use correct solution data. We're only testing
# that the plumbing works. This could change if JuMP gains the ability to verify
# feasibility independently of a solver.

module TestGenerateAndSolve

using LinearAlgebra
using JuMP
using Test

function test_generate_solve_LP()
    m = Model()
    @variable(m, x <= 2.0)
    @variable(m, y >= 0.0)
    @objective(m, Min, -x)
    c = @constraint(m, x + y <= 1)
    set_name(c, "c")
    modelstring = """
    variables: x, y
    minobjective: -1.0*x
    x <= 2.0
    y >= 0.0
    c: x + y <= 1.0
    """
    model = MOI.Utilities.Model{Float64}()
    MOI.Utilities.loadfromstring!(model, modelstring)
    MOI.Test.util_test_models_equal(
        backend(m).model_cache,
        model,
        ["x", "y"],
        ["c"],
        [("x", MOI.LessThan(2.0)), ("y", MOI.GreaterThan(0.0))],
    )
    set_optimizer(
        m,
        () -> MOI.Utilities.MockOptimizer(
            MOI.Utilities.Model{Float64}();
            eval_objective_value = false,
        ),
    )
    optimize!(m)
    mock = unsafe_backend(m)
    MOI.set(mock, MOI.TerminationStatus(), MOI.OPTIMAL)
    MOI.set(mock, MOI.RawStatusString(), "solver specific string")
    MOI.set(mock, MOI.ObjectiveValue(), -1.0)
    MOI.set(mock, MOI.ResultCount(), 1)
    MOI.set(mock, MOI.PrimalStatus(), MOI.FEASIBLE_POINT)
    MOI.set(mock, MOI.DualStatus(), MOI.FEASIBLE_POINT)
    MOI.set(mock, MOI.VariablePrimal(), optimizer_index(x), 1.0)
    MOI.set(mock, MOI.VariablePrimal(), optimizer_index(y), 0.0)
    MOI.set(mock, MOI.ConstraintDual(), optimizer_index(c), -1.0)
    MOI.set(mock, MOI.ConstraintDual(), optimizer_index(UpperBoundRef(x)), 0.0)
    MOI.set(mock, MOI.ConstraintDual(), optimizer_index(LowerBoundRef(y)), 1.0)
    MOI.set(mock, MOI.SimplexIterations(), Int64(1))
    MOI.set(mock, MOI.BarrierIterations(), Int64(1))
    MOI.set(mock, MOI.NodeCount(), Int64(1))
    @test has_values(m)
    @test MOI.OPTIMAL == @inferred termination_status(m)
    @test "solver specific string" == raw_status(m)
    @test MOI.FEASIBLE_POINT == @inferred primal_status(m)
    @test 1.0 == @inferred value(x)
    @test 0.0 == @inferred value(y)
    @test 1.0 == @inferred value(x + y)
    @test 1.0 == @inferred value(c)
    @test -1.0 == objective_value(m)
    @test -1.0 == @inferred dual_objective_value(m)
    @test MOI.FEASIBLE_POINT == @inferred dual_status(m)
    @test -1.0 == @inferred dual(c)
    @test 0.0 == @inferred dual(UpperBoundRef(x))
    @test 1.0 == @inferred dual(LowerBoundRef(y))
    @test 1 == simplex_iterations(m)
    @test 1 == barrier_iterations(m)
    @test 1 == node_count(m)
    return
end

function test_generate_solve_LP_direct_mode()
    mock = MOI.Utilities.MockOptimizer(
        MOI.Utilities.Model{Float64}();
        eval_objective_value = false,
    )
    m = direct_model(mock)
    @variable(m, x <= 2.0)
    @variable(m, y >= 0.0)
    @objective(m, Min, -x)
    c = @constraint(m, x + y <= 1)
    MOI.set(mock, MOI.TerminationStatus(), MOI.OPTIMAL)
    MOI.set(mock, MOI.RawStatusString(), "solver specific string")
    MOI.set(mock, MOI.ObjectiveValue(), -1.0)
    MOI.set(mock, MOI.ResultCount(), 1)
    MOI.set(mock, MOI.PrimalStatus(), MOI.FEASIBLE_POINT)
    MOI.set(mock, MOI.DualStatus(), MOI.FEASIBLE_POINT)
    MOI.set(mock, MOI.VariablePrimal(), optimizer_index(x), 1.0)
    MOI.set(mock, MOI.VariablePrimal(), optimizer_index(y), 0.0)
    MOI.set(mock, MOI.ConstraintDual(), optimizer_index(c), -1.0)
    MOI.set(mock, MOI.ConstraintDual(), optimizer_index(UpperBoundRef(x)), 0.0)
    MOI.set(mock, MOI.ConstraintDual(), optimizer_index(LowerBoundRef(y)), 1.0)
    MOI.set(mock, MOI.SimplexIterations(), Int64(1))
    MOI.set(mock, MOI.BarrierIterations(), Int64(1))
    MOI.set(mock, MOI.NodeCount(), Int64(1))
    optimize!(m)
    @test has_values(m)
    @test MOI.OPTIMAL == @inferred termination_status(m)
    @test "solver specific string" == raw_status(m)
    @test MOI.FEASIBLE_POINT == @inferred primal_status(m)
    @test 1.0 == @inferred value(x)
    @test 0.0 == @inferred value(y)
    @test 1.0 == @inferred value(x + y)
    @test -1.0 == objective_value(m)
    @test MOI.FEASIBLE_POINT == @inferred dual_status(m)
    @test -1.0 == @inferred dual(c)
    @test 0.0 == @inferred dual(UpperBoundRef(x))
    @test 1.0 == @inferred dual(LowerBoundRef(y))
    @test 1 == simplex_iterations(m)
    @test 1 == barrier_iterations(m)
    @test 1 == node_count(m)
    return
end

function test_generate_solve_IP()
    m = Model(
        () -> MOI.Utilities.MockOptimizer(
            MOI.Utilities.Model{Float64}();
            eval_objective_value = false,
        ),
    )
    @variable(m, x == 1.0, Int)
    @variable(m, y, Bin)
    @objective(m, Max, x)
    modelstring = """
    variables: x, y
    maxobjective: x
    x == 1.0
    x in Integer()
    y in ZeroOne()
    """
    model = MOI.Utilities.Model{Float64}()
    MOI.Utilities.loadfromstring!(model, modelstring)
    MOI.Test.util_test_models_equal(
        backend(m).model_cache,
        model,
        ["x", "y"],
        String[],
        [("x", MOI.EqualTo(1.0)), ("x", MOI.Integer()), ("y", MOI.ZeroOne())],
    )
    MOI.Utilities.attach_optimizer(m)
    mock = unsafe_backend(m)
    MOI.set(mock, MOI.TerminationStatus(), MOI.OPTIMAL)
    MOI.set(mock, MOI.RawStatusString(), "solver specific string")
    MOI.set(mock, MOI.ObjectiveValue(), 1.0)
    MOI.set(mock, MOI.ResultCount(), 1)
    MOI.set(mock, MOI.PrimalStatus(), MOI.FEASIBLE_POINT)
    MOI.set(mock, MOI.VariablePrimal(), optimizer_index(x), 1.0)
    MOI.set(mock, MOI.VariablePrimal(), optimizer_index(y), 0.0)
    MOI.set(mock, MOI.DualStatus(), MOI.NO_SOLUTION)
    MOI.set(mock, MOI.SimplexIterations(), Int64(1))
    MOI.set(mock, MOI.BarrierIterations(), Int64(1))
    MOI.set(mock, MOI.NodeCount(), Int64(1))
    MOI.set(mock, MOI.RelativeGap(), 0.0)
    optimize!(m)
    @test has_values(m)
    @test MOI.OPTIMAL == @inferred termination_status(m)
    @test "solver specific string" == raw_status(m)
    @test MOI.FEASIBLE_POINT == @inferred primal_status(m)
    @test 1.0 == @inferred value(x)
    @test 0.0 == @inferred value(y)
    @test 1.0 == objective_value(m)
    @test 1 == simplex_iterations(m)
    @test 1 == barrier_iterations(m)
    @test 1 == node_count(m)
    @test 0.0 == @inferred relative_gap(m)
    @test !has_duals(m)
    return
end

function test_generate_solve_QCQP()
    m = Model()
    @variable(m, x)
    @variable(m, y)
    @objective(m, Min, x^2)
    @constraint(m, c1, 2x * y <= 1)
    @constraint(m, c2, y^2 == x^2)
    @constraint(m, c3, 2x + 3y * x >= 2)
    modelstring = """
    variables: x, y
    minobjective: 1*x*x
    c1: 2*x*y <= 1.0
    c2: 1*y*y + -1*x*x == 0.0
    c3: 2x + 3*y*x >= 2.0
    """
    model = MOI.Utilities.Model{Float64}()
    MOI.Utilities.loadfromstring!(model, modelstring)
    MOI.Test.util_test_models_equal(
        backend(m).model_cache,
        model,
        ["x", "y"],
        ["c1", "c2", "c3"],
    )
    set_optimizer(
        m,
        () -> MOI.Utilities.MockOptimizer(
            MOI.Utilities.Model{Float64}();
            eval_objective_value = false,
        ),
    )
    optimize!(m)
    mock = unsafe_backend(m)
    MOI.set(mock, MOI.TerminationStatus(), MOI.OPTIMAL)
    MOI.set(mock, MOI.RawStatusString(), "solver specific string")
    MOI.set(mock, MOI.ObjectiveValue(), -1.0)
    MOI.set(mock, MOI.ResultCount(), 1)
    MOI.set(mock, MOI.PrimalStatus(), MOI.FEASIBLE_POINT)
    MOI.set(mock, MOI.DualStatus(), MOI.FEASIBLE_POINT)
    MOI.set(mock, MOI.VariablePrimal(), optimizer_index(x), 1.0)
    MOI.set(mock, MOI.VariablePrimal(), optimizer_index(y), 0.0)
    MOI.set(mock, MOI.ConstraintDual(), optimizer_index(c1), -1.0)
    MOI.set(mock, MOI.ConstraintDual(), optimizer_index(c2), 2.0)
    MOI.set(mock, MOI.ConstraintDual(), optimizer_index(c3), 3.0)
    MOI.set(mock, MOI.SimplexIterations(), Int64(1))
    MOI.set(mock, MOI.BarrierIterations(), Int64(1))
    MOI.set(mock, MOI.NodeCount(), Int64(1))
    @test has_values(m)
    @test MOI.OPTIMAL == @inferred termination_status(m)
    @test "solver specific string" == raw_status(m)
    @test MOI.FEASIBLE_POINT == @inferred primal_status(m)
    @test 1.0 == @inferred value(x)
    @test 0.0 == @inferred value(y)
    @test -1.0 == objective_value(m)
    @test 5.0 == @inferred dual_objective_value(m)
    @test MOI.FEASIBLE_POINT == @inferred dual_status(m)
    @test -1.0 == @inferred dual(c1)
    @test 2.0 == @inferred dual(c2)
    @test 3.0 == @inferred dual(c3)
    @test 2.0 == @inferred value(2 * x + 3 * y * x)
    @test 1 == simplex_iterations(m)
    @test 1 == barrier_iterations(m)
    @test 1 == node_count(m)
    return
end

function test_generate_solve_SOC()
    m = Model()
    @variables(m, begin
        x
        y
        z
    end)
    @objective(m, Max, 1.0 * x)
    @constraint(m, varsoc, [x, y, z] in SecondOrderCone())
    # Equivalent to `[x+y,z,1.0] in SecondOrderCone()`
    @constraint(m, affsoc, [x + y, z, 1.0] in MOI.SecondOrderCone(3))
    @constraint(m, rotsoc, [x + 1, y, z] in RotatedSecondOrderCone())
    modelstring = """
    variables: x, y, z
    maxobjective: 1.0*x
    varsoc: [x,y,z] in SecondOrderCone(3)
    affsoc: [x+y,z,1.0] in SecondOrderCone(3)
    rotsoc: [x+1,y,z] in RotatedSecondOrderCone(3)
    """
    model = MOI.Utilities.Model{Float64}()
    MOI.Utilities.loadfromstring!(model, modelstring)
    MOI.Test.util_test_models_equal(
        backend(m).model_cache,
        model,
        ["x", "y", "z"],
        ["varsoc", "affsoc", "rotsoc"],
    )
    mock = MOI.Utilities.MockOptimizer(
        MOI.Utilities.Model{Float64}();
        eval_objective_value = false,
        eval_variable_constraint_dual = false,
    )
    MOI.Utilities.reset_optimizer(m, mock)
    MOI.Utilities.attach_optimizer(m)
    MOI.set(mock, MOI.TerminationStatus(), MOI.OPTIMAL)
    MOI.set(mock, MOI.RawStatusString(), "solver specific string")
    MOI.set(mock, MOI.ResultCount(), 1)
    MOI.set(mock, MOI.PrimalStatus(), MOI.FEASIBLE_POINT)
    MOI.set(mock, MOI.DualStatus(), MOI.FEASIBLE_POINT)
    MOI.set(mock, MOI.VariablePrimal(), optimizer_index(x), 1.0)
    MOI.set(mock, MOI.VariablePrimal(), optimizer_index(y), 0.0)
    MOI.set(mock, MOI.VariablePrimal(), optimizer_index(z), 0.0)
    MOI.set(
        mock,
        MOI.ConstraintDual(),
        optimizer_index(varsoc),
        [-1.0, -2.0, -3.0],
    )
    MOI.set(
        mock,
        MOI.ConstraintDual(),
        optimizer_index(affsoc),
        [1.0, 2.0, 3.0],
    )
    MOI.set(mock, MOI.SimplexIterations(), Int64(1))
    MOI.set(mock, MOI.BarrierIterations(), Int64(1))
    MOI.set(mock, MOI.NodeCount(), Int64(1))
    optimize!(m)
    @test has_values(m)
    @test MOI.OPTIMAL == @inferred termination_status(m)
    @test "solver specific string" == raw_status(m)
    @test MOI.FEASIBLE_POINT == @inferred primal_status(m)
    @test 1.0 == @inferred value(x)
    @test 0.0 == @inferred value(y)
    @test 0.0 == @inferred value(z)
    @test has_duals(m)
    @test [-1.0, -2.0, -3.0] == @inferred dual(varsoc)
    @test [1.0, 2.0, 3.0] == @inferred dual(affsoc)
    @test 1 == simplex_iterations(m)
    @test 1 == barrier_iterations(m)
    @test 1 == node_count(m)
    return
end

function test_generate_solve_SDP()
    m = Model()
    @variable(m, x[1:2, 1:2], Symmetric)
    set_name(x[1, 1], "x11")
    set_name(x[1, 2], "x12")
    set_name(x[2, 2], "x22")
    @objective(m, Max, tr(x))
    var_psd = @constraint(m, x in PSDCone())
    set_name(var_psd, "var_psd")
    sym_psd = @constraint(m, Symmetric(x - [1.0 0.0; 0.0 1.0]) in PSDCone())
    set_name(sym_psd, "sym_psd")
    con_psd = @constraint(m, x >= [1.0 0.0; 0.0 1.0], PSDCone())
    set_name(con_psd, "con_psd")
    modelstring = """
    variables: x11, x12, x22
    maxobjective: 1.0*x11 + 1.0*x22
    var_psd: [x11,x12,x22] in PositiveSemidefiniteConeTriangle(2)
    sym_psd: [x11 + -1.0,x12,x22 + -1.0] in PositiveSemidefiniteConeTriangle(2)
    con_psd: [x11 + -1.0,x12,x12,x22 + -1.0] in PositiveSemidefiniteConeSquare(2)
    """
    model = MOI.Utilities.Model{Float64}()
    MOI.Utilities.loadfromstring!(model, modelstring)
    MOI.Test.util_test_models_equal(
        backend(m).model_cache,
        model,
        ["x11", "x12", "x22"],
        ["var_psd", "sym_psd", "con_psd"],
    )
    mock = MOI.Utilities.MockOptimizer(
        MOI.Utilities.Model{Float64}();
        eval_objective_value = false,
        eval_variable_constraint_dual = false,
    )
    MOI.Utilities.reset_optimizer(m, mock)
    MOI.Utilities.attach_optimizer(m)
    MOI.set(mock, MOI.TerminationStatus(), MOI.OPTIMAL)
    MOI.set(mock, MOI.RawStatusString(), "solver specific string")
    MOI.set(mock, MOI.ResultCount(), 1)
    MOI.set(mock, MOI.PrimalStatus(), MOI.FEASIBLE_POINT)
    MOI.set(mock, MOI.DualStatus(), MOI.FEASIBLE_POINT)
    MOI.set(mock, MOI.VariablePrimal(), optimizer_index(x[1, 1]), 1.0)
    MOI.set(mock, MOI.VariablePrimal(), optimizer_index(x[1, 2]), 2.0)
    MOI.set(mock, MOI.VariablePrimal(), optimizer_index(x[2, 2]), 4.0)
    MOI.set(
        mock,
        MOI.ConstraintDual(),
        optimizer_index(var_psd),
        [1.0, 2.0, 3.0],
    )
    MOI.set(
        mock,
        MOI.ConstraintDual(),
        optimizer_index(sym_psd),
        [4.0, 5.0, 6.0],
    )
    MOI.set(
        mock,
        MOI.ConstraintDual(),
        optimizer_index(con_psd),
        [7.0, 8.0, 9.0, 10.0],
    )
    MOI.set(mock, MOI.SimplexIterations(), Int64(1))
    MOI.set(mock, MOI.BarrierIterations(), Int64(1))
    MOI.set(mock, MOI.NodeCount(), Int64(1))
    optimize!(m)
    @test MOI.OPTIMAL == @inferred termination_status(m)
    @test "solver specific string" == raw_status(m)
    @test MOI.FEASIBLE_POINT == @inferred primal_status(m)
    @test has_values(m)
    @test [1.0 2.0; 2.0 4.0] == value.(x)
    @test value(x) isa Symmetric
    @test [1.0 2.0; 2.0 4.0] == @inferred value(x)
    @test value(var_psd) isa Symmetric
    @test [1.0 2.0; 2.0 4.0] == @inferred value(var_psd)
    @test value(sym_psd) isa Symmetric
    @test [0.0 2.0; 2.0 3.0] == @inferred value(sym_psd)
    @test value(con_psd) isa Matrix
    @test [0.0 2.0; 2.0 3.0] == @inferred value(con_psd)
    @test has_duals(m)
    @test dual(var_psd) isa Symmetric
    @test [1.0 2.0; 2.0 3.0] == @inferred dual(var_psd)
    @test dual(sym_psd) isa Symmetric
    @test [4.0 5.0; 5.0 6.0] == @inferred dual(sym_psd)
    @test dual(con_psd) isa Matrix
    @test [7.0 9.0; 8.0 10.0] == @inferred dual(con_psd)
    @test 1 == simplex_iterations(m)
    @test 1 == barrier_iterations(m)
    @test 1 == node_count(m)
    return
end

function test_generate_solve_unsupported_nonlinear_problems()
    model =
        Model(() -> MOI.Utilities.MockOptimizer(MOI.Utilities.Model{Float64}()))
    # We put some additional constraints here to confuse
    # JuMP._uses_new_nonlinear_interface.
    @variable(model, x >= 0)
    @constraint(model, 2 * x <= 1)
    @NLobjective(model, Min, sin(x))
    err = ErrorException(
        "The solver does not support nonlinear problems " *
        "(that is, NLobjective and NLconstraint).",
    )
    @test_throws err optimize!(model)
    return
end

function test_generate_solve_lower_bound_already_set()
    model = Model() do
        return MOI.Bridges.Constraint.Scalarize{Float64}(
            MOI.Utilities.MockOptimizer(MOI.Utilities.Model{Float64}()),
        )
    end
    @variable(model, x >= -1)
    @constraint(model, [x] in Zeros())
    @test_throws(
        ErrorException(
            """
            The model contains a variable for which multiple lower bounds have been
            set. Each variable can have at most one lower bound.

            This error can occur if you have added a bounded variable in `@variable`
            and then added a constraint on the variable using `@constraint`.

            To fix, change the `@constraint` call by replacing `x` with `1.0 * x`.

            Here is an example:
            ```julia
            model = Model()
            @variable(model, x >= 0)
            @constraint(model, [x] in Nonnegatives())
            # Rewrite as
            @constraint(model, [1.0 * x] in Nonnegatives())
            ```
            """,
        ),
        optimize!(model),
    )
    return
end

function test_generate_solve_upper_bound_already_set()
    model = Model() do
        return MOI.Bridges.Constraint.Scalarize{Float64}(
            MOI.Utilities.MockOptimizer(MOI.Utilities.Model{Float64}()),
        )
    end
    @variable(model, x <= 1)
    @constraint(model, [x] in Zeros())
    @test_throws(
        ErrorException(
            """
            The model contains a variable for which multiple upper bounds have been
            set. Each variable can have at most one upper bound.

            This error can occur if you have added a bounded variable in `@variable`
            and then added a constraint on the variable using `@constraint`.

            To fix, change the `@constraint` call by replacing `x` with `1.0 * x`.

            Here is an example:
            ```julia
            model = Model()
            @variable(model, x <= 1)
            @constraint(model, [x] in Nonpositives())
            # Rewrite as
            @constraint(model, [1.0 * x] in Nonpositives())
            ```
            """,
        ),
        optimize!(model),
    )
    return
end

function test_generate_solve_ResultCount()
    m = Model()
    @variable(m, x >= 0.0)
    @variable(m, y >= 0.0)
    @objective(m, Max, x + y)
    @constraint(m, c1, x <= 2)
    @constraint(m, c2, x + y <= 1)
    model = MOI.Utilities.Model{Float64}()
    MOI.Utilities.loadfromstring!(
        model,
        """
variables: x, y
maxobjective: x + y
x >= 0.0
y >= 0.0
x <= 2.0
c2: x + y <= 1.0
""",
    )
    set_optimizer(
        m,
        () -> MOI.Utilities.MockOptimizer(
            MOI.Utilities.Model{Float64}();
            eval_objective_value = false,
        ),
    )
    optimize!(m)
    mock = unsafe_backend(m)
    MOI.set(mock, MOI.TerminationStatus(), MOI.OPTIMAL)
    MOI.set(mock, MOI.ResultCount(), 2)
    aff_expr = @expression(m, x + y)
    quad_expr = @expression(m, x * y)
    nl_expr = @NLexpression(m, log(x + y))
    @test result_count(m) == 2
    MOI.set(mock, MOI.PrimalStatus(1), MOI.FEASIBLE_POINT)
    MOI.set(mock, MOI.DualStatus(1), MOI.FEASIBLE_POINT)
    MOI.set(mock, MOI.ObjectiveValue(1), 1.0)
    MOI.set(mock, MOI.DualObjectiveValue(1), 1.0)
    MOI.set(mock, MOI.VariablePrimal(1), optimizer_index(x), 1.0)
    MOI.set(mock, MOI.VariablePrimal(1), optimizer_index(y), 0.0)
    MOI.set(mock, MOI.ConstraintDual(1), optimizer_index(c1), 0.0)
    MOI.set(mock, MOI.ConstraintDual(1), optimizer_index(c2), -1.0)
    @test MOI.OPTIMAL == @inferred termination_status(m)
    @test MOI.FEASIBLE_POINT == @inferred primal_status(m, result = 1)
    @test MOI.FEASIBLE_POINT == @inferred dual_status(m, result = 1)
    @test 1.0 == objective_value(m; result = 1)
    @test 1.0 == @inferred dual_objective_value(m, result = 1)
    @test 1.0 == @inferred value(x, result = 1)
    @test 0.0 == @inferred value(y, result = 1)
    @test 1.0 == @inferred value(aff_expr, result = 1)
    @test 0.0 == @inferred value(quad_expr, result = 1)
    @test 0.0 == @inferred value(nl_expr, result = 1)
    @test 1.0 == @inferred value(c2, result = 1)
    @test 0.0 == @inferred dual(c1, result = 1)
    @test -1.0 == @inferred dual(c2, result = 1)
    @test 0.0 == @inferred dual(LowerBoundRef(x), result = 1)
    @test 0.0 == @inferred dual(LowerBoundRef(y), result = 1)
    MOI.set(mock, MOI.PrimalStatus(2), MOI.FEASIBLE_POINT)
    MOI.set(mock, MOI.DualStatus(2), MOI.FEASIBLE_POINT)
    MOI.set(mock, MOI.ObjectiveValue(2), 1.0)
    MOI.set(mock, MOI.DualObjectiveValue(2), 1.0)
    MOI.set(mock, MOI.VariablePrimal(2), optimizer_index(x), 0.0)
    MOI.set(mock, MOI.VariablePrimal(2), optimizer_index(y), 1.0)
    MOI.set(mock, MOI.ConstraintDual(2), optimizer_index(c1), 0.0)
    MOI.set(mock, MOI.ConstraintDual(2), optimizer_index(c2), -1.0)
    @test MOI.FEASIBLE_POINT == @inferred primal_status(m, result = 2)
    @test MOI.FEASIBLE_POINT == @inferred dual_status(m, result = 2)
    @test 1.0 == objective_value(m; result = 2)
    @test 1.0 == @inferred dual_objective_value(m, result = 2)
    @test 0.0 == @inferred value(x, result = 2)
    @test 1.0 == @inferred value(y, result = 2)
    @test 1.0 == @inferred value(aff_expr, result = 2)
    @test 0.0 == @inferred value(quad_expr, result = 2)
    @test 0.0 == @inferred value(nl_expr, result = 2)
    @test 1.0 == @inferred value(c2, result = 2)
    @test 0.0 == @inferred dual(c1, result = 2)
    @test -1.0 == @inferred dual(c2, result = 2)
    @test 0.0 == @inferred dual(LowerBoundRef(x), result = 2)
    @test 0.0 == @inferred dual(LowerBoundRef(y), result = 2)
    @test MOI.NO_SOLUTION == @inferred primal_status(m, result = 3)
    @test MOI.NO_SOLUTION == @inferred dual_status(m, result = 3)
    @test_throws MOI.ResultIndexBoundsError objective_value(m, result = 3)
    @test_throws MOI.ResultIndexBoundsError dual_objective_value(m, result = 3)
    @test_throws MOI.ResultIndexBoundsError value(x, result = 3)
    @test_throws MOI.ResultIndexBoundsError value(aff_expr, result = 3)
    @test_throws MOI.ResultIndexBoundsError value(quad_expr, result = 3)
    @test_throws MOI.ResultIndexBoundsError value(nl_expr, result = 3)
    @test_throws MOI.ResultIndexBoundsError value(c2, result = 3)
    @test_throws MOI.ResultIndexBoundsError dual(c1, result = 3)
    @test_throws MOI.ResultIndexBoundsError dual(c2, result = 3)
    @test_throws MOI.ResultIndexBoundsError dual(LowerBoundRef(x), result = 3)
    @test_throws MOI.ResultIndexBoundsError dual(LowerBoundRef(y), result = 3)
    return
end

function test_generate_solve_vector_objective()
    model = Model() do
        return MOI.Utilities.MockOptimizer(MOI.Utilities.Model{Float64}())
    end
    @variable(model, x >= 0.0)
    @variable(model, y >= 1.0)
    @objective(model, Min, [x, y])
    optimize!(model)
    mock = unsafe_backend(model)
    MOI.set(mock, MOI.TerminationStatus(), MOI.OPTIMAL)
    MOI.set(mock, MOI.RawStatusString(), "Optimal")
    MOI.set(mock, MOI.ResultCount(), 1)
    MOI.set(mock, MOI.PrimalStatus(), MOI.FEASIBLE_POINT)
    MOI.set(mock, MOI.DualStatus(), MOI.NO_SOLUTION)
    MOI.set(mock, MOI.ObjectiveBound(), [0.0, 1.0])
    MOI.set(mock, MOI.VariablePrimal(), optimizer_index(x), 0.0)
    MOI.set(mock, MOI.VariablePrimal(), optimizer_index(y), 1.0)
    @test sprint(print, solution_summary(model)) == """
    solution_summary(; result = 1, verbose = false)
    ├ solver_name          : Mock
    ├ Termination
    │ ├ termination_status : OPTIMAL
    │ ├ result_count       : 1
    │ ├ raw_status         : Optimal
    │ └ objective_bound    : [0.00000e+00,1.00000e+00]
    └ Solution (result = 1)
      ├ primal_status        : FEASIBLE_POINT
      ├ dual_status          : NO_SOLUTION
      └ objective_value      : [0.00000e+00,1.00000e+00]"""
    return
end

function test_optimize_not_called_direct()
    mock = MOI.Utilities.MockOptimizer(MOI.Utilities.Model{Float64}())
    model = direct_model(mock)
    @test_throws OptimizeNotCalled objective_value(model)
    return
end

function test_compute_conflict()
    mock = MOI.Utilities.MockOptimizer(MOI.Utilities.Model{Float64}())
    model = direct_model(mock)
    @test !mock.compute_conflict_called
    compute_conflict!(model)
    @test mock.compute_conflict_called
    return
end

function test__optimizer_for_precompile_tools()
    @test JuMP._optimizer_for_precompile_tools() isa MOI.Utilities.MockOptimizer
    return
end

function test_value_symmetric()
    inner = MOI.Utilities.MockOptimizer(
        MOI.Utilities.UniversalFallback(MOI.Utilities.Model{Float64}()),
    )
    model = direct_model(inner)
    @variable(model, x[i in 1:2, j in 1:2], Symmetric, start = i + j)
    Y = LinearAlgebra.Symmetric(1.0 .* x)
    Z = LinearAlgebra.I - Y #element Z.data[2, 1] is #undef
    @test_throws OptimizeNotCalled value(x)
    @test_throws OptimizeNotCalled value(Y)
    @test value(start_value, x) == LinearAlgebra.Symmetric([2 3; 3 4])
    @test value(start_value, Y) == LinearAlgebra.Symmetric([2 3; 3 4])
    @test value(start_value, Z) == LinearAlgebra.Symmetric([-1 -3; -3 -3])
    optimize!(model)
    MOI.set(inner, MOI.TerminationStatus(), MOI.OPTIMAL)
    MOI.set.(inner, MOI.VariablePrimal(), index.(x), [3 4; 4 5])
    @test value(x) == LinearAlgebra.Symmetric([3 4; 4 5])
    @test value(Y) == LinearAlgebra.Symmetric([3 4; 4 5])
    @test value(Z) == LinearAlgebra.Symmetric([-2 -4; -4 -4])
    @test_throws MOI.ResultIndexBoundsError value(x; result = 2)
    @test_throws MOI.ResultIndexBoundsError value(Y; result = 2)
    return
end

function test_value_hermitian()
    inner = MOI.Utilities.MockOptimizer(
        MOI.Utilities.UniversalFallback(MOI.Utilities.Model{Float64}()),
    )
    model = direct_model(inner)
    H = [1 2+3im; 2-3im 4]
    @variable(model, x[i in 1:2, j in 1:2], Hermitian, start = H[i, j])
    Z = LinearAlgebra.I - x
    @test_throws OptimizeNotCalled value(x)
    @test value(start_value, x) == LinearAlgebra.Hermitian(H)
    @test value(start_value, Z) == LinearAlgebra.Hermitian(LinearAlgebra.I - H)
    optimize!(model)
    MOI.set(inner, MOI.TerminationStatus(), MOI.OPTIMAL)
    y = index.(all_variables(model))
    MOI.set.(inner, MOI.VariablePrimal(), y, [3, 4, 5, 6])
    @test value(x) == LinearAlgebra.Hermitian([3 4+6im; 4-6im 5])
    @test value(Z) == LinearAlgebra.Hermitian([-2 -4-6im; 4+6im -4])
    @test_throws MOI.ResultIndexBoundsError value(x; result = 2)
    return
end

function test_value_array()
    inner = MOI.Utilities.MockOptimizer(
        MOI.Utilities.UniversalFallback(MOI.Utilities.Model{Float64}()),
    )
    model = direct_model(inner)
    @variable(model, x[i in 1:2], start = 3 - i)
    @variable(
        model,
        y[i in 1:2, j in 1:2] in ComplexPlane(),
        start = i + im * j
    )
    @test_throws OptimizeNotCalled value(x)
    @test_throws OptimizeNotCalled value(y)
    @test value(start_value, x) == [2, 1]
    @test value(start_value, y) == [1+im 1+2im; 2+im 2+2im]
    optimize!(model)
    MOI.set(inner, MOI.TerminationStatus(), MOI.OPTIMAL)
    idx = index.(all_variables(model))
    MOI.set.(inner, MOI.VariablePrimal(), idx, [1, 2, 3, 4, 5, 6, 7, 8, 9, 10])
    @test value(x) == [1, 2]
    @test value(y) == [3+4im 7+8im; 5+6im 9+10im]
    @test_throws MOI.ResultIndexBoundsError value(x; result = 2)
    @test_throws MOI.ResultIndexBoundsError value(y; result = 2)
    return
end

end  # module

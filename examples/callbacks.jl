#  Copyright 2017, Iain Dunning, Joey Huchette, Miles Lubin, and contributors
#  This Source Code Form is subject to the terms of the Mozilla Public
#  License, v. 2.0. If a copy of the MPL was not distributed with this
#  file, You can obtain one at http://mozilla.org/MPL/2.0/.
#############################################################################
# JuMP
# An algebraic modeling language for Julia
# See http://github.com/jump-dev/JuMP.jl
#############################################################################

using GLPK
using JuMP
using Random
using Test

"""
    example_lazy_constraint()

An example using a lazy constraint callback.
"""
function example_lazy_constraint()
    model = Model(GLPK.Optimizer)
    @variable(model, 0 <= x <= 2.5, Int)
    @variable(model, 0 <= y <= 2.5, Int)
    @objective(model, Max, y)
    lazy_called = false
    function my_callback_function(cb_data)
        lazy_called = true
        x_val = callback_value(cb_data, x)
        y_val = callback_value(cb_data, y)
        if y_val - x_val > 1 + 1e-6
            con = @build_constraint(y - x <= 1)
            MOI.submit(model, MOI.LazyConstraint(cb_data), con)
        elseif y_val + x_val > 3 + 1e-6
            con = @build_constraint(y - x <= 1)
            MOI.submit(model, MOI.LazyConstraint(cb_data), con)
        end
    end
    MOI.set(model, MOI.LazyConstraintCallback(), my_callback_function)
    optimize!(model)
    @test termination_status(model) == MOI.OPTIMAL
    @test primal_status(model) == MOI.FEASIBLE_POINT
    @test lazy_called
    @test value(x) == 1
    @test value(y) == 2
end

example_lazy_constraint()

"""
    example_user_cut_constraint()

An example using a user-cut callback.
"""
function example_user_cut_constraint()
    Random.seed!(1)
    N = 30
    item_weights, item_values = rand(N), rand(N)

    model = Model(GLPK.Optimizer)
    @variable(model, x[1:N], Bin)
    @constraint(model, sum(item_weights[i] * x[i] for i = 1:N) <= 10)
    @objective(model, Max, sum(item_values[i] * x[i] for i = 1:N))

    callback_called = false
    function my_callback_function(cb_data)
        callback_called = true
        # TODO(odow): remove Ref once GLPK supports broadcasting over cb_data.
        x_vals = callback_value.(Ref(cb_data), x)
        accumulated = sum(item_weights[i] for i=1:N if x_vals[i] > 1e-4)
        n_terms = sum(1 for i=1:N if x_vals[i] > 1e-4)
        if accumulated > 10
            con = @build_constraint(
                sum(x[i] for i = 1:N if x_vals[i] > 0.5) <= n_terms - 1
            )
            MOI.submit(model, MOI.UserCut(cb_data), con)
        end
    end
    MOI.set(model, MOI.UserCutCallback(), my_callback_function)
    optimize!(model)
    @test termination_status(model) == MOI.OPTIMAL
    @test primal_status(model) == MOI.FEASIBLE_POINT
    @test callback_called
end

example_user_cut_constraint()

"""
    example_heuristic_solution()

An example using a heuristic solution callback.
"""
function example_heuristic_solution()
    Random.seed!(1)
    N = 30
    item_weights, item_values = rand(N), rand(N)

    model = Model(GLPK.Optimizer)
    @variable(model, x[1:N], Bin)
    @constraint(model, sum(item_weights[i] * x[i] for i = 1:N) <= 10)
    @objective(model, Max, sum(item_values[i] * x[i] for i = 1:N))

    callback_called = false
    function my_callback_function(cb_data)
        callback_called = true
        # TODO(odow): remove Ref once GLPK supports broadcasting over cb_data.
        x_vals = callback_value.(Ref(cb_data), x)
        @test MOI.submit(
            model, MOI.HeuristicSolution(cb_data), x, floor.(x_vals)
        ) in (MOI.HEURISTIC_SOLUTION_ACCEPTED, MOI.HEURISTIC_SOLUTION_REJECTED)
    end
    MOI.set(model, MOI.HeuristicCallback(), my_callback_function)
    optimize!(model)
    @test termination_status(model) == MOI.OPTIMAL
    @test primal_status(model) == MOI.FEASIBLE_POINT
    @test callback_called
end

example_heuristic_solution()

"""
    example_solver_dependent_callback()

An example using a solver_dependent callback.
"""
function example_solver_dependent_callback()
    model = Model(GLPK.Optimizer)
    @variable(model, 0 <= x <= 2.5, Int)
    @variable(model, 0 <= y <= 2.5, Int)
    @objective(model, Max, y)
    lazy_called = false
    function my_callback_function(cb_data)
        lazy_called = true
        reason = GLPK.ios_reason(cb_data.tree)
        if reason != GLPK.IROWGEN
            return
        end
        x_val = callback_value(cb_data, x)
        y_val = callback_value(cb_data, y)
        if y_val - x_val > 1 + 1e-6
            con = @build_constraint(y - x <= 1)
            MOI.submit(model, MOI.LazyConstraint(cb_data), con)
        elseif y_val + x_val > 3 + 1e-6
            con = @build_constraint(y - x <= 1)
            MOI.submit(model, MOI.LazyConstraint(cb_data), con)
        end
    end
    MOI.set(model, GLPK.CallbackFunction(), my_callback_function)
    optimize!(model)
    @test termination_status(model) == MOI.OPTIMAL
    @test primal_status(model) == MOI.FEASIBLE_POINT
    @test lazy_called
    @test value(x) == 1
    @test value(y) == 2
end

example_solver_dependent_callback()

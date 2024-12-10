# Copyright 2017, Iain Dunning, Joey Huchette, Miles Lubin, and contributors    #src
# This Source Code Form is subject to the terms of the Mozilla Public License   #src
# v.2.0. If a copy of the MPL was not distributed with this file, You can       #src
# obtain one at https://mozilla.org/MPL/2.0/.                                   #src

# # [Callbacks](@id callbacks_tutorial)

# The purpose of the tutorial is to demonstrate the various solver-independent
# and solver-dependent callbacks that are supported by JuMP.

# The tutorial uses the following packages:

using JuMP
import GLPK
import Random
import Test  #src

# !!! info
#     This tutorial uses the [MathOptInterface](@ref moi_documentation) API.
#     By default, JuMP exports the `MOI` symbol as an alias for the
#     MathOptInterface.jl package. We recommend making this more explicit in
#     your code by adding the following lines:
#     ```julia
#     import MathOptInterface as MOI
#     ```

# ## Lazy constraints

# An example using a lazy constraint callback.

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
        println("Called from (x, y) = ($x_val, $y_val)")
        status = callback_node_status(cb_data, model)
        if status == MOI.CALLBACK_NODE_STATUS_FRACTIONAL
            println(" - Solution is integer infeasible!")
        elseif status == MOI.CALLBACK_NODE_STATUS_INTEGER
            println(" - Solution is integer feasible!")
        else
            @assert status == MOI.CALLBACK_NODE_STATUS_UNKNOWN
            println(" - I don't know if the solution is integer feasible :(")
        end
        if y_val - x_val > 1 + 1e-6
            con = @build_constraint(y - x <= 1)
            println("Adding $(con)")
            MOI.submit(model, MOI.LazyConstraint(cb_data), con)
        elseif y_val + x_val > 3 + 1e-6
            con = @build_constraint(y + x <= 3)
            println("Adding $(con)")
            MOI.submit(model, MOI.LazyConstraint(cb_data), con)
        end
    end
    MOI.set(model, MOI.LazyConstraintCallback(), my_callback_function)
    optimize!(model)
    Test.@test termination_status(model) == OPTIMAL    #src
    Test.@test primal_status(model) == FEASIBLE_POINT  #src
    Test.@test lazy_called    #src
    Test.@test value(x) == 1  #src
    Test.@test value(y) == 2  #src
    println("Optimal solution (x, y) = ($(value(x)), $(value(y)))")
    return
end

example_lazy_constraint()

# ## User-cuts

# An example using a user-cut callback.

function example_user_cut_constraint()
    Random.seed!(1)
    N = 30
    item_weights, item_values = rand(N), rand(N)
    model = Model(GLPK.Optimizer)
    @variable(model, x[1:N], Bin)
    @constraint(model, sum(item_weights[i] * x[i] for i in 1:N) <= 10)
    @objective(model, Max, sum(item_values[i] * x[i] for i in 1:N))
    callback_called = false
    function my_callback_function(cb_data)
        callback_called = true
        x_vals = callback_value.(Ref(cb_data), x)
        accumulated = sum(item_weights[i] for i in 1:N if x_vals[i] > 1e-4)
        println("Called with accumulated = $(accumulated)")
        n_terms = sum(1 for i in 1:N if x_vals[i] > 1e-4)
        if accumulated > 10
            con = @build_constraint(
                sum(x[i] for i in 1:N if x_vals[i] > 0.5) <= n_terms - 1
            )
            println("Adding $(con)")
            MOI.submit(model, MOI.UserCut(cb_data), con)
        end
    end
    MOI.set(model, MOI.UserCutCallback(), my_callback_function)
    optimize!(model)
    Test.@test termination_status(model) == OPTIMAL  #src
    Test.@test primal_status(model) == FEASIBLE_POINT  #src
    Test.@test callback_called  #src
    @show callback_called
    return
end

example_user_cut_constraint()

# ## Heuristic solutions

# An example using a heuristic solution callback.

function example_heuristic_solution()
    Random.seed!(1)
    N = 30
    item_weights, item_values = rand(N), rand(N)
    model = Model(GLPK.Optimizer)
    @variable(model, x[1:N], Bin)
    @constraint(model, sum(item_weights[i] * x[i] for i in 1:N) <= 10)
    @objective(model, Max, sum(item_values[i] * x[i] for i in 1:N))
    callback_called = false
    function my_callback_function(cb_data)
        callback_called = true
        x_vals = callback_value.(Ref(cb_data), x)
        ret =
            MOI.submit(model, MOI.HeuristicSolution(cb_data), x, floor.(x_vals))
        println("Heuristic solution status = $(ret)")
        Test.@test ret in (                     #src
            MOI.HEURISTIC_SOLUTION_ACCEPTED,    #src
            MOI.HEURISTIC_SOLUTION_REJECTED,    #src
        )                                       #src
    end
    MOI.set(model, MOI.HeuristicCallback(), my_callback_function)
    optimize!(model)
    Test.@test termination_status(model) == OPTIMAL  #src
    Test.@test primal_status(model) == FEASIBLE_POINT  #src
    Test.@test callback_called  #src
    return
end

example_heuristic_solution()

# ## GLPK solver-dependent callback

# An example using GLPK's solver-dependent callback.

function example_solver_dependent_callback()
    model = Model(GLPK.Optimizer)
    @variable(model, 0 <= x <= 2.5, Int)
    @variable(model, 0 <= y <= 2.5, Int)
    @objective(model, Max, y)
    lazy_called = false
    function my_callback_function(cb_data)
        lazy_called = true
        reason = GLPK.glp_ios_reason(cb_data.tree)
        println("Called from reason = $(reason)")
        if reason != GLPK.GLP_IROWGEN
            return
        end
        x_val = callback_value(cb_data, x)
        y_val = callback_value(cb_data, y)
        if y_val - x_val > 1 + 1e-6
            con = @build_constraint(y - x <= 1)
            println("Adding $(con)")
            MOI.submit(model, MOI.LazyConstraint(cb_data), con)
        elseif y_val + x_val > 3 + 1e-6
            con = @build_constraint(y - x <= 1)
            println("Adding $(con)")
            MOI.submit(model, MOI.LazyConstraint(cb_data), con)
        end
    end
    MOI.set(model, GLPK.CallbackFunction(), my_callback_function)
    optimize!(model)
    Test.@test termination_status(model) == OPTIMAL  #src
    Test.@test primal_status(model) == FEASIBLE_POINT  #src
    Test.@test lazy_called  #src
    Test.@test value(x) == 1  #src
    Test.@test value(y) == 2  #src
    return
end

example_solver_dependent_callback()

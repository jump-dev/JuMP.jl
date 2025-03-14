# Copyright 2017, Iain Dunning, Joey Huchette, Miles Lubin, and contributors    #src
# This Source Code Form is subject to the terms of the Mozilla Public License   #src
# v.2.0. If a copy of the MPL was not distributed with this file, You can       #src
# obtain one at https://mozilla.org/MPL/2.0/.                                   #src

# # [Callbacks](@id callbacks_tutorial)

# The purpose of the tutorial is to demonstrate the various solver-independent
# and solver-dependent callbacks that are supported by JuMP.

# The tutorial uses the following packages:

using JuMP
import Gurobi
import Random
import Test

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
    model = Model(Gurobi.Optimizer)
    set_silent(model)
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
        return
    end
    set_attribute(model, MOI.LazyConstraintCallback(), my_callback_function)
    optimize!(model)
    assert_is_solved_and_feasible(model)
    Test.@test lazy_called
    Test.@test value(x) == 1
    Test.@test value(y) == 2
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
    model = Model(Gurobi.Optimizer)
    set_silent(model)
    ## Turn off "Cuts" parameter so that our new one must be called. In real
    ## models, you should leave "Cuts" turned on.
    set_attribute(model, "Cuts", 0)
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
    set_attribute(model, MOI.UserCutCallback(), my_callback_function)
    optimize!(model)
    assert_is_solved_and_feasible(model)
    Test.@test callback_called
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
    model = Model(Gurobi.Optimizer)
    set_silent(model)
    ## Turn off "Heuristics" parameter so that our new one must be called. In
    ## real models, you should leave "Heuristics" turned on.
    set_attribute(model, "Heuristics", 0)
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
        Test.@test ret in (
            MOI.HEURISTIC_SOLUTION_ACCEPTED,
            MOI.HEURISTIC_SOLUTION_REJECTED,
        )
    end
    set_attribute(model, MOI.HeuristicCallback(), my_callback_function)
    optimize!(model)
    assert_is_solved_and_feasible(model)
    Test.@test callback_called
    return
end

example_heuristic_solution()

# ## Gurobi solver-dependent callback

# An example using Gurobi's solver-dependent callback.

function example_solver_dependent_callback()
    model = direct_model(Gurobi.Optimizer())
    @variable(model, 0 <= x <= 2.5, Int)
    @variable(model, 0 <= y <= 2.5, Int)
    @objective(model, Max, y)
    cb_calls = Cint[]
    function my_callback_function(cb_data, cb_where::Cint)
        ## You can reference variables outside the function as normal
        push!(cb_calls, cb_where)
        ## You can select where the callback is run
        if cb_where == Gurobi.GRB_CB_MIPNODE
            ## You can query a callback attribute using GRBcbget
            resultP = Ref{Cint}()
            Gurobi.GRBcbget(
                cb_data,
                cb_where,
                Gurobi.GRB_CB_MIPNODE_STATUS,
                resultP,
            )
            if resultP[] != Gurobi.GRB_OPTIMAL
                return  # Solution is something other than optimal.
            end
        elseif cb_where != Gurobi.GRB_CB_MIPSOL
            return
        end
        ## Before querying `callback_value`, you must call:
        Gurobi.load_callback_variable_primal(cb_data, cb_where)
        x_val = callback_value(cb_data, x)
        y_val = callback_value(cb_data, y)
        ## You can submit solver-independent MathOptInterface attributes such as
        ## lazy constraints, user-cuts, and heuristic solutions.
        if y_val - x_val > 1 + 1e-6
            con = @build_constraint(y - x <= 1)
            MOI.submit(model, MOI.LazyConstraint(cb_data), con)
        elseif y_val + x_val > 3 + 1e-6
            con = @build_constraint(y + x <= 3)
            MOI.submit(model, MOI.LazyConstraint(cb_data), con)
        end
        ## You can terminate the callback as follows:
        Gurobi.GRBterminate(backend(model))
        return
    end
    ## You _must_ set this parameter if using lazy constraints.
    set_attribute(model, "LazyConstraints", 1)
    set_attribute(model, Gurobi.CallbackFunction(), my_callback_function)
    optimize!(model)
    Test.@test termination_status(model) == MOI.INTERRUPTED
    return
end

example_solver_dependent_callback()

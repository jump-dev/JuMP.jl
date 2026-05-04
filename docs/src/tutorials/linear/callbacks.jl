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
import Ipopt
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

# Here is an example of using a lazy constraint callback with Gurobi. For more
# information about lazy constraints, see the [Lazy constraints](@ref manual_lazy_constraints)
# section of the manual.

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

# Here is an example of using a user cut callback with Gurobi. For more
# information about user cuts, see the [User cuts](@ref manual_user_cuts)
# section of the manual.

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

# Here is an example of using a heuristic solution callback with Gurobi. For
# more information about heuristic solutions, see the
# [Heuristic solutions](@ref manual_heuristic_solutions) section of the manual.

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

# ## Solver-dependent callback

# Some solvers expose solver-dependent callbacks. The syntax of the callback
# depends on the solver. Typically, this requires you to interact with the
# low-level C API of the solver. If a solver supports a solver-dependent
# callback this will be documented in the README of the solver wrapper.

# Here's an example of Gurobi's:

function example_solver_dependent_callback_gurobi()
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

example_solver_dependent_callback_gurobi()

# And here is an example of Ipopt's:

function example_solver_dependent_callback_ipopt()
    model = Model(Ipopt.Optimizer)
    set_silent(model)
    @variable(model, x >= 1)
    @objective(model, Min, x + 0.5)
    x_vals = Float64[]
    function my_callback(
        alg_mod::Cint,
        iter_count::Cint,
        obj_value::Float64,
        inf_pr::Float64,
        inf_du::Float64,
        mu::Float64,
        d_norm::Float64,
        regularization_size::Float64,
        alpha_du::Float64,
        alpha_pr::Float64,
        ls_trials::Cint,
    )
        push!(x_vals, callback_value(model, x))
        Test.@test isapprox(obj_value, 1.0 * x_vals[end] + 0.5, atol = 1e-1)
        ## return `true` to keep going, or `false` to terminate the optimization.
        return iter_count < 1
    end
    set_attribute(model, Ipopt.CallbackFunction(), my_callback)
    optimize!(model)
    termination_status(model)
    Test.@test length(x_vals) == 2
    return x_vals
end

example_solver_dependent_callback_ipopt()

# ### Using solver-dependent callbacks

# If you want to write a package that is solver-independent, but you also want
# to use a solver-dependent callback, write a [package extension](https://pkgdocs.julialang.org/v1/creating-packages/#Conditional-loading-of-code-in-packages-(Extensions)).

# To proceed, assume that we have a package named `MyPackage`:

module MyPackage
import JuMP
add_callback(model) = add_callback(model, typeof(JuMP.unsafe_backend(model)))
add_callback(model, ::Type{T}) where {T} = error("Unsupported optimizer: $T")
function build_model(optimizer)
    model = JuMP.Model(optimizer)
    add_callback(model)
    return model
end
end  # MyPackage

# This package defines a public `build_model` function. Inside the function it
# calls `add_callback`, which defaults to throwing an error.

# Now, assume we want to add callbacks for Gurobi and Ipopt.

# In `/ext/MyPackageGurobiExt.jl` define:

module MyPackageGurobiExt
## In a real example, change `..MyPackage` to `MyPackage`. The dots are needed
## only because of how we've structured this documentation.
import ..MyPackage
import Gurobi
import JuMP
function MyPackage.add_callback(model::JuMP.Model, ::Type{Gurobi.Optimizer})
    function callback(cb_data, cb_where)
        if rand() < 0.5
            Gurobi.GRBterminate(JuMP.backend(model))
        end
        return
    end
    JuMP.set_attribute(model, Gurobi.CallbackFunction(), callback)
    return
end
end  # MyPackageGurobiExt

# and in `/ext/MyPackageIpoptExt.jl` define:

module MyPackageIpoptExt
## In a real example, change `..MyPackage` to `MyPackage`. The dots are needed
## only because of how we've structured this documentation.
import ..MyPackage
import Ipopt
import JuMP
function MyPackage.add_callback(model::JuMP.Model, ::Type{Ipopt.Optimizer})
    function callback(args...)
        return rand() < 0.5
    end
    JuMP.set_attribute(model, Ipopt.CallbackFunction(), callback)
    return
end
end  # MyPackageIpoptExt

# Finally, change the `Project.toml` of `MyPackage` to include an `[extensions]`
# section:
#
# ```toml
# [weakdeps]
# Gurobi = "2e9cd046-0924-5485-92f1-d5272153d98b"
# Ipopt = "b6b21f68-93f8-5de0-b562-5493be1d77c9"
#
# [extensions]
# MyPackageGurobiExt = "Gurobi"
# MyPackageIpoptExt = "Ipopt"
# ```

# Now, `build_model` can be called with `Gurobi.Optimizer`:

MyPackage.build_model(Gurobi.Optimizer)

# and with `Ipopt.Optimizer`:

MyPackage.build_model(Ipopt.Optimizer)

# And you didn't need to add Gurobi or Ipopt as dependencies to MyPackage.

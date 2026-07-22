# Copyright (c) 2019 Arpit Bhatia and contributors                               #src
#                                                                                #src
# Permission is hereby granted, free of charge, to any person obtaining a copy   #src
# of this software and associated documentation files (the "Software"), to deal  #src
# in the Software without restriction, including without limitation the rights   #src
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell      #src
# copies of the Software, and to permit persons to whom the Software is          #src
# furnished to do so, subject to the following conditions:                       #src
#                                                                                #src
# The above copyright notice and this permission notice shall be included in all #src
# copies or substantial portions of the Software.                                #src
#                                                                                #src
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR     #src
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,       #src
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE    #src
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER         #src
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,  #src
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE  #src
# SOFTWARE.                                                                      #src

# # [Benders decomposition](@id benders_decomposition_classical)

# **This tutorial was originally contributed by Shuvomoy Das Gupta.**

# This tutorial describes how to implement [Benders decomposition](https://en.wikipedia.org/wiki/Benders_decomposition)
# in JuMP. It uses the following packages:

using JuMP
import GLPK
import HiGHS
import Printf
import Test  #src

# ## Theory

# Benders decomposition is a useful algorithm for solving convex optimization
# problems with a large number of variables. It works best when a larger problem
# can be decomposed into two (or more) smaller problems that are individually
# much easier to solve.
#
# This tutorial demonstrates Benders decomposition on the following
# mixed-integer linear program:
# ```math
# \begin{aligned}
# \text{min}        \ & c_1(x) + c_2(y) \\
# \text{subject to} \ & f_1(x) \in S_1 \\
#                     & f_2(y) \in S_2 \\
#                     & f_3(x, y) \in S_3 \\
#                     & x \in \mathbb{Z}^m \\
#                     & y \in \mathbb{R}^n \\
# \end{aligned}
# ```
# where the functions $f$ and $c$ are linear, and the sets $S$ are inequality
# sets like $\ge l$, $\le u$, or $= b$.
#
# Any mixed integer programming problem can be written in the form above.

# If there are relatively few integer variables, and many more continuous
# variables, then it may be beneficial to decompose the problem into a small
# problem containing only integer variables and a linear program containing
# only continuous variables. Hopefully, the linear program will be much easier
# to solve in isolation than in the full mixed-integer linear program.

# For example, if we knew a feasible solution for ``\bar{x}``, we could obtain a
# solution for ``y`` by solving:
# ```math
# \begin{aligned}
# V_2(\bar{x}) = & \text{min}  \ & c_2(y)\\
#          & \text{subject to} \ & f_2(y) \in S_2 \\
#          &                     & f_3(x, y) \in S_3 \\
#          &                     & x = \bar{x} & \ [\pi] \\
#          &                     & y \in \mathbb{R}^n \\
# \end{aligned}
# ```
# Note that we have included a "copy" of the `x` variable to simplify computing
# $\pi$, which is the dual of $V_2$ with respect to $\bar{x}$.
#
# Because this model is a linear program, it is easy to solve.

# Replacing the ``c_2(y)`` component of the objective in our original problem
# with ``V_2`` yields:
# ```math
# \begin{aligned}
# V_1 = & \text{min}        \ & c_1(x) + V_2(x) \\
# & \text{subject to} \ & f_1(x) \in S_1 \\
# &                     & x \in \mathbb{Z}^m.
# \end{aligned}
# ```
# This problem looks a lot simpler to solve because it involves only $x$ and a
# subset of the constraints, but we need to do something else with ``V_2`` first.

# Because ``\bar{x}`` is a constant that appears on the right-hand side of the
# constraints, ``V_2`` is a convex function with respect to ``\bar{x}``, and the
# dual variable ``\pi`` is a subgradient of ``V_2(x)`` with respect to ``x``.
# Therefore, if we have a candidate solution ``x_k``, then we can solve
# ``V_2(x_k)`` and obtain a feasible dual vector ``\pi_k``. Using these values,
# we can construct a first-order Taylor-series approximation of ``V_2`` about
# the point ``x_k``:
# ```math
# V_2(x) \ge V_2(x_k) + \pi_k^\top (x - x_k).
# ```
# By convexity, we know that this inequality holds for all ``x``, and we call
# these inequalities _cuts_.

# Benders decomposition is an iterative technique that replaces ``V_2(x)`` with
# a new decision variable ``\theta``, and approximates it from below using cuts:
# ```math
# \begin{aligned}
# V_1^K = & \text{min}        \ & c_1(x) + \theta      \\
#         & \text{subject to} \ & f_1(x) \in S_1       \\
#         &                     & x \in \mathbb{Z}^m   \\
#         &                     & \theta \ge M         \\
#         &                     & \theta \ge V_2(x_k) + \pi_k^\top(x - x_k) & \quad \forall k = 1,\ldots,K.
# \end{aligned}
# ```
# This integer program is called the _first-stage_ subproblem.

# To generate cuts, we solve ``V_1^K`` to obtain a candidate first-stage
# solution ``x_k``, then we use that solution to solve ``V_2(x_k)``. Then, using
# the optimal objective value and dual solution from ``V_2``, we add a new cut
# to form ``V_1^{K+1}`` and repeat.

# ### Bounds

# Due to convexity, we know that ``V_2(x) \ge \theta`` for all ``x``. Therefore,
# the optimal objective value of ``V_1^K`` provides a valid _lower_ bound on the
# objective value of the full problem. In addition, if we take a feasible
# solution for ``x`` from the first-stage problem, then ``c_1(x) + V_2(x)``
# is a valid _upper_ bound on the objective value of the full problem.

# Benders decomposition uses the lower and upper bounds to determine when it has
# found the global optimal solution.

# ## Monolithic problem

# As an example problem, we consider the following variant of
# [The max-flow problem](@ref), in which there is a binary variable to decide
# whether to open each arc for a cost of 0.1 unit, and we can open at most 11
# arcs:

G = [
    0 3 2 2 0 0 0 0
    0 0 0 0 5 1 0 0
    0 0 0 0 1 3 1 0
    0 0 0 0 0 1 0 0
    0 0 0 0 0 0 0 4
    0 0 0 0 0 0 0 2
    0 0 0 0 0 0 0 4
    0 0 0 0 0 0 0 0
]
n = size(G, 1)
model = Model(HiGHS.Optimizer)
set_silent(model)
@variable(model, x[1:n, 1:n], Bin)
@variable(model, y[1:n, 1:n] >= 0)
@constraint(model, sum(x) <= 11)
@constraint(model, [i = 1:n, j = 1:n], y[i, j] <= G[i, j] * x[i, j])
@constraint(model, [i = 2:n-1], sum(y[i, :]) == sum(y[:, i]))
@objective(model, Min, 0.1 * sum(x) - sum(y[1, :]))
optimize!(model)
Test.@test is_solved_and_feasible(model)  #src
solution_summary(model)

# The optimal objective value is -5.1:

Test.@test isapprox(objective_value(model), -5.1; atol = 1e-4)  #src
objective_value(model)

# and the optimal flows are:

function optimal_flows(x)
    return [(i, j) => x[i, j] for i in 1:n for j in 1:n if x[i, j] > 0]
end

monolithic_solution = optimal_flows(value.(y))

# ## [Iterative method](@id benders_iterative)

# !!! warning
#     This is a basic implementation for pedagogical purposes. We haven't
#     discussed any of the computational tricks that are required to build a
#     performant implementation for large-scale problems. See [In-place iterative method](@ref)
#     for one improvement that helps computation time.

# We start by formulating the first-stage subproblem. It includes the `x`
# variables, and the constraints involving only `x`, and the terms in the
# objective containing only `x`. We also need an initial lower bound on the
# cost-to-go variable `θ`. One valid lower bound is to assume that we do not pay
# for opening arcs, and there is flow all the arcs.

M = -sum(G)
model = Model(HiGHS.Optimizer)
set_silent(model)
@variable(model, x[1:n, 1:n], Bin)
@variable(model, θ >= M)
@constraint(model, sum(x) <= 11)
@objective(model, Min, 0.1 * sum(x) + θ)
model

# For the next step, we need a function that takes a first-stage candidate
# solution `x` and returns the optimal solution from the second-stage
# subproblem:

function solve_subproblem(x_bar)
    model = Model(HiGHS.Optimizer)
    set_silent(model)
    @variable(model, x[i in 1:n, j in 1:n] == x_bar[i, j])
    @variable(model, y[1:n, 1:n] >= 0)
    @constraint(model, [i = 1:n, j = 1:n], y[i, j] <= G[i, j] * x[i, j])
    @constraint(model, [i = 2:n-1], sum(y[i, :]) == sum(y[:, i]))
    @objective(model, Min, -sum(y[1, :]))
    optimize!(model)
    @assert is_solved_and_feasible(model; dual = true)
    return (obj = objective_value(model), y = value.(y), π = reduced_cost.(x))
end

# Note that `solve_subproblem` returns a `NamedTuple` of the objective value,
# the optimal primal solution for `y`, and the optimal dual solution for `π`,
# which we obtained from the [`reduced_cost`](@ref) of the `x` variables.

# We're almost ready for our optimization loop, but first, here's a helpful
# function for logging:

function print_iteration(k, args...)
    f(x) = Printf.@sprintf("%12.4e", x)
    println(lpad(k, 9), " ", join(f.(args), " "))
    return
end

# We also need to put a limit on the number of iterations before termination:

MAXIMUM_ITERATIONS = 100

# And a way to check if the lower and upper bounds are close-enough to
# terminate:

ABSOLUTE_OPTIMALITY_GAP = 1e-6

# Now we're ready to iterate Benders decomposition:

println("Iteration  Lower Bound  Upper Bound          Gap")
for k in 1:MAXIMUM_ITERATIONS
    optimize!(model)
    @assert is_solved_and_feasible(model)
    lower_bound = objective_value(model)
    x_k = value.(x)
    ret = solve_subproblem(x_k)
    upper_bound = (objective_value(model) - value(θ)) + ret.obj
    gap = abs(upper_bound - lower_bound) / abs(upper_bound)
    print_iteration(k, lower_bound, upper_bound, gap)
    if gap < ABSOLUTE_OPTIMALITY_GAP
        println("Terminating with the optimal solution")
        break
    end
    cut = @constraint(model, θ >= ret.obj + sum(ret.π .* (x .- x_k)))
    @info "Adding the cut $(cut)"
end

# Finally, we can obtain the optimal solution:

optimize!(model)
@assert is_solved_and_feasible(model)
x_optimal = value.(x)
optimal_ret = solve_subproblem(x_optimal)
iterative_solution = optimal_flows(optimal_ret.y)

# which is the same as the monolithic solution:

iterative_solution == monolithic_solution

# and it has the same objective value:

objective_value(model)

# ## Callback method

# The [Iterative method](@ref benders_iterative) section implemented Benders
# decomposition using a loop. In each iteration, we re-solved the first-stage
# subproblem to generate a candidate solution. However, modern MILP solvers such
# as CPLEX, Gurobi, and GLPK provide lazy constraint callbacks which allow us to
# add new cuts _while the solver is running_. This can be more efficient than an
# iterative method because we can avoid repeating work such as solving the root
# node of the first-stage MILP at each iteration.

# !!! tip
#     We use GLPK for this model because HiGHS does not support lazy constraints.
#     For more information on callbacks, read the page
#     [Solver-independent callbacks](@ref callbacks_manual).

# As before, we construct the same first-stage subproblem:

lazy_model = Model(GLPK.Optimizer)
set_silent(lazy_model)
@variable(lazy_model, x[1:n, 1:n], Bin)
@variable(lazy_model, θ >= M)
@constraint(lazy_model, sum(x) <= 11)
@objective(lazy_model, Min, 0.1 * sum(x) + θ)
lazy_model

# What differs is that we write a callback function instead of a loop:

number_of_subproblem_solves = 0
function my_callback(cb_data)
    status = callback_node_status(cb_data, lazy_model)
    if status != MOI.CALLBACK_NODE_STATUS_INTEGER
        ## Only add the constraint if `x` is an integer feasible solution
        return
    end
    x_k = callback_value.(cb_data, x)
    θ_k = callback_value(cb_data, θ)
    global number_of_subproblem_solves += 1
    ret = solve_subproblem(x_k)
    if θ_k < (ret.obj - 1e-6)
        ## Only add the constraint if θ_k violates the constraint
        cut = @build_constraint(θ >= ret.obj + sum(ret.π .* (x .- x_k)))
        MOI.submit(lazy_model, MOI.LazyConstraint(cb_data), cut)
    end
    return
end

set_attribute(lazy_model, MOI.LazyConstraintCallback(), my_callback)

# Now when we optimize!, our callback is run:

optimize!(lazy_model)
@assert is_solved_and_feasible(lazy_model)

# For this model, the callback algorithm required more solves of the subproblem:

number_of_subproblem_solves

# But for larger problems, you can expect the callback algorithm to be more
# efficient than the iterative algorithm.

# Finally, we can obtain the optimal solution:

x_optimal = value.(x)
optimal_ret = solve_subproblem(x_optimal)
callback_solution = optimal_flows(optimal_ret.y)

# which is the same as the monolithic solution:

Test.@test callback_solution == monolithic_solution  #src
callback_solution == monolithic_solution

# ## In-place iterative method

# Our implementation of the iterative method has a problem: every time we need
# to solve the subproblem, we must rebuild it from scratch. This is expensive,
# and it can be the bottleneck in the solution process. We can improve our
# implementation by using re-using the subproblem between solves.

# First, we create our first-stage problem as usual:

model = Model(HiGHS.Optimizer)
set_silent(model)
@variable(model, x[1:n, 1:n], Bin)
@variable(model, θ >= M)
@constraint(model, sum(x) <= 11)
@objective(model, Min, 0.1 * sum(x) + θ)
model

# Then, instead of building the subproblem in a function, we build it once here:

subproblem = Model(HiGHS.Optimizer)
set_silent(subproblem)
@variable(subproblem, x_copy[i in 1:n, j in 1:n])
@variable(subproblem, y[1:n, 1:n] >= 0)
@constraint(subproblem, [i = 1:n, j = 1:n], y[i, j] <= G[i, j] * x_copy[i, j])
@constraint(subproblem, [i = 2:n-1], sum(y[i, :]) == sum(y[:, i]))
@objective(subproblem, Min, -sum(y[1, :]))
subproblem

# Our function to solve the subproblem is also slightly different because we
# need to fix the value of the `x_copy` variables to the value of `x` from the
# first-stage problem:

function solve_subproblem(model, x)
    fix.(model[:x_copy], x)
    optimize!(model)
    @assert is_solved_and_feasible(model; dual = true)
    return (
        obj = objective_value(model),
        y = value.(model[:y]),
        π = reduced_cost.(model[:x_copy]),
    )
end

# Now we're ready to iterate our in-place Benders decomposition:

println("Iteration  Lower Bound  Upper Bound          Gap")
for k in 1:MAXIMUM_ITERATIONS
    optimize!(model)
    @assert is_solved_and_feasible(model)
    lower_bound = objective_value(model)
    x_k = value.(x)
    ret = solve_subproblem(subproblem, x_k)
    upper_bound = (objective_value(model) - value(θ)) + ret.obj
    gap = abs(upper_bound - lower_bound) / abs(upper_bound)
    print_iteration(k, lower_bound, upper_bound, gap)
    if gap < ABSOLUTE_OPTIMALITY_GAP
        println("Terminating with the optimal solution")
        break
    end
    cut = @constraint(model, θ >= ret.obj + sum(ret.π .* (x .- x_k)))
    @info "Adding the cut $(cut)"
end

# Finally, we can obtain the optimal solution:

optimize!(model)
@assert is_solved_and_feasible(model)
x_optimal = value.(x)
optimal_ret = solve_subproblem(subproblem, x_optimal)
inplace_solution = optimal_flows(optimal_ret.y)

# which is the same as the monolithic solution:

Test.@test inplace_solution == monolithic_solution  #src
inplace_solution == monolithic_solution

# ## Feasibility cuts

# So far, we have discussed only Benders optimality cuts. However, for some
# first-stage values of `x`, the subproblem might be infeasible. The solution is
# to add a Benders feasibility cut:
# ```math
# v_k + u_k^\top (x - x_k) \le 0
# ```
# where $u_k$ is a dual unbounded ray of the subproblem and $v_k$ is the
# intercept of the unbounded ray.

# As a variation of our example which leads to infeasibilities, we add a
# constraint that `sum(y) >= 1`. This means we need a choice of first-stage `x`
# for which at least one unit can flow.

# The first-stage problem remains the same:

model = Model(HiGHS.Optimizer)
set_silent(model)
@variable(model, x[1:n, 1:n], Bin)
@variable(model, θ >= M)
@constraint(model, sum(x) <= 11)
@objective(model, Min, 0.1 * sum(x) + θ)
model

# But the subproblem has a new constraint that `sum(y) >= 1`:

subproblem = Model(HiGHS.Optimizer)
set_silent(subproblem)
## We need to turn presolve off so that HiGHS will return an infeasibility
## certificate.
set_attribute(subproblem, "presolve", "off")
@variable(subproblem, x_copy[i in 1:n, j in 1:n])
@variable(subproblem, y[1:n, 1:n] >= 0)
@constraint(subproblem, sum(y) >= 1)  # <--- THIS IS NEW
@constraint(subproblem, [i = 1:n, j = 1:n], y[i, j] <= G[i, j] * x_copy[i, j])
@constraint(subproblem, [i = 2:n-1], sum(y[i, :]) == sum(y[:, i]))
@objective(subproblem, Min, -sum(y[1, :]))
subproblem

# The function to solve the subproblem now checks for feasibility, and returns
# the dual objective value and an dual unbounded ray if the subproblem is
# infeasible:

function solve_subproblem_with_feasibility(model, x)
    fix.(model[:x_copy], x)
    optimize!(model)
    if is_solved_and_feasible(model; dual = true)
        return (
            is_feasible = true,
            obj = objective_value(model),
            y = value.(model[:y]),
            π = reduced_cost.(model[:x_copy]),
        )
    end
    return (
        is_feasible = false,
        v = dual_objective_value(model),
        u = reduced_cost.(model[:x_copy]),
    )
end

# Now we're ready to iterate our in-place Benders decomposition:

println("Iteration  Lower Bound  Upper Bound          Gap")
for k in 1:MAXIMUM_ITERATIONS
    optimize!(model)
    @assert is_solved_and_feasible(model)
    lower_bound = objective_value(model)
    x_k = value.(x)
    ret = solve_subproblem_with_feasibility(subproblem, x_k)
    if ret.is_feasible
        ## Benders Optimality Cuts
        upper_bound = (objective_value(model) - value(θ)) + ret.obj
        gap = abs(upper_bound - lower_bound) / abs(upper_bound)
        print_iteration(k, lower_bound, upper_bound, gap)
        if gap < ABSOLUTE_OPTIMALITY_GAP
            println("Terminating with the optimal solution")
            break
        end
        @constraint(model, θ >= ret.obj + sum(ret.π .* (x .- x_k)))
    else
        ## Benders Feasibility Cuts
        cut = @constraint(model, ret.v + sum(ret.u .* (x .- x_k)) <= 0)
        @info "Adding the feasibility cut $(cut)"
    end
end

# Finally, we can obtain the optimal solution:

optimize!(model)
@assert is_solved_and_feasible(model)
x_optimal = value.(x)
optimal_ret = solve_subproblem(subproblem, x_optimal)
feasible_inplace_solution = optimal_flows(optimal_ret.y)

# which is the same as the monolithic solution (because `sum(y) >= 1` in the
# monolithic solution):

Test.@test feasible_inplace_solution == monolithic_solution  #src
feasible_inplace_solution == monolithic_solution

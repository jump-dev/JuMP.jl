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
import Printf
import Test  #src

# ## Theory

# Benders decomposition is a useful algorithm for solving convex optimization
# problems with a large number of variables. It works best when a larger problem
# can be decomposed into two (or more) smaller problems that are individually
# much easier to solve. This tutorial demonstrates Benders decomposition on the
# following mixed-integer linear program:
# ```math
# \begin{aligned}
# \text{min}        \ & c_1^\top x+c_2^\top y   \\
# \text{subject to} \ & A_1 x+ A_2 y \le b      \\
#                     & x \ge 0                 \\
#                     & y \ge 0                 \\
#                     & x \in \mathbb{Z}^n
# \end{aligned}
# ```
# where $b \in \mathbb{R}^m$, $A_1 \in \mathbb{R}^{m \times n}$,
# $A_2 \in \mathbb{R}^{m \times p}$ and $\mathbb{Z}$ is the set of integers.
#
# Any mixed integer programming problem can be written in the form above.

# If there are relatively few integer variables, and many more continuous
# variables, then it may be beneficial to decompose the problem into a small
# problem containing only integer variables and a linear program containing
# only continuous variables. Hopefully, the linear program will be much easier
# to solve in isolation than in the full mixed-integer linear program.

# For example, if we knew a feasible solution for ``x``, we could obtain a
# solution for ``y`` by solving:
# ```math
# \begin{aligned}
# V_2(x) = & \text{min}          \ & c_2^\top y                        \\
#          & \text{subject to}   \ & A_2 y \le b - A_1 x & \quad [\pi] \\
#          &                       & y \ge 0,
# \end{aligned}
# ```
# where ``\pi`` is the dual variable associated with the constraints. Because
# this is a linear program, it is easy to solve.

# Replacing the ``c_2^\top y`` component of the objective in our original
# problem with ``V_2`` yields:
# ```math
# \begin{aligned}
# \text{min}        \ & c_1^\top x + V_2(x) \\
# \text{subject to} \ & x \ge 0             \\
#                     & x \in \mathbb{Z}^n
# \end{aligned}
# ```
# This problem looks a lot simpler to solve, but we need to do something else
# with ``V_2`` first.

# Because ``x`` is a constant that appears on the right-hand side of the
# constraints, ``V_2`` is a convex function with respect to ``x``, and the dual
# variable ``\pi`` can be multiplied by ``-A_1`` to obtain a subgradient of
# ``V_2(x)`` with respect to ``x``. Therefore, if we have a candidate solution
# ``x_k``, then we can solve ``V_2(x_k)`` and obtain a feasible dual vector
# ``\pi_k``. Using these values, we can construct a first-order Taylor-series
# approximation of ``V_2`` about the point ``x_k``:
# ```math
# V_2(x) \ge V_2(x_k) + -\pi_k^\top A_1 (x - x_k).
# ```
# By convexity, we know that this inequality holds for all ``x``, and we call
# these inequalities _cuts_.

# Benders decomposition is an iterative technique that replaces ``V_2(x)`` with
# a new decision variable ``\theta``, and approximates it from below using cuts:
# ```math
# \begin{aligned}
# V_1^K = & \text{min}         \ & c_1^\top x + \theta  \\
#         &  \text{subject to} \ & x \ge 0              \\
#         &                    \ & x \in \mathbb{Z}^n   \\
#         &                    \ & \theta \ge M         \\
#         &                    \ & \theta \ge V_2(x_k) + \pi_k^\top(x - x_k) & \quad \forall k = 1,\ldots,K.
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
# solution for ``x`` from the first-stage problem, then ``c_1^\top x + V_2(x)``
# is a valid _upper_ bound on the objective value of the full problem.

# Benders decomposition uses the lower and upper bounds to determine when it has
# found the global optimal solution.

# ## Input data

# As an example for this tutorial, we use the input data is from page 139 of
# Garfinkel, R. & Nemhauser, G. L. Integer programming. (Wiley, 1972).

c_1 = [1, 4]
c_2 = [2, 3]
dim_x = length(c_1)
dim_y = length(c_2)
b = [-2; -3]
A_1 = [1 -3; -1 -3]
A_2 = [1 -2; -1 -1]
M = -1000;

# ## [Iterative method](@id benders_iterative)

# !!! warning
#     This is a basic implementation for pedagogical purposes. We haven't
#     discussed Benders feasibility cuts, or any of the computational tricks
#     that are required to build a performant implementation for large-scale
#     problems. See [In-place iterative method](@ref) for one improvement that
#     helps computation time.

# We start by formulating the first-stage subproblem:

model = Model(GLPK.Optimizer)
@variable(model, x[1:dim_x] >= 0, Int)
@variable(model, θ >= M)
@objective(model, Min, c_1' * x + θ)
print(model)

# For the next step, we need a function that takes a first-stage candidate
# solution `x` and returns the optimal solution from the second-stage
# subproblem:

function solve_subproblem(x)
    model = Model(GLPK.Optimizer)
    @variable(model, y[1:dim_y] >= 0)
    con = @constraint(model, A_2 * y .<= b - A_1 * x)
    @objective(model, Min, c_2' * y)
    optimize!(model)
    @assert termination_status(model) == OPTIMAL
    return (obj = objective_value(model), y = value.(y), π = dual.(con))
end

# Note that `solve_subproblem` returns a `NamedTuple` of the objective value,
# the optimal primal solution for `y`, and the optimal dual solution for `π`.

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
    lower_bound = objective_value(model)
    x_k = value.(x)
    ret = solve_subproblem(x_k)
    upper_bound = c_1' * x_k + ret.obj
    gap = (upper_bound - lower_bound) / upper_bound
    print_iteration(k, lower_bound, upper_bound, gap)
    if gap < ABSOLUTE_OPTIMALITY_GAP
        println("Terminating with the optimal solution")
        break
    end
    cut = @constraint(model, θ >= ret.obj + -ret.π' * A_1 * (x .- x_k))
    @info "Adding the cut $(cut)"
end

# Finally, we can obtain the optimal solution

optimize!(model)
Test.@test value.(x) == [0.0, 1.0]  #src
x_optimal = value.(x)

#-

optimal_ret = solve_subproblem(x_optimal)
Test.@test optimal_ret.y == [0.0, 0.0]  #src
y_optimal = optimal_ret.y

# ## Callback method

# The [Iterative method](@ref benders_iterative) section implemented Benders
# decomposition using a loop. In each iteration, we re-solved the first-stage
# subproblem to generate a candidate solution. However, modern MILP solvers such
# as CPLEX, Gurobi, and GLPK provide lazy constraint callbacks which allow us to
# add new cuts _while the solver is running_. This can be more efficient than an
# iterative method because we can avoid repeating work such as solving the root
# node of the first-stage MILP at each iteration.

# !!! tip
#     For more information on callbacks, read the page
#     [Solver-independent callbacks](@ref callbacks_manual).

# As before, we construct the same first-stage subproblem:

lazy_model = Model(GLPK.Optimizer)
@variable(lazy_model, x[1:dim_x] >= 0, Int)
@variable(lazy_model, θ >= M)
@objective(lazy_model, Min, θ)
print(lazy_model)

# What differs is that we write a callback function instead of a loop:

k = 0

"""
    my_callback(cb_data)

A callback that implements Benders decomposition. Note how similar it is to the
inner loop of the iterative method.
"""
function my_callback(cb_data)
    global k += 1
    x_k = callback_value.(cb_data, x)
    θ_k = callback_value(cb_data, θ)
    lower_bound = c_1' * x_k + θ_k
    ret = solve_subproblem(x_k)
    upper_bound = c_1' * x_k + c_2' * ret.y
    gap = (upper_bound - lower_bound) / upper_bound
    print_iteration(k, lower_bound, upper_bound, gap)
    if gap < ABSOLUTE_OPTIMALITY_GAP
        println("Terminating with the optimal solution")
        return
    end
    cut = @build_constraint(θ >= ret.obj + -ret.π' * A_1 * (x .- x_k))
    MOI.submit(model, MOI.LazyConstraint(cb_data), cut)
    return
end

set_attribute(lazy_model, MOI.LazyConstraintCallback(), my_callback)

# Now when we optimize!, our callback is run:

optimize!(lazy_model)

# Note how this problem also takes 4 iterations to converge, but the sequence
# of bounds is different compared to the iterative method.

# Finally, we can obtain the optimal solution:

Test.@test value.(x) == [0.0, 1.0]  #src
x_optimal = value.(x)

#-

optimal_ret = solve_subproblem(x_optimal)
Test.@test optimal_ret.y == [0.0, 0.0]  #src
y_optimal = optimal_ret.y

# ## In-place iterative method

# Our implementation of the iterative method has a problem: every time we need
# to solve the subproblem, we must rebuild it from scratch. This is expensive,
# and it can be the bottleneck in the solution process. We can improve our
# implementation by using re-using the subproblem between solves.

# First, we create our first-stage problem as usual:

model = Model(GLPK.Optimizer)
@variable(model, x[1:dim_x] >= 0, Int)
@variable(model, θ >= M)
@objective(model, Min, c_1' * x + θ)
print(model)

# Then, instead of building the subproblem in a function, we build it once here:

subproblem = Model(GLPK.Optimizer)
@variable(subproblem, x_copy[1:dim_x])
@variable(subproblem, y[1:dim_y] >= 0)
@constraint(subproblem, A_1 * x_copy + A_2 * y .<= b)
@objective(subproblem, Min, c_2' * y)
print(subproblem)

# This formulation is slightly different. We have included a copy of the x
# variables, `x_copy`, and used `x_copy` in the left-hand side of the
# constraints.

# Our function to solve the subproblem is also slightly different. First, we
# need to fix the value of the `x_copy` variables to the value of `x` from the
# first-stage problem, and second, we compute the dual using the
# [`reduced_cost`](@ref) of `x_copy`, not the dual of the linear constraints:

function solve_subproblem(model, x)
    fix.(model[:x_copy], x)
    optimize!(model)
    @assert termination_status(model) == OPTIMAL
    return (
        obj = objective_value(model),
        y = value.(model[:y]),
        π = reduced_cost.(model[:x_copy]),
    )
end

# Now we're ready to iterate our in-place Benders decomposition, but this time
# the cut computation is slightly different. Because we used
# [`reduced_cost`](@ref) on the `x_copy` variables, the value of `π` is a valid
# subgradient on the objective of `subproblem` with respect to `x`. Therefore,
# we don't need to multiply the duals by `-A_1`, and so our cut uses `ret.π'`
# instead of `-ret.π' * A_1`:

println("Iteration  Lower Bound  Upper Bound          Gap")
for k in 1:MAXIMUM_ITERATIONS
    optimize!(model)
    lower_bound = objective_value(model)
    x_k = value.(x)
    ret = solve_subproblem(subproblem, x_k)
    upper_bound = c_1' * x_k + ret.obj
    gap = (upper_bound - lower_bound) / upper_bound
    print_iteration(k, lower_bound, upper_bound, gap)
    if gap < ABSOLUTE_OPTIMALITY_GAP
        println("Terminating with the optimal solution")
        break
    end
    cut = @constraint(model, θ >= ret.obj + ret.π' * (x .- x_k))
    @info "Adding the cut $(cut)"
end

# Finally, we can obtain the optimal solution:

optimize!(model)
Test.@test value.(x) == [0.0, 1.0]  #src
x_optimal = value.(x)

#-

optimal_ret = solve_subproblem(subproblem, x_optimal)
Test.@test optimal_ret.y == [0.0, 0.0]  #src
y_optimal = optimal_ret.y

# For larger problems, the benefit of re-using the same subproblem and not
# needing to multiply the duals by `A_1` in the cut coefficient usually
# outweights the cost of needing a copy of the `x` variables in the subproblem.

# As a secondary benefit, because we no longer need an explicit representation
# of `A_1` in the cut, we can build the `model` and `subproblem` formulations
# using arbitrary JuMP syntax; they do not need to be in matrix form.

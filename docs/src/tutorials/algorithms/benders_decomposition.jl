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

# **Originally Contributed by**: Shuvomoy Das Gupta

# This tutorial describes how to implement [Benders decomposition](https://en.wikipedia.org/wiki/Benders_decomposition)
# in JuMP using an iterative method. It uses the following packages:

using JuMP
import GLPK
import Printf
import Test  #src

# !!! info
#     For an alternative approach using callbacks, read
#     [Benders decomposition (via callbacks)](@ref benders_decomposition_lazy).

# Benders decomposition is a useful algorithm for solving convex optimization
# problems with a large number of variables. It works best when a larger problem
# can be decomposed into two (or more) smaller problems that are invidually much
# easier to solve. This tutorial demonstrates Benders decomposition on the
# following mixed-integer linear program:
# ```math
# \begin{aligned}
# \text{min} \quad &&c_1^\top x+c_2^\top y \\
# \text{subject to} \quad &&A_1 x+ A_2 y \le b \\
# & x \ge 0 \\
# & y \ge 0 \\
# & x \in \mathbb{Z}^n
# \end{aligned}
# ```
# where $b \in \mathbb{R}^m$, $A_1 \in \mathbb{R}^{m \times n}$,
# $A_2 \in \mathbb{R}^{m \times p}$ and $\mathbb{Z}$ is the set of integers.
#
# Any mixed integer programming problem can be written in the form above.

# If there are relatively few integer variables, and many more continuous
# variables, then it may be beneficial to decompose the problem into a small
# problem containing only integer variables, and a linear program containing
# only continous variables. Hopefully, the linear program will be much easier to
# solve in isolation that in the full mixed-integer linear program.

# For example, if we knew a feasible solution for ``x``, we could obtain a
# solution for ``y`` by solving:
# ```math
# \begin{aligned}
# V_2(x) = \text{min} \quad && c_2^\top y \\
# \text{subject to} \quad && A_2 y \le b - A_1 x & [\pi] \\
# & y \ge 0.
# \end{aligned}
# ```
# Since this is a linear program it is easy to solve.

# Replacing the ``c_2^\top y`` component of the objective in our original
# problem with ``V_2`` yields:
# ```math
# \begin{aligned}
# \text{min} \quad && c_1^\top x + V_2(x) \\
# \text{subject to} & x \ge 0 \\
# & x \in \mathbb{Z}^n
# \end{aligned}
# ```
# This problem looks a lot simpler to solve, but we need to do something else
# with ``V_2`` first.

# Because ``x`` is a constant that appears on the right-hand side of the
# constraints, ``V_2`` is a convex function with respect to ``x``, and the dual
# variable ``\pi`` is a subgradient of ``V_2(x)``. Therefore, if we have some
# candidate solution ``x_k``, then we can solve ``V_2(x_k)`` and obtain a
# feasible dual vector ``\pi_k``. Using these values, we can constraint a
# first-order Taylor-series approximation of ``V_2`` about the point ``x_k``:
# ```math
# V_2(x) \ge V_2(x_k) + -\pi_k^\top A_1 (x - x_k).
# ```
# By convexity, we know that this inequality holds for all ``x``, and we call
# these inequalities _cuts_.

# Benders decomposition is an iterative technique that replaces ``V_2(x)`` with
# a new decision variable ``theta``, and approximates it from below using cuts:
# ```math
# \begin{aligned}
# V_1^K = \text{min} \quad && c_1^\top x + \theta \\
# \text{subject to} & x \ge 0 \\
# & x \in \mathbb{Z}^n \\
# & \theta \le M \\
# & \theta \le V_2(x_k) + \pi_k^\top(x - x_k) & \forall k = 1,\ldots,K.
# \end{aligned}
# ```
# This optimization is called the _first-stage_ subproblem.

# ## Input data

# The input data is from page 139 of Garfinkel, R. & Nemhauser, G. L. Integer
# programming. (Wiley, 1972).

c_1 = [1, 4]
c_2 = [2, 3]
dim_x = length(c_1)
dim_y = length(c_2)
b = [-2; -3]
A_1 = [1 -3; -1 -3]
A_2 = [1 -2; -1 -1]
M = -1000;

# ## Implementation

# !!! warning
#     This is a basic implementation for pedagogical purposes. We haven't
#     discussed Benders feasibility cuts, or any of the computational tricks
#     that are required to build a performative implementation for large-scale
#     problems.

# Here is the first-stage subproblem:

model = Model(GLPK.Optimizer)
@variable(model, x[1:dim_x] >= 0, Int)
@variable(model, θ)
@constraint(model, c_1' * x + θ >= M)
@objective(model, Min, c_1' * x + θ)
print(model)

# Here is a function that takes a first-stage candidate solution `x` and returns
# the optimal solution from the second-stage subproblem:

function solve_subproblem(x)
    model = Model(GLPK.Optimizer)
    @variable(model, y[1:dim_y] >= 0)
    con = @constraint(model, A_2 * y .<= b - A_1 * x)
    @objective(model, Min, c_2' * y)
    optimize!(model)
    @assert termination_status(model) == OPTIMAL
    return (obj = objective_value(model), y = value.(y), π = dual.(con))
end

# We're almost ready for our optimization look, but first, here's a helpful
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
    upper_bound = objective_value(model) - value(θ) + ret.obj
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

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

# # [Benders decomposition (via callbacks)](@id benders_decomposition_lazy)

# **Originally Contributed by**: Mathieu Besançon

# This notebook describes how to implement the Benders decomposition in JuMP
# using lazy constraints. We keep the same notation and problem form as the
# first notebook *Benders decomposition*.

# ```math
# \begin{aligned}
# & \text{maximize} \quad &&c_1^T x+c_2^T v \\
# & \text{subject to} \quad &&A_1 x+ A_2 v \preceq b \\
# & &&x \succeq 0, x \in \mathbb{Z}^n \\
# & &&v \succeq 0, v \in \mathbb{R}^p \\
# \end{aligned}
# ```

# where $b \in \mathbb{R}^m$, $A_1 \in \mathbb{R}^{m \times n}$,
# $A_2 \in \mathbb{R}^{m \times p}$ and $\mathbb{Z}$ is the set of integers.
#
# Here the symbol $\succeq$ ($\preceq$) stands for element-wise greater (less)
# than or equal to. Any mixed integer programming problem can be written in the
# form above.

# For a detailed explanation on the Benders decomposition algorithm, see the
# introduction notebook.

# ## Lazy constraints

# Some optimization solvers allow users to interact with them during the
# solution process by providing user-defined functions which are triggered under
# certain conditions.

# The generic term for these functions is **callback**.

# In integer optimization, the main callback types are lazy constraints,
# user-cuts and heuristic solutions. See the [Callbacks](@ref callbacks_manual)
# section for an introduction on how to use them.
#
# Some callbacks define a new constraint which is only activated when necessary,
# i.e., when a current solution does not respect them. It can avoid building an
# optimization model with too many constraints up-front.
#
# This is the case for Benders decomposition, since the sub-problem defines an
# exponential number of primal vertices and therefore dual cuts.
#
# A detailed explanation on the distinction between user-cuts and lazy
# constraints is also available on [Paul Rubin's blog](https://orinanobworld.blogspot.com/2012/08/user-cuts-versus-lazy-constraints.html).
# He also [describes this approach to Benders Decomposition](https://orinanobworld.blogspot.ca/2011/10/benders-decomposition-then-and-now.html).

# We use the data from the original notebook and change the solution algorithm
# to leverage lazy constraints:

# * Step 1 (Initialization)
# * Step 2 (defining the subproblem model)
# * Step 3 (registering the lazy constraint of the subproblem)

# ## Input data

# The input data is from page 139, Integer programming by Garfinkel and
# Nemhauser[[1]](#c1).

c1 = [-1; -4]
c2 = [-2; -3]

dim_x = length(c1)
dim_u = length(c2)

b = [-2; -3]

A1 = [
    1 -3
    -1 -3
]
A2 = [
    1 -2
    -1 -1
]

M = 1000;

# ## Implementation

using JuMP
import GLPK
import Test

# !!! info
#     This tutorial uses the [MathOptInterface](@ref moi_documentation) API.
#     By default, JuMP exports the `MOI` symbol as an alias for the
#     MathOptInterface.jl package. We recommend making this more explicit in
#     your code by adding the following lines:
#     ```julia
#     import MathOptInterface
#     const MOI = MathOptInterface
#     ```

# Subproblem creation

function build_subproblem()
    sub_problem_model = Model(GLPK.Optimizer)
    @variable(sub_problem_model, u[1:dim_u] >= 0)
    @constraint(
        sub_problem_model,
        constr_ref_subproblem[j = 1:size(A2, 2)],
        A2[:, j]' * u >= c2[j],
    )
    return (sub_problem_model, u)
end

# ### Master Problem Description

master_problem_model = Model(GLPK.Optimizer)
@variable(master_problem_model, 0 <= x[1:dim_x] <= 1e6, Int)
@variable(master_problem_model, t <= 1e6)
@objective(master_problem_model, Max, t)
print(master_problem_model)

#-

(sub_problem_model, u) = build_subproblem();

# Track the calls to the callback

iter_num = 0

# Define lazy constraints

function benders_lazy_constraint_callback(cb_data)
    global iter_num
    iter_num += 1
    println("Iteration number = ", iter_num)

    x_current = callback_value.(Ref(cb_data), x)
    fm_current = callback_value(cb_data, t)

    c_sub = b - A1 * x_current
    @objective(sub_problem_model, Min, c1' * x_current + c_sub' * u)
    optimize!(sub_problem_model)

    print("\nThe current subproblem model is \n")
    print(sub_problem_model)

    t_status_sub = termination_status(sub_problem_model)
    p_status_sub = primal_status(sub_problem_model)

    fs_x_current = objective_value(sub_problem_model)

    u_current = value.(u)

    γ = b' * u_current

    if p_status_sub == FEASIBLE_POINT && fs_x_current ≈ fm_current # we are done
        @info("No additional constraint from the subproblem")
    end

    if p_status_sub == FEASIBLE_POINT && fs_x_current < fm_current
        println(
            "\nThere is a suboptimal vertex, add the corresponding constraint",
        )
        cv = A1' * u_current - c1
        new_optimality_cons = @build_constraint(t + cv' * x <= γ)
        MOI.submit(
            master_problem_model,
            MOI.LazyConstraint(cb_data),
            new_optimality_cons,
        )
    end

    if t_status_sub == INFEASIBLE_OR_UNBOUNDED
        println(
            "\nThere is an  extreme ray, adding the corresponding constraint",
        )
        ce = A1' * u_current
        new_feasibility_cons = @build_constraint(dot(ce, x) <= γ)
        MOI.submit(
            master_problem_model,
            MOI.LazyConstraint(cb_data),
            new_feasibility_cons,
        )
    end
end

MOI.set(
    master_problem_model,
    MOI.LazyConstraintCallback(),
    benders_lazy_constraint_callback,
)

optimize!(master_problem_model)

t_status = termination_status(master_problem_model)
p_status = primal_status(master_problem_model)

if p_status == INFEASIBLE_POINT
    println("The problem is infeasible :-(")
end

if t_status == INFEASIBLE_OR_UNBOUNDED
    fm_current = M
    x_current = M * ones(dim_x)
end

if p_status == FEASIBLE_POINT
    fm_current = value(t)
    x_current = value.(x)
end

Test.@test value(t) ≈ -4 #hide

println(
    "Status of the master problem is ",
    t_status,
    "\nwith fm_current = ",
    fm_current,
    "\nx_current = ",
    x_current,
)

# ### References

# ```@raw html
# <a id='c1'></a>
# ```

# 1. Garfinkel, R. & Nemhauser, G. L. Integer programming. (Wiley, 1972).

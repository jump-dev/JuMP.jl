# Copyright 2017, Iain Dunning, Joey Huchette, Miles Lubin, and contributors    #src
# This Source Code Form is subject to the terms of the Mozilla Public License   #src
# v.2.0. If a copy of the MPL was not distributed with this file, You can       #src
# obtain one at https://mozilla.org/MPL/2.0/.                                   #src

# # Dualization

# The purpose of this tutorial is to explain how to use [Dualization.jl](@ref) to
# improve the performance of some conic optimization models. There are two
# important takeaways:
#
#  1. JuMP reformulates problems to meet the input requirements of the
#     solver, potentially increasing the problem size by adding slack variables
#     and constraints.
#  2. Solving the dual of a conic model can be more efficient than solving the
#     primal.

# This tutorial uses the following packages

using JuMP
import Dualization
import SCS

# ## Background

# Conic optimization solvers typically accept one of two input formulations.

# The first is the _standard_ conic form:
# ```math
# \begin{align}
#     \min_{x \in \mathbb{R}^n} \; & c^\top x \\
#               \;\;\text{s.t.} \; & A x = b  \\
#                               & x \in \mathcal{K}
# \end{align}
# ```
# in which we have a set of linear equality constraints $Ax = b$ and the
# variables belong to a cone $\mathcal{K}$.
#
# The second is the _geometric_ conic form:
# ```math
# \begin{align}
#     \min_{x \in \mathbb{R}^n} \; & c^\top x \\
#               \;\;\text{s.t.} \; & A x - b \in \mathcal{K}
# \end{align}
# ```
# in which an affine function $Ax - b$ belongs to a cone $\mathcal{K}$ and the
# variables are free.

# It is trivial to convert between these two representations, for example, to go
# from the geometric conic form to the standard conic form we introduce slack
# variables $y$:
# ```math
# \begin{align}
#     \min_{x \in \mathbb{R}^n} \; & c^\top x \\
#               \;\;\text{s.t.} \; & \begin{bmatrix}A & -I\end{bmatrix} [x; y] = b \\
#                               & [x; y] \in \mathbb{R}^n \times \mathcal{K}
# \end{align}
# ```
# and to go from the standard conic form to the geometric conic form, we can
# rewrite the equality constraint as a function belonging to the `{0}` cone:
# ```math
# \begin{align}
#     \min_{x \in \mathbb{R}^n} \; & c^\top x \\
#               \;\;\text{s.t.} & [A; I] x - [b; 0] \in \{0\} \times \mathcal{K}
# \end{align}
# ```

# From a theoretical perspective, the two formulations are equivalent, and if
# you implement a model in the standard conic form and pass it to a geometric
# conic form solver (or vice versa), then JuMP will automatically reformulate
# the problem into the correct formulation.

# From a practical perspective though, the geometric to standard reformulation is
# particularly problematic, because the additional slack variables and
# constraints can make the problem much larger and therefore harder to solve.

# You should also note many problems contain a mix of conic constraints and
# variables, and so they do not neatly fall into one of the two formulations. In
# these cases, JuMP reformulates only the variables and constraints as necessary
# to convert the problem into the desired form.

# ## Primal and dual formulations

# Duality plays a large role in conic optimization. For a detailed description
# of conic duality, see [Duality](@ref).

# A useful observation is that if the primal problem is in standard conic form,
# then the dual problem is in geometric conic form, and vice versa. Moreover, the
# primal and dual may have a different number of variables and constraints,
# although which one is smaller depends on the problem. Therefore, instead of
# reformulating the problem from one form to the other, it can be more
# efficient to solve the dual instead of the primal.

# To demonstrate, we use a variation of the [Maximum cut via SDP](@ref) example.

# The primal formulation (in standard conic form) is:

model_primal = Model()
@variable(model_primal, X[1:2, 1:2], PSD)
@objective(model_primal, Max, sum([1 -1; -1 1] .* X))
@constraint(model_primal, primal_c[i = 1:2], 1 - X[i, i] == 0)
print(model_primal)

# This problem has three scalar decision variables (the matrix `X` is symmetric),
# two scalar equality constraints, and a constraint that `X` is positive
# semidefinite.

# The dual of `model_primal` is:

model_dual = Model()
@variable(model_dual, y[1:2])
@objective(model_dual, Min, sum(y))
@constraint(model_dual, dual_c, [y[1]-1 1; 1 y[2]-1] in PSDCone())
print(model_dual)

# This problem has two scalar decision variables, and a 2x2 positive
# semidefinite matrix constraint.

# !!! tip
#     If you haven't seen conic duality before, try deriving the dual problem
#     based on the description in [Duality](@ref). You'll need to know that the
#     dual cone of [`PSDCone`](@ref) is the [`PSDCone`](@ref).

# When we solve `model_primal` with `SCS.Optimizer`, SCS reports three variables
# (`variables n: 3`), five rows in the constraint matrix (`constraints m: 5`),
# and five non-zeros in the matrix (`nnz(A): 5`):

set_optimizer(model_primal, SCS.Optimizer)
optimize!(model_primal)
Base.Libc.flush_cstdio()  #hide

# (There are five rows in the constraint matrix because SCS expects problems in
# geometric conic form, and so JuMP has reformulated the `X, PSD` variable
# constraint into the affine constraint `X .+ 0 in PSDCone()`.)

# The solution we obtain is:

value.(X)

#-

dual.(primal_c)

#-

objective_value(model_primal)

# When we solve `model_dual` with `SCS.Optimizer`, SCS reports two variables
# (`variables n: 2`), three rows in the constraint matrix (`constraints m: 3`),
# and two non-zeros in the matrix (`nnz(A): 2`):

set_optimizer(model_dual, SCS.Optimizer)
optimize!(model_dual)
Base.Libc.flush_cstdio()  #hide

# and the solution we obtain is:

dual.(dual_c)

#-

value.(y)

#-

objective_value(model_dual)

# The problems are small enough that it isn't meaningful to compare the solve
# times, but in general, we should expect the dual problem to solve faster
# because it contains fewer variables and constraints. The difference is
# particularly noticeable on large-scale optimization problems.

# ## `dual_optimizer`

# Manually deriving the conic dual is difficult and error-prone. The package
# [Dualization.jl](@ref) provides the `Dualization.dual_optimizer` meta-solver,
# which wraps any MathOptInterface-compatible solver in an interface that
# automatically formulates the dual of an input problem, solves the dual
# problem, and then reports the primal solution back to the user.

# To demonstrate, we use `Dualization.dual_optimizer` to solve `model_primal`:

set_optimizer(model_primal, Dualization.dual_optimizer(SCS.Optimizer))
optimize!(model_primal)
Base.Libc.flush_cstdio()  #hide

# The performance is the same as if we solved `model_dual`, and the correct
# solution is returned to `X`:

value.(X)

#-

dual.(primal_c)

# Moreover, if we use `dual_optimizer` on `model_dual`, then we get the same
# performance as if we had solved `model_primal`:

set_optimizer(model_dual, Dualization.dual_optimizer(SCS.Optimizer))
optimize!(model_dual)
Base.Libc.flush_cstdio()  #hide

#-

dual.(dual_c)

#-

value.(y)

# ## When to use `dual_optimizer`

# Because it can make the problem larger or smaller, depending on the problem
# and the choice of solver, there is no definitive rule on when you should
# use `dual_optimizer`. However, you should try `dual_optimizer` if your conic
# optimization problem takes a long time to solve, or if you need to repeatedly
# solve similarly structured problems with different data. In some cases solving
# the dual instead of the primal can make a large difference.

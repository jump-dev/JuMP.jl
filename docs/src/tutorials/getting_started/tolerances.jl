# Copyright (c) 2024 Oscar Dowson and contributors                               #src
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

# # Tolerances and numerical issues

# Optimization solvers can seem like magic black boxes that take in an algebraic
# formulation of a problem and return a solution. It is tempting to treat their
# solutions at face value, since we often have little ability to verify that the
# solution is in fact optimal. However, like all numerical algorithms that use
# floating point arithmetic, optimization solvers use tolerances to check
# whether a solution satisfies the constraints. In the best case, the solution
# satisfies the original constraints to
# [machine precision](https://en.wikipedia.org/wiki/Machine_epsilon). In most
# cases, the solution satisfies the constraints to some very small tolerance
# that has no noticeable impact on the quality of the optimal solution. In the
# worst case, the solver can return a "wrong" solution, or fail to find one even
# if it exists. (In the last case, the solution is wrong only in the sense of
# user expectation. It will satisfy the solution to the tolerances that are
# provided.)

# The purpose of this tutorial is to explain the various types of tolerances
# that are used in optimization solvers and what you can reasonably expect from
# a solution.

# There are a few sources of additional information:
#  * Ambros Gleixner has an excellent YouTube talk
#    [Numerics in LP & MIP Solvers](https://youtu.be/rKcdF4Fgl-g?feature=shared).
#  * Gurobi has a series of articles in their documentation called
#    [Guidelines for Numerical Issues](https://www.gurobi.com/documentation/current/refman/guidelines_for_numerical_i.html)

# !!! tip
#     This tutorial is more advanced than the other "Getting started" tutorials.
#     It's in the "Getting started" section to give you an early preview of how
#     tolerances affect JuMP models. However, if you are new to JuMP, you may
#     want to briefly skim the tutorial, and come back to it once you have
#     written a few JuMP models.

# ## Required packages

# This tutorial uses the following packages:

using JuMP
import HiGHS
import SCS

# ## Background

# Optimization solvers use tolerances to check the feasibility of constraints.

# There are four main types of tolerances:
#
# 1. primal feasibility: controls how feasibility of the primal solution is
#    measured
# 2. dual feasibility: controls how feasibility of the dual solution is measured
# 3. integrality: controls how feasibility of the binary and integer variables
#    are measured
# 4. optimality: controls how close the primal and dual solutions must be.
#
# Solvers may use absolute tolerances, relative tolerances, or some mixture of
# both. The definition and default value of each tolerance is solver-dependent.

# The dual feasibility tolerance is much the same as the primal feasibility
# tolerance, only that operates on the space of dual solutions instead of the
# primal. HiGHS has `dual_feasibility_tolerance`, but some solvers have only a
# single feasibility tolerance that uses the same value for both.

# The optimality tolerance is a more technical tolerance that is used to test
# the equivalence of the primal and dual objectives in the KKT system if you are
# solving a continuous problem via interior point. HiGHS has
# `ipm_optimality_tolerance`, but some solvers will not have such a tolerance.
# Note that the optimality tolerance is different to the relative MIP gap that
# controls early termination of a MIP solution during branch-and-bound.

# Because the dual and optimality tolerances are less used, this tutorial
# focuses on the primal feasibility and integrality tolerances.

# ## Primal feasibility

# The primal feasibility tolerance controls how primal constraints are
# evaluated. For example, the constraint ``2x = 1`` is actually implemented as
# ``|2x - 1| \le \varepsilon``, where ``\varepsilon`` is a small
# solver-dependent primal feasibility tolerance that is typically on the order
# of `1e-8`.

# Here's an example in practice. This model should be infeasible, since `x` must
# be non-negative, but there is also an equality constraint that `x` is equal to
# a small negative number. Yet when we solve this problem, we get:

model = Model(HiGHS.Optimizer)
set_silent(model)
@variable(model, x >= 0)
@constraint(model, x == -1e-8)
optimize!(model)
@assert is_solved_and_feasible(model)  #src
is_solved_and_feasible(model)

#-

value(x)

# In other words, HiGHS thinks that the solution `x = 0` satisfies the
# constraint `x == -1e-8`. The value of ``\varepsilon`` in HiGHS is controlled
# by the `primal_feasibility_tolerance` option. If we set this to a smaller
# value, HiGHS will now correctly deduce that the problem is infeasible:

set_attribute(model, "primal_feasibility_tolerance", 1e-10)
optimize!(model)
@assert !is_solved_and_feasible(model)  #src
is_solved_and_feasible(model)

# ### Realistic example

# Here's a more realistic example, which was reported in the [SCS.jl](https://github.com/jump-dev/SCS.jl/issues/297)
# repository:

n, ε = 13, 0.0234
N = 2^n
model = Model(SCS.Optimizer)
@variable(model, x[1:N] >= 0)
@objective(model, Min, x[1])
@constraint(model, sum(x) == 1)
z = [(-1)^((i & (1 << j)) >> j) for j in 0:n-1, i in 0:N-1]
@constraint(model, z * x .>= 1 - ε)
optimize!(model)

# SCS reports that it solved the problem to optimality:

@assert is_solved_and_feasible(model)  #src
is_solved_and_feasible(model)

# and that the solution for `x[1]` is nearly zero:

value(x[1])

# However, the analytic solution for `x[1]` is:

1 - n * ε / 2

# The answer is very wrong, and there is no indication from the solver that
# anything untoward happened. What's going on?

# One useful debugging tool is [`primal_feasibility_report`](@ref):

report = primal_feasibility_report(model)

# `report` is a dictionary which maps constraints to the violation.  The largest
# violation is approximately `1e-5`:

maximum(values(report))

# This makes sense, because the default primal feasibility tolerance for SCS is
# `1e-4`.

# Most of the entries are lower bound constraints on the variables. Here are all
# the variables which violate their lower bound:

violated_variables = filter(xi -> value(xi) < 0, x)

# The first one is:

y = first(violated_variables)

# It has a primal value of:

value(y)

# which matches the value in the feasibility report:

report[LowerBoundRef(y)]

# Despite the small primal feasibility tolerance and the small actual violations
# of the constraints, our optimal solution is very far from the theoretical
# optimal.

# We can "fix" our model by decreasing `eps_abs` and `eps_rel`, which SCS uses
# to control the absolute and relative feasibility tolerances. Now the solver
# finds the correct solution:

set_attribute(model, "eps_abs", 1e-5)
set_attribute(model, "eps_rel", 1e-5)
optimize!(model)

#-

@assert is_solved_and_feasible(model)  #src
@assert is_solved_and_feasible(model)
value(x[1])

# ### Why you shouldn't use a small tolerance

# There is no direct relationship between the size of feasibility tolerance and
# the quality of the solution.

# You might surmise from this section that you should set the tightest
# feasibility tolerance possible. However, tighter tolerances come at the cost
# of increased solve time.

# For example, SCS is a first-order solver. This means it uses only local
# gradient information at update each iteration. SCS took 100 iterations to
# solve the problem with the default tolerance of `1e-4`, and 550 iterations to
# solve the problem with `1e-5`. SCS may not be able to find a solution to our
# problem with a tighter tolerance in a reasonable amount of time.

# ## Integrality

# Integrality tolerances control how the solver decides if a variable satisfies
# an integrality of binary constraint. The tolerance is typically defined as:
# ``|x - \lfloor x + 0.5 \rfloor | \le \varepsilon``, which you can read as the
# absolute distance to the nearest integer.

# Here's a simple example:

model = Model(HiGHS.Optimizer)
set_silent(model)
@variable(model, x == 1 + 1e-6, Int)
optimize!(model)
@assert is_solved_and_feasible(model)  #src
is_solved_and_feasible(model)

# HiGHS found an optimal solution, and the value of `x` is:

@assert isapprox(value(x), 1.000001)  #src
value(x)

# In other words, HiGHS thinks that the solution `x = 1.000001` satisfies the
# constraint that `x` must be an integer. [`primal_feasibility_report`](@ref)
# shows that indeed, the integrality constraint is violated:

primal_feasibility_report(model)

# The value of ``\varepsilon`` in HiGHS is controlled by the
# `mip_feasibility_tolerance` option. The default is `1e-6`. If we set the
# attribute to a smaller value, HiGHS will now correctly deduce that the problem
# is infeasible:

set_attribute(model, "mip_feasibility_tolerance", 1e-10)
optimize!(model)
@assert !is_solved_and_feasible(model)  #src
is_solved_and_feasible(model)

# ### Realistic example

# Integrality tolerances are particularly important when you have big-M type
# constraints. Small non-integer values in the integer variables can cause
# "leakage" flows even when the big-M switch is "off." Consider this example:

M = 1e6
model = Model()
@variable(model, x >= 0)
@variable(model, y, Bin)
@constraint(model, x <= M * y)
print(model)

# This model has a feasible solution (to tolerances) of `(x, y) = (1, 1e-6)`.
# There can be a non-zero value of `x` even when `y` is (approximately) `0`.

primal_feasibility_report(model, Dict(x => 1.0, y => 1e-6))

# ### Rounding the solution

# Integrality tolerances are why JuMP does not return `Int` for `value(x)` of an
# integer variable or `Bool` for `value(x)` of a binary variable.

# In most cases, it is safe to post-process the solution using
# `y_int = round(Int, value(y))`. However, in some cases "fixing" the
# integrality like this can cause violations in primal feasibility that exceed
# the primal feasibility tolerance. For example, if we rounded our
# `(x, y) = (1, 1e-6)` solution to `(x, y) = (1, 0)`, then the constraint
# `x <= M * y` is now violated by a value of `1.0`, which is much greater than a
# typical feasibility tolerance of `1e-8`.

primal_feasibility_report(model, Dict(x => 1.0, y => 0.0))

# ### Why you shouldn't use a small tolerance

# Just like primal feasibility tolerances, using a smaller value for the
# integrality tolerance can lead to greatly increased solve times.

# ## Contradictory results

# The distinction between feasible and infeasible can be surprisingly nuanced.
# Solver A might decide the problem is feasible while solver B might decide it
# is infeasible. Different algorithms _within_ solver A (like simplex and
# barrier) may also come to different conclusions. Even changing settings like
# turning presolve on and off can make a difference.

# Here is an example where HiGHS reports the problem is infeasible, but there
# exists a feasible (to tolerance) solution:

model = Model(HiGHS.Optimizer)
set_silent(model)
@variable(model, x >= 0)
@variable(model, y >= 0)
@constraint(model, x + 1e8 * y == -1)
optimize!(model)
@assert !is_solved_and_feasible(model)  #src
is_solved_and_feasible(model)

# The feasible solution `(x, y) = (0.0, -1e-8)` has a maximum primal violation
# of `1e-8` which is the HiGHS feasibility tolerance:

primal_feasibility_report(model, Dict(x => 0.0, y => -1e-8))

# This happens because there are two basic solutions. The first is infeasible at
# `(x, y) = (-1, 0)` and the second is feasible `(x, y) = (0, -1e-8)`. Different
# algorithms may terminate at either of these bases.

# Another example is a variation on our integrality example, but this time, there
# are constraints that `x >= 1` and `y <= 0.5`:

M = 1e6
model = Model(HiGHS.Optimizer)
set_silent(model)
@variable(model, x >= 1)
@variable(model, y, Bin)
@constraint(model, y <= 0.5)
@constraint(model, x <= M * y)
optimize!(model)
@assert !is_solved_and_feasible(model)  #src
is_solved_and_feasible(model)

# HiGHS reports the problem is infeasible, but there is a feasible (to
# tolerance) solution of:

primal_feasibility_report(model, Dict(x => 1.0, y => 1e-6))

# This happens because the presolve routine deduces that the `y <= 0.5`
# constraint forces the binary variable `y` to take the value `0`. Substituting
# the value for `y` into the last constraint, presolve may also deduce that
# `x <= 0`, which violates the bound of `x >= 1` and so the problem is
# infeasible.

# We can work around this by providing HiGHS with the feasible starting
# solution:

set_start_value(x, 1)
set_start_value(y, 1e-6)

# Now HiGHS will report that the problem is feasible:

optimize!(model)
@assert is_solved_and_feasible(model)  #src
is_solved_and_feasible(model)

# ### Contradictory results are not a bug in the solver

# These contradictory examples are not bugs in the HiGHS solver. They are an
# expected result of the interaction between the tolerances and the solution
# algorithm. There will always be models in the gray boundary at the edge of
# feasibility, for which the question of feasibility is not a clear true or
# false.

# ## Problem scaling

# Problem scaling refers to the absolute magnitudes of the data in your problem.
# The data is any numbers in the objective, the constraints, or the variable
# bounds.
#
# We say that a problem is poorly scaled if there are very small ($< 10^{-3}$)
# or very large ($> 10^6$) coefficients in the problem, or if the ratio of the
# largest to smallest coefficient is large.

# Numerical issues related to the feasibility tolerances most commonly arise
# because of poor problem scaling. The next examples assume a primal feasibility
# tolerance of 1e-8, but actual tolerances may vary from one solver to another.

# ### Small magnitudes

# If the problem data is too small, then the feasibility tolerance can be too
# easily satisfied. For example, consider:

model = Model()
@variable(model, x)
@constraint(model, 1e-8 * x == 1e-4)

# This should have the solution that $x = 10^4$, but because the feasibility
# tolerance of this constraint is $|10^{-4} - 10^{-8} * x| < 10^{-8}$, it
# actually permits any value of `x` between 9999 and 10,001, which is a larger
# range of feasible values than you might have expected.

# ### Large magnitudes

# If the problem data is too large, then the feasibility tolerance can be too
# difficult to satisfy.

model = Model()
@variable(model, x)
@constraint(model, 1e12 * x == 1e4)

# This should have the solution that $x = 10^{-8}$, but because the feasibility
# tolerance of this constraint is $|10^{12}x - 10^{4}| < 10^{-8}$, it
# actually permits any value of `x` in $10^{-8} \pm 10^{-20}$, which is a
# smaller range of feasible values than you might have expected.

# ### Large magnitude ratios

# If the ratio of the smallest to the largest magnitude is too large, then the
# tolerances or small changes in the input data can lead to large changes in the
# optimal solution. We have already seen an example with the integrality
# tolerance, but we can exacerbate the behavior by putting a small coefficient
# on `x`:

model = Model()
@variable(model, x >= 0)
@variable(model, y, Bin)
@constraint(model, 1e-6x <= 1e6 * y)

# This problem has a feasible (to tolerance) solution of:

primal_feasibility_report(model, Dict(x => 1_000_000.01, y => 1e-6))

# If you intended the constraint to read that if `x` is non-zero then `y = 1`,
# this solution might be unexpected.

# ### Recommended values

# There are no hard rules that you must follow, and the interaction between
# tolerances, problem scaling, and the solution is problem dependent. You should
# always check the solution returned by the solver to check it makes sense for
# your application.

# With that caveat in mind, a general rule of thumb to follow is:

# Try to keep the ratio of the smallest to largest coefficient less than $10^6$
# in any row and column, and try to keep most values between $10^{-3}$ and
# $10^6$.

# ### Choosing the correct units

# The best way to fix problem scaling is by changing the units of your variables
# and constraints. Here's an example. Suppose we are choosing the level of
# capacity investment in a new power plant. We can install up to 1 GW of
# capacity at a cost of $1.78/W, and we have a budget of $200 million.

model = Model()
@variable(model, 0 <= x_capacity_W <= 10^9)
@constraint(model, 1.78 * x_capacity_W <= 200e6)

# This constraint violates the recommendations because there are values greater
# than $10^6$, and the ratio of the coefficients in the constraint is $10^8$.

# One fix is the convert our capacity variable from Watts to Megawatts. This
# yields:

model = Model()
@variable(model, 0 <= x_capacity_MW <= 10^3)
@constraint(model, 1.78e6 * x_capacity_MW <= 200e6)

# We can improve our model further by dividing the constraint by $10^6$ to
# change the units from dollars to million dollars.

model = Model()
@variable(model, 0 <= x_capacity_MW <= 10^3)
@constraint(model, 1.78 * x_capacity_MW <= 200)

# This problem is equivalent to the original problem, but it has much better
# problem scaling.

# As a general rule, to fix problem scaling you must simultaneously scale both
# variables and constraints. It is usually not sufficient to scale variables or
# constraints in isolation.

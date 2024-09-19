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

# # Tolerances

# Optimization solvers can seem like magic black boxes that take in an algebraic
# formulation of a problem and return a solution. It is tempting to treat their
# solutions at face value, since we often have little ability to verify that the
# solution is in fact optimal. However, like all numerical algorithms that use
# floating point arithmetic, optimization solvers use tolerances to check
# whether a solution satisfies the constraints. In the best case, the solution
# satisfies the original constraints to machine precision. In most cases, the
# solution satisfies the constraints to some very small tolerance that has no
# noticeable impact on the quality of the optimal solution. In the worst case,
# the solver can return a "wrong" solution, or fail to find one even if it
# exists. (In the last case, the solution is wrong only in the sense of user
# expectation. It will satisfy the solution to the tolerances that are
# provided.)

# The purpose of this tutorial is to explain the various types of tolerances
# that are used in optimization solvers and what you can reasonably expect from
# a solution.

# !!! tip
#     This tutorial is more advanced than the other "Getting started" tutorials.
#     It's in the "Getting started" section to give you an early preview of how
#     tolerances affect JuMP models. However, if you are new to JuMP, you may want to
#     briefly skim the tutorial, and come back to it once you have written a few
#     JuMP models.

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
# integrality tolerance and lead to greatly increased solve times.

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
# algorthims may terminate at either of these bases.

# Another example is a variation on our integrality eample, but this time, there
# is are constraint that `x >= 1` and `y <= 0.5`:

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

# This happens becauuse the presolve routine deduces that the `y <= 0.5`
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

# ### Contradictory are not a bug in the solver

# These contradictory examples are not bugs in the HiGHS solver. They are an
# expected result of the interaction between the tolerances and the solution
# algorithm. There will always be models in the gray boundary at the edge of
# feasibility, for which the question of feasibility is not a clear true or
# false.

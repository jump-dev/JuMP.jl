# Copyright 2017, Iain Dunning, Joey Huchette, Miles Lubin, and contributors    #src
# This Source Code Form is subject to the terms of the Mozilla Public License   #src
# v.2.0. If a copy of the MPL was not distributed with this file, You can       #src
# obtain one at https://mozilla.org/MPL/2.0/.                                   #src

# # Basis matrices

# This tutorial explains how to query basis of a linear program.

# ## Setup

# This tutorial uses the following packages:

using JuMP
import HiGHS

# ## Standard form example

# Consider the following example, which is from the Wikipedia article on
# [Basic feasible solutions](https://en.wikipedia.org/wiki/Basic_feasible_solution):
#
# ```math
# \begin{aligned}
# \max \;     & 0 \\
# \text{s.t.} & x_1 + 5x_2 + 3x_3 + 4x_4 + 6x_5 = 14 \\
#             & x_2 + 3x_3 + 5x_4 + 6x_5 = 7 \\
#             & x_i \ge 0, \forall i = 1,\ldots,5.
# \end{aligned}
# ```
#
# The `A` matrix is:

A = [1 5 3 4 6; 0 1 3 5 6]

# and the right-hand side `b` vector is:

b = [14, 7]

# We can create an optimize the problem in the standard form:

n = size(A, 2)
model = Model(HiGHS.Optimizer)
@variable(model, x[1:n] >= 0)
@constraint(model, A * x == b)
optimize!(model)
@assert is_solved_and_feasible(model)

# This has a solution:

value.(x)

# Query the basis status using [`MOI.VariableBasisStatus`](@ref):

MOI.get(model, MOI.VariableBasisStatus(), x[1])

# the result is a [`MOI.BasisStatusCode`](@ref). Query all of the basis statuses
# with the broadcast `MOI.get.(`:

MOI.get.(model, MOI.VariableBasisStatus(), x)

# The values are either [`MOI.NONBASIC_AT_LOWER`](@ref) or [`MOI.BASIC`](@ref).
# All of the [`MOI.NONBASIC_AT_LOWER`](@ref) have a value at their lower bound.
# The [`MOI.BASIC`](@ref) variables correspond to the columns of the optimal
# basis.

# We can get the columns using:

indices =
    [MOI.get(model, MOI.VariableBasisStatus(), x[i]) == MOI.BASIC for i in 1:5]

# Filter the basis matrix from `A`:

B = A[:, indices]

# Since the basis matrix is non-singular, solving the system `Bx = b` must yield
# the optimal primal solution:

B \ b

#-

value.(x[indices])

# ## A more complicated example

# Often, you may want to work with the basis of a model that is not in a nice
# standard form. For example:

model = Model(HiGHS.Optimizer)
@variable(model, x >= 0)
@variable(model, 0 <= y <= 3)
@variable(model, z <= 1)
@objective(model, Min, 12x + 20y - z)
@constraint(model, c1, 6x + 8y >= 100)
@constraint(model, c2, 7x + 12y >= 120)
@constraint(model, c3, x + y <= 20)
optimize!(model)
@assert is_solved_and_feasible(model)
solution_summary(model; verbose = true)

# A common way to query the basis status of every variable is:

v_basis = Dict(
    xi => MOI.get(model, MOI.VariableBasisStatus(), xi) for
    xi in all_variables(model)
)

# Despite having three constraints, there are only two basis variables. Since
# the basis matrix must be square, where is the other basis variable?

# The answer is that solvers will reformulate inequality constraints:
# ```math
# \begin{aligned}
# \text{s.t.} & l \le A x \le u
# \end{aligned}
# ```
# into the system:
# ```math
# \begin{aligned}
# \text{s.t.} & A x - Is = 0
#             & l \le s \le u
# \end{aligned}
# ```
# Thus, for every inequality constraint there is a slack variable `s`.

# Query the basis status of the slack variables associated with a constraint
# using [`MOI.ConstraintBasisStatus`](@ref):

c_basis = Dict(
    ci => MOI.get(model, MOI.ConstraintBasisStatus(), ci) for ci in
    all_constraints(model; include_variable_in_set_constraints = false)
)

# Thus, the basis is formed by `x`, `y`, and the slack associated with `c3`.

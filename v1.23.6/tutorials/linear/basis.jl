# Copyright 2017, Iain Dunning, Joey Huchette, Miles Lubin, and contributors    #src
# This Source Code Form is subject to the terms of the Mozilla Public License   #src
# v.2.0. If a copy of the MPL was not distributed with this file, You can       #src
# obtain one at https://mozilla.org/MPL/2.0/.                                   #src

# # Basis matrices

# This tutorial explains how to query the basis of a linear program.

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
# \text{s.t.}\; & 1x_1 + 5x_2 + 3x_3 + 4x_4 + 6x_5 = 14 \\
#             & 0x_1 + 1x_2 + 3x_3 + 5x_4 + 6x_5 = 7 \\
#             & x_i \ge 0, \forall i = 1,\ldots,5.
# \end{aligned}
# ```
#
# The `A` matrix is:

A = [1 5 3 4 6; 0 1 3 5 6]

# and the right-hand side `b` vector is:

b = [14, 7]

# We can create and optimize the problem in the standard form:

n = size(A, 2)
model = Model(HiGHS.Optimizer)
set_silent(model)
@variable(model, x[1:n] >= 0)
@constraint(model, A * x == b)
optimize!(model)
@assert is_solved_and_feasible(model)

# This has a solution:

value.(x)

# Query the basis status of a variable using [`MOI.VariableBasisStatus`](@ref):

get_attribute(x[1], MOI.VariableBasisStatus())

# the result is a [`MOI.BasisStatusCode`](@ref). Query all of the basis statuses
# with the broadcast `get_attribute.(`:

get_attribute.(x, MOI.VariableBasisStatus())

# For this problem, the values are either [`MOI.NONBASIC_AT_LOWER`](@ref) or [`MOI.BASIC`](@ref).
# All of the [`MOI.NONBASIC_AT_LOWER`](@ref) variables have a value at their
# lower bound. The [`MOI.BASIC`](@ref) variables correspond to the columns of
# the optimal basis.

# Get the columns using:

indices = get_attribute.(x, MOI.VariableBasisStatus()) .== MOI.BASIC

# Filter the basis matrix from `A`:

B = A[:, indices]

# Since the basis matrix is non-singular, solving the system `Bx = b` must yield
# the optimal primal solution of the basic variables:

B \ b

#-

value.(x[indices])

# ## A more complicated example

# Often, you may want to work with the basis of a model that is not in a nice
# standard form. For example:

model = Model(HiGHS.Optimizer)
set_silent(model)
@variable(model, x >= 0)
@variable(model, 0 <= y <= 3)
@variable(model, z <= 1)
@objective(model, Min, 12x + 20y - z)
@constraint(model, c1, 6x + 8y >= 100)
@constraint(model, c2, 7x + 12y >= 120)
@constraint(model, c3, x + y <= 20)
optimize!(model)
@assert is_solved_and_feasible(model)

# A common way to query the basis status of every variable is:

v_basis = Dict(
    xi => get_attribute(xi, MOI.VariableBasisStatus()) for
    xi in all_variables(model)
)

# Despite the model having three constraints, there are only two basic
# variables. Since the basis matrix must be square, where is the other basic
# variable?

# The answer is that solvers will reformulate inequality constraints:
# ```math
# A x \le b
# ```
# into the system:
# ```math
# A x + Is = b
# ```
# Thus, for every inequality constraint there is a slack variable `s`.

# Query the basis status of the slack variables associated with a constraint
# using [`MOI.ConstraintBasisStatus`](@ref):

c_basis = Dict(
    ci => get_attribute(ci, MOI.ConstraintBasisStatus()) for ci in
    all_constraints(model; include_variable_in_set_constraints = false)
)

# Thus, the basis is formed by `x`, `y`, and the slack associated with `c3`.

# A simple way to get the `A` matrix of an unstructured linear program is with
# [`lp_matrix_data`](@ref):

matrix = lp_matrix_data(model)
matrix.A

# You can check the permutation of the rows and columns using

matrix.variables

# and

matrix.affine_constraints

# We can construct the slack column associated with `c3` as:

s_column = zeros(size(matrix.A, 1))
s_column[3] = 1.0

# The full basis matrix is therefore:

B = hcat(matrix.A[:, [1, 2]], s_column)

# [`lp_matrix_data`](@ref) returns separate vectors for the lower and upper row
# bounds. Convert to a single right-hand side vector by taking the finite
# elements:

b = ifelse.(isfinite.(matrix.b_lower), matrix.b_lower, matrix.b_upper)

# Solving the Basis system as before yields:

B \ b

# which is the value of `x`, `y`, and the slack associated with `c3`.

# ## Identifying degenerate variables

# Another common task is identifying degenerate variables. A degenerate variable
# is a basic variable that has an optimal value at its lower or upper bound.

# Here is a function that computes whether a variable is degenerate:

function is_degenerate(x)
    if get_attribute(x, MOI.VariableBasisStatus()) == MOI.BASIC
        return (has_lower_bound(x) && ≈(value(x), lower_bound(x))) ||
               (has_upper_bound(x) && ≈(value(x), upper_bound(x)))
    end
    return false
end

# A simple example of a linear program with a degenerate solution is:

A, b, c = [1 1; 0 1], [1, 1], [1, 1]
model = Model(HiGHS.Optimizer);
set_silent(model)
@variable(model, x[1:2] >= 0)
@objective(model, Min, c' * x)
@constraint(model, A * x == b)
optimize!(model)
degenerate_variables = filter(is_degenerate, all_variables(model))

# The solution is degenerate because:

value(x[1])

# and

get_attribute(x[1], MOI.VariableBasisStatus())

# Copyright 2017, Iain Dunning, Joey Huchette, Miles Lubin, and contributors    #src
# This Source Code Form is subject to the terms of the Mozilla Public License   #src
# v.2.0. If a copy of the MPL was not distributed with this file, You can       #src
# obtain one at https://mozilla.org/MPL/2.0/.                                   #src

# # The `VectorNonlinearOracle` set

# Many nonlinear solvers support constraints of the form:
# ```math
# l \le f(x) \le u
# ```
# where $l$, $u$, and $f(x)$ are vectors. The purpose of this tutorial is to
# explain how to add this constraint using the [`MOI.VectorNonlinearOracle`](@ref)
# set.

# ## Required packages

# This tutorial uses the following packages:

using JuMP
import Ipopt
import Test

# ## Example

# As an example, we use the model from [Quadratically constrained programs](@ref).
# It can be formulated as:

model = Model()
@variable(model, x[1:3])
@objective(model, Max, x[1])
@constraint(model, sum(x) == 1)
@constraint(model, x[1]^2 + x[2]^2 - x[3]^2 <= 0)
@constraint(model, x[1]^2 - x[2] * x[3] <= 0)
print(model)

# In math, the constraints are of the form:
# ```math
# \begin{bmatrix}1\\-\infty\\-\infty\end{bmatrix}
# \le
# \begin{bmatrix}
#     x_1 + x_2 + x_3       \\
#     x_1^2 + x_2^2 - x_3^2 \\
#     x_1^2 - x_2 x_3
# \end{bmatrix}
# \le
# \begin{bmatrix} 1\\0\\0 \end{bmatrix}
# ```

# ## Inputs to the `VectorNonlinearOracle` set

# To create a [`MOI.VectorNonlinearOracle`](@ref) set, we need a few components.

# First, we need a function with the signature `(ret::AbstractVector, x::AbstractVector)`
# that evaluates $f(x)$ and stores the result in `ret`:

function eval_f(ret::AbstractVector, x::AbstractVector)
    ret[1] = sum(x)
    ret[2] = x[1]^2 + x[2]^2 - x[3]^2
    ret[3] = x[1]^2 - x[2] * x[3]
    return
end

# Here's an example of it in action:

x = [1.0, 2.0, 3.0]
f_ret = zeros(3)
eval_f(f_ret, x)
f_ret

# Next, we need the sparse Jacobian of $f$. In math, the Jacobian of $f$ is
# ```math
# \nabla f(x) = \begin{bmatrix}
#     1 & 1 & 1 \\
#     2x_1 & 2x_2 & -2x_3 \\
#     2x_1 & - x_3 & -x_2
# \end{bmatrix}
# ```
# In code, the Jacobian has two components: a vector of the `(row, column)`
# tuples indicating the structural non-zeros of the Jacobian, and a function
# with the signature `(ret::AbstractVector, x::AbstractVector)` that evaluates
# the non-zeros and stores the result in `ret`. Note that our example Jacobian
# is dense. If it was sparse, we could omit some of the structure tuples to end
# up with a vector with fewer than 9 elements.

jacobian_structure = [
    ## The first constraint
    (1, 1),
    (1, 2),
    (1, 3),
    ## The second constraint
    (2, 1),
    (2, 2),
    (2, 3),
    ## The third constraint
    (3, 1),
    (3, 2),
    (3, 3),
]

function eval_jacobian(ret::AbstractVector, x::AbstractVector)
    ## The first constraint
    ret[1] = 1.0
    ret[2] = 1.0
    ret[3] = 1.0
    ## The second constraint
    ret[4] = 2.0 * x[1]
    ret[5] = 2.0 * x[2]
    ret[6] = -2.0 * x[3]
    ## The third constraint
    ret[7] = 2.0 * x[1]
    ret[8] = -1.0 * x[3]
    ret[9] = -1.0 * x[2]
    return
end

# Here's an example of it in action:

J_ret = zeros(9)
eval_jacobian(J_ret, x)
J_ret

# Now we're ready to create the [`MOI.VectorNonlinearOracle`](@ref) set. In
# addition to `eval_f`, `jacobian_structure`, and `eval_jacobian`, we need to
# provide the input dimension of `x` (in this case, `3`), and the $l$ and $u$
# vectors:

set = MOI.VectorNonlinearOracle(;
    dimension = 3,
    l = [1.0, -Inf, -Inf],
    u = [1.0, 0.0, 0.0],
    eval_f,
    jacobian_structure,
    eval_jacobian,
)

# ## Using the set in a JuMP model

# Now we can create a JuMP model using the `x in set` syntax:

model = Model()
@variable(model, x[1:3])
@objective(model, Max, x[1])
@constraint(model, x in set)
print(model)

# Now we can solve and check the solution:

set_optimizer(model, Ipopt.Optimizer)
set_silent(model)
optimize!(model)
assert_is_solved_and_feasible(model)
Test.@test objective_value(model) ≈ 0.32699 atol = 1e-5
value(x)

# ## Hessians

# To improve performance, we can also compute and pass the explicit Hessian of
# the Lagrangian
# ```math
# L(x, u) = \nabla^2 u^\top f(x)
# ```
# In math, the Hessian of our example is:
# ```math
# \nabla^2 L(x, u) = \begin{bmatrix}
#     2(u_2 + u_3) & \cdot & \cdot \\
#     \cdot & 2u_2 & -u_3 \\
#     \cdot & -u_3 & -2u_2
# \end{bmatrix}
# ```
# (If you need the practice, try deriving it.)

# Like the Jacobian, we need the sparsity structure and a function to compute
# the non-zeros. Importantly, because the Hessian is symmetric, we need to pass
# only the upper triangular values.

hessian_lagrangian_structure = [(1, 1), (2, 2), (2, 3), (3, 3)]

function eval_hessian_lagrangian(
    ret::AbstractVector,
    x::AbstractVector,
    u::AbstractVector,
)
    ret[1] = 2.0 * (u[2] + u[3])
    ret[2] = 2.0 * u[2]
    ret[3] = -1.0 * u[3]
    ret[4] = -2.0 * u[2]
    return
end

# Here's an example of it in action:

H_ret = zeros(4)
u = [1.0, 2.0, 3.0]
eval_hessian_lagrangian(H_ret, x, u)
H_ret

# Putting it all together:

model = Model(Ipopt.Optimizer)
set_silent(model)
@variable(model, x[1:3])
@objective(model, Max, x[1])
set = MOI.VectorNonlinearOracle(;
    dimension = 3,
    l = [1.0, -Inf, -Inf],
    u = [1.0, 0.0, 0.0],
    eval_f,
    jacobian_structure,
    eval_jacobian,
    hessian_lagrangian_structure,
    eval_hessian_lagrangian,
)
@constraint(model, x in set)
optimize!(model)
assert_is_solved_and_feasible(model)
Test.@test objective_value(model) ≈ 0.32699 atol = 1e-5
value(x)

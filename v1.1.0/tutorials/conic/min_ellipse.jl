# Copyright 2017, Iain Dunning, Joey Huchette, Miles Lubin, and contributors    #src
# This Source Code Form is subject to the terms of the Mozilla Public License   #src
# v.2.0. If a copy of the MPL was not distributed with this file, You can       #src
# obtain one at https://mozilla.org/MPL/2.0/.                                   #src

# # Minimal ellipses

# This example comes from section 8.4.1 of the book *Convex Optimization* by [Boyd
# and Vandenberghe (2004)](https://web.stanford.edu/~boyd/cvxbook/).

# Given a set of ``m`` ellipses of the form 
# ```math
# E(A, b, c) = \{ x : x' A x + 2 b' x + c \leq 0 \},
# ```
# we find the ellipse of smallest area that encloses the given ellipses.

# It is convenient to parameterize the minimal enclosing ellipse as
# ```math
# \{ x : || Px + q || \leq 1 \}.
# ```

# Then the optimal ``P`` and ``q`` are given by the convex semidefinite program
# ```math
# \begin{aligned}
# \text{maximize } & \quad \log(\det(P))  \\
# \text{subject to } & \quad \tau_i \geq 0, & i = 1, \ldots, m, \\
# & \quad\begin{bmatrix}
#     P^2 - \tau_i A_i      &  P q - \tau_i b_i   &    0     \\
#     (P q - \tau_i b_i)^T  &   -1 - \tau_i c_i   &  (P q)^T \\
#                 0  &          (P q)   & - P^2    \\
# \end{bmatrix}
# \preceq 0 \text{ (PSD) } &  i=1, \ldots, m
# \end{aligned}
# ```
# with helper variables ``\tau``. 
# 
# The program can be solved by using a variable representing ``P^2`` (`Psqr` in the Julia code),
# a vector of variables ``\tilde{q}`` (`q_tilde`) in place of ``P q`` and the variables ``\tau`` (`tau[i]`).

# This tutorial uses the following packages:

using JuMP
using SCS
using Plots
using Test

# ## Set-up
# First, define the ``m`` input ellipses (here ``m = 6``),
# parameterized as ``x^T A_i x + 2 b_i^T x + c \leq 0``:
As = [
    [1.2576 -0.3873; -0.3873 0.3467],
    [1.4125 -2.1777; -2.1777 6.7775],
    [1.7018 0.8141; 0.8141 1.7538],
    [0.9742 -0.7202; -0.7202 1.5444],
    [0.6798 -0.1424; -0.1424 0.6871],
    [0.1796 -0.1423; -0.1423 2.6181],
];
bs = [
    [0.2722, 0.1969],
    [-1.228, -0.0521],
    [-0.4049, 1.5713],
    [0.0265, 0.5623],
    [-0.4301, -1.0157],
    [-0.3286, 0.557],
];
cs = [0.1831, 0.3295, 0.2077, 0.2362, 0.3284, 0.4931];

# We visualise the ellipses using the Plots package:
pl = plot(size = (600, 600))
thetas = range(0, 2pi + 0.05, step = 0.05)
for (A, b, c) in zip(As, bs, cs)
    sqrtA = sqrt(A)
    b_tilde = sqrtA \ b
    alpha = b' * (A \ b) - c
    rhs = hcat(
        sqrt(alpha) * cos.(thetas) .- b_tilde[1],
        sqrt(alpha) * sin.(thetas) .- b_tilde[2],
    )
    ellipse = sqrtA \ rhs'
    plot!(pl, ellipse[1, :], ellipse[2, :], label = nothing, c = :navy)
end
plot(pl)

# ## Build the model
# Now let's build the initial model, using the change-of-variables
# `Psqr` = ``P^2`` and `q_tilde` = ``P q``:
model = Model(SCS.Optimizer)
m = length(As)
n, _ = size(first(As))
@variable(model, tau[1:m] ≥ 0)
@variable(model, Psqr[1:n, 1:n], PSD)
@variable(model, q_tilde[1:n])
@variable(model, logdetP);

# Next, create the PSD constraints and objective:
for (A, b, c, t) in zip(As, bs, cs, tau)
    if !(isreal(A) && transpose(A) == A)
        @error "Input matrices need to be real, symmetric matrices."
    end
    @constraint(
        model,
        -[
            #! format: off
            Psqr-t*A           q_tilde-t*b zeros(n, n)
            (q_tilde - t * b)' -1-t*c      q_tilde'
            zeros(n, n)        q_tilde     -Psqr
            #! format: on
        ] in PSDCone()
    )
end

@constraint(
    model,
    [logdetP; [Psqr[i, j] for i in 1:n for j in i:n]] in MOI.RootDetConeTriangle(n)
);
@objective(model, Max, logdetP);

# Note that here the root-determinant cone is used for constructing the objective function.
# While the more consistent choice for the mathematical formulation is to use
# `MOI.LogDetConeTriangle(n)` instead, `MOI.RootDetConeTriangle(n)` will produce equivalent
# optimal solutions and is found to be more efficient for the SCS solver for this example.

# Now, solve the program:
optimize!(model)
@test termination_status(model) == OPTIMAL;
@test primal_status(model) == FEASIBLE_POINT;

# ## Results
# After solving the model to optimality
# we can recover the original solution parameterization as
P = sqrt(value.(Psqr))
q = P \ value.(q_tilde)

# We can test that we get the expected results to within approximation tolerance.
@test isapprox(P, [0.4237 -0.0396; -0.0396 0.3163], atol = 1e-2);
@test isapprox(q, [-0.3960, -0.0214], atol = 1e-2);

# Finally, overlaying the solution in the plot we see the minimal area enclosing ellipsoid.
plot!(
    pl,
    [tuple(P \ ([cos(theta), sin(theta)] - q)...) for theta in thetas],
    c = :crimson,
    label = nothing,
)
plot(pl)

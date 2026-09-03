# Copyright 2017, Iain Dunning, Joey Huchette, Miles Lubin, and contributors    #src
# This Source Code Form is subject to the terms of the Mozilla Public License   #src
# v.2.0. If a copy of the MPL was not distributed with this file, You can       #src
# obtain one at https://mozilla.org/MPL/2.0/.                                   #src

# # Minimal ellipses

# This example comes from section 8.4.1 of the book *Convex Optimization* by
# [Boyd and Vandenberghe (2004)](https://web.stanford.edu/~boyd/cvxbook/).

# ## Formulation

# Given a set of ``m`` ellipses of the form:
# ```math
# E(A, b, c) = \{ x : x^\top A x + 2 b^\top x + c \leq 0 \},
# ```
# the minimal ellipse problem finds an ellipse with the minimum area that
# encloses the given ellipses.

# It is convenient to parameterize the minimal enclosing ellipse as
# ```math
# \{ x : || Px + q || \leq 1 \}.
# ```

# Then the optimal ``P`` and ``q`` are given by the convex semidefinite program;
# ```math
# \begin{aligned}
# \text{maximize } & \quad \log(\det(P))  \\
# \text{subject to } & \quad \tau_i \geq 0, & i = 1, \ldots, m \\
# & \quad\begin{bmatrix}
#     P^2 - \tau_i A_i        &  P q - \tau_i b_i   &    0        \\
#     (P q - \tau_i b_i)^\top &   -1 - \tau_i c_i   &  (P q)^\top \\
#                 0           &          (P q)      & - P^2       \\
# \end{bmatrix}
# \preceq 0 \text{ (PSD) } &  i=1, \ldots, m
# \end{aligned}
# ```
# with helper variables ``\tau``.

# ## Required packages

# This tutorial uses the following packages:

using JuMP
import LinearAlgebra
import Plots
import SCS
import Test  #src

# ## Data

# First, define the ``m`` input ellipses (here ``m = 6``), parameterized as
# ``x^T A_i x + 2 b_i^T x + c \leq 0``:

struct Ellipse
    A::Matrix{Float64}
    b::Vector{Float64}
    c::Float64
    function Ellipse(A::Matrix{Float64}, b::Vector{Float64}, c::Float64)
        @assert isreal(A) && LinearAlgebra.issymmetric(A)
        return new(A, b, c)
    end
end

ellipses = [
    Ellipse([1.2576 -0.3873; -0.3873 0.3467], [0.2722, 0.1969], 0.1831),
    Ellipse([1.4125 -2.1777; -2.1777 6.7775], [-1.228, -0.0521], 0.3295),
    Ellipse([1.7018 0.8141; 0.8141 1.7538], [-0.4049, 1.5713], 0.2077),
    Ellipse([0.9742 -0.7202; -0.7202 1.5444], [0.0265, 0.5623], 0.2362),
    Ellipse([0.6798 -0.1424; -0.1424 0.6871], [-0.4301, -1.0157], 0.3284),
    Ellipse([0.1796 -0.1423; -0.1423 2.6181], [-0.3286, 0.557], 0.4931),
];

# We visualise the ellipses using the Plots package:

function plot_ellipse(plot, ellipse::Ellipse)
    A, b, c = ellipse.A, ellipse.b, ellipse.c
    θ = range(0, 2pi + 0.05; step = 0.05)
    ## Some linear algebra to convert θ into (x,y) coordinates.
    x_y = √A \ (√(b' * (A \ b) - c) .* hcat(cos.(θ), sin.(θ)) .- (√A \ b)')'
    Plots.plot!(plot, x_y[1, :], x_y[2, :]; label = nothing, c = :navy)
    return
end

plot = Plots.plot(; size = (600, 600))
for ellipse in ellipses
    plot_ellipse(plot, ellipse)
end
plot

# ## Build the model

# Now let's build the model, using the change-of-variables `P²` = ``P^2`` and
# `P_q` = ``P q``. We'll recover the true value of `P` and `q` after the solve.

model = Model(SCS.Optimizer)
## We need to use a tighter tolerance for this example, otherwise the bounding
## ellipse won't actually be bounding...
set_optimizer_attribute(model, "eps_rel", 1e-6)
set_silent(model)
m, n = length(ellipses), size(first(ellipses).A, 1)
@variable(model, τ[1:m] >= 0)
@variable(model, P²[1:n, 1:n], PSD)
@variable(model, P_q[1:n])

for (i, ellipse) in enumerate(ellipses)
    A, b, c = ellipse.A, ellipse.b, ellipse.c
    X = [
        #! format: off
        (P² - τ[i] * A)   (P_q - τ[i] * b) zeros(n, n)
        (P_q - τ[i] * b)' (-1 - τ[i] * c)  P_q'
        zeros(n, n)       P_q              -P²
        #! format: on
    ]
    @constraint(model, LinearAlgebra.Symmetric(X) <= 0, PSDCone())
end

# We cannot directly represent the objective ``\log(\det(P))``, so we introduce
# the conic reformulation:

@variable(model, log_det_P)
@constraint(model, [log_det_P; 1; vec(P²)] in MOI.LogDetConeSquare(n))
@objective(model, Max, log_det_P)

# Now, solve the program:

optimize!(model)
Test.@test termination_status(model) == OPTIMAL    #src
Test.@test primal_status(model) == FEASIBLE_POINT  #src
solution_summary(model)

# ## Results

# After solving the model to optimality we can recover the solution in terms of
# ``P`` and ``q``:

P = sqrt(value.(P²))
q = P \ value.(P_q)

# Finally, overlaying the solution in the plot we see the minimal area enclosing
# ellipsoid:

Test.@test isapprox(P, [0.4237 -0.0396; -0.0396 0.3163]; atol = 1e-2)  #src
Test.@test isapprox(q, [-0.3960, -0.0214]; atol = 1e-2)                #src

Plots.plot!(
    plot,
    [tuple(P \ [cos(θ) - q[1], sin(θ) - q[2]]...) for θ in 0:0.05:(2pi+0.05)];
    c = :crimson,
    label = nothing,
)

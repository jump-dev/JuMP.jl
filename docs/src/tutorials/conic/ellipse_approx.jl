# Copyright 2017, Iain Dunning, Joey Huchette, Miles Lubin, and contributors    #src
# This Source Code Form is subject to the terms of the Mozilla Public License   #src
# v.2.0. If a copy of the MPL was not distributed with this file, You can       #src
# obtain one at https://mozilla.org/MPL/2.0/.                                   #src

# # Ellipsoid approximation

# This example considers the problem of computing _extremal ellipsoids_:
# finding ellipsoids that best approximate a given set.
# Our example will focus on outer approximations of finite sets
# of points.

# This example comes from section 4.9 "Applications VII: Extremal Ellipsoids"
# of the book *Lectures on Modern Convex Optimization* by
# [Ben-Tal and Nemirovski (2001)](http://epubs.siam.org/doi/book/10.1137/1.9780898718829).

# ## Problem formulation

# Suppose that we are given set ``S`` of ``m`` points in ``n``-dimensional space:
# ```math
# \mathcal{S} = \{ x_1, \ldots, x_m \} \subset \mathbb{R}^n
# ```
# Our goal is to determine an optimal vector ``c \in  \mathbb{R}^n`` and
# an optimal ``n \times n`` real symmetric matrix ``D`` such that the ellipse
# ```math
# E(D, c) = \{ x : (x - c)^\top D ( x - c) \leq 1 \},
# ```
# contains ``\mathcal{S}`` and such that this ellipse has
# the smallest possible volume.

# The optimal ``D`` and ``c`` are given by the convex semidefinite program
# ```math
# \begin{aligned}
# \text{maximize }   && \quad (\det(Z))^{\frac{1}{n}}  & \\
# \text{subject to } && \quad Z \; & \succeq \; 0 & \text{ (PSD) }, & \\
# && \quad\begin{bmatrix}
#     s  &  z^\top   \\
#     z  &  Z        \\
# \end{bmatrix}
#  \; & \succeq \; 0 & \text{ (PSD) }, & \\
# && x_i^\top Z x_i - 2x_i^\top z + s \; & \leq \; 0 &  i=1, \ldots, m &
# \end{aligned}
# ```
# with matrix variable ``Z``, vector variable ``z`` and real variables ``t, s``.
# The optimal solution ``(t_*, Z_*, z_*, s_*)`` gives the solution
# ```math
# D = Z_*, \quad c = Z_*^{-1} z_*.
# ```

# ## Required packages

# This tutorial uses the following packages:

using JuMP
import LinearAlgebra
import Random
import Plots
import SCS
import Test  #src

# ## Data

# We need some ``m`` points 

function generate_point_cloud(
    m;       # number of 2-dimensional points
    a = 10,  # scaling in x direction
    b = 2,   # scaling in y direction
    rho = deg2rad(30),  # rotation of points around origin
    random_seed = 1,
)
    rng = Random.MersenneTwister(random_seed)
    P = randn(rng, Float64, m, 2)
    Phi = [a*cos(rho) a*sin(rho); -b*sin(rho) b*cos(rho)]
    S = P * Phi
    return S
end

# For the sake of this example, let's take ``m = 600``:
S = generate_point_cloud(600)

# We will visualise the points using the Plots package:

function plot_point_cloud(plot, S; r = 1.1 * maximum(abs.(S)), colour = :green)
    Plots.scatter!(
        plot,
        S[:, 1],
        S[:, 2];
        xlim = (-r, r),
        ylim = (-r, r),
        label = nothing,
        c = colour,
        shape = :x,
    )
    return
end

plot = Plots.plot(; size = (600, 600))
plot_point_cloud(plot, S)
plot

# ## JuMP formulation

# Now let's build the JuMP model. We'll be able to compute
# ``D`` and ``c`` after the solve.

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

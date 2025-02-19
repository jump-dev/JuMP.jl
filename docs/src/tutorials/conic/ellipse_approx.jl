# Copyright 2017, Iain Dunning, Joey Huchette, Miles Lubin, and contributors    #src
# This Source Code Form is subject to the terms of the Mozilla Public License   #src
# v.2.0. If a copy of the MPL was not distributed with this file, You can       #src
# obtain one at https://mozilla.org/MPL/2.0/.                                   #src

# # Example: ellipsoid approximation

# This tutorial considers the problem of computing _extremal ellipsoids_:
# finding ellipsoids that best approximate a given set. As an extension, we show
# how to use JuMP to inspect the bridges that were used, and how to explore
# alternative formulations.

# The model comes from Section 4.9 of [BenTal2001](@cite).

# For a related example, see also the [Example: minimal ellipses](@ref) tutorial.

# ## Required packages

# This tutorial uses the following packages:

using JuMP
import Clarabel
import LinearAlgebra
import Plots
import Random
import Test

# ## Problem formulation

# Suppose that we are given a set ``\mathcal{S}`` consisting of ``m`` points in
# ``n``-dimensional space:
# ```math
# \mathcal{S} = \{ x_1, \ldots, x_m \} \subset \mathbb{R}^n
# ```
# Our goal is to determine an optimal vector ``c \in  \mathbb{R}^n`` and
# an optimal ``n \times n`` real symmetric matrix ``D`` such that the ellipse:
# ```math
# E(D, c) = \{ x : (x - c)^\top D ( x - c) \leq 1 \},
# ```
# contains ``\mathcal{S}`` and has the smallest possible volume.

# The optimal ``D`` and ``c`` are given by the optimization problem:
# ```math
# \begin{aligned}
# \max        \quad & t                                                                    \\
# \text{s.t.} \quad & Z \succeq 0 \\
#                   & \begin{bmatrix} s & z^\top \\ z & Z \\ \end{bmatrix} \succeq 0 \\
#                   & x_i^\top Z x_i - 2x_i^\top z + s \leq 1    & \quad i=1, \ldots, m \\
#                   & t \le \sqrt[n]{\det(Z)},
# \end{aligned}
# ```
# where ``D = Z_*`` and ``c = Z_*^{-1} z_*``.

# ## Data

# We first need to generate some points to work with.

function generate_point_cloud(
    m;            # number of 2-dimensional points
    a = 10,       # scaling in x direction
    b = 2,        # scaling in y direction
    rho = π / 6,  # rotation of points around origin
    random_seed = 1,
)
    rng = Random.MersenneTwister(random_seed)
    P = randn(rng, Float64, m, 2)
    Phi = [a*cos(rho) a*sin(rho); -b*sin(rho) b*cos(rho)]
    S = P * Phi
    return S
end

# For the sake of this example, let's take ``m = 100``:
S = generate_point_cloud(100);

# We will visualise the points (and ellipse) using the Plots package:

r = 1.1 * maximum(abs.(S))
plot = Plots.scatter(
    S[:, 1],
    S[:, 2];
    xlim = (-r, r),
    ylim = (-r, r),
    label = nothing,
    c = :green,
    shape = :x,
    size = (600, 600),
)

# ## JuMP formulation

# Now let's build and the JuMP model. We'll compute ``D`` and ``c`` after the
# solve.

model = Model(Clarabel.Optimizer)
set_silent(model)
m, n = size(S)
@variable(model, z[1:n])
@variable(model, Z[1:n, 1:n], PSD)
@variable(model, s)
@variable(model, t)
@constraint(model, [s z'; z Z] >= 0, PSDCone())
@constraint(
    model,
    [i in 1:m],
    S[i, :]' * Z * S[i, :] - 2 * S[i, :]' * z + s <= 1,
)
@constraint(model, [t; vec(Z)] in MOI.RootDetConeSquare(n))
@objective(model, Max, t)
optimize!(model)
assert_is_solved_and_feasible(model)
solution_summary(model)

# ## Results

# After solving the model to optimality we can recover the solution in terms of
# ``D`` and ``c``:

D = value.(Z)

#-

c = D \ value.(z)

# We can check that each point lies inside the ellipsoid, by checking if the
# largest normalized radius is less than 1:

largest_radius = maximum(map(x -> (x - c)' * D * (x - c), eachrow(S)))

# Finally, overlaying the solution in the plot we see the minimal volume
# approximating ellipsoid:

P = sqrt(D)
q = -P * c
data = [tuple(P \ [cos(θ) - q[1], sin(θ) - q[2]]...) for θ in 0:0.05:(2pi+0.05)]
Plots.plot!(plot, data; c = :crimson, label = nothing)

# ## Alternative formulations

# The formulation of `model` uses [`MOI.RootDetConeSquare`](@ref). However,
# because Clarabel does not natively support this cone, JuMP automatically
# reformulates the problem into an equivalent problem that Clarabel _does_ support.
# You can see the reformulation that JuMP chose using [`print_active_bridges`](@ref):

print_active_bridges(model)

# There's a lot going on here, but the first bullet is:
# ```raw
# * Unsupported objective: MOI.VariableIndex
# |  bridged by:
# |   MOIB.Objective.FunctionConversionBridge{Float64}
# |  may introduce:
# |   * Supported objective: MOI.ScalarAffineFunction{Float64}
# ```
# This says that Clarabel does not support a `MOI.VariableIndex` objective
# function, and that JuMP used a [`MOI.Bridges.Objective.FunctionConversionBridge`](@ref)
# to convert it into a `MOI.ScalarAffineFunction{Float64}` objective function.

# We can leave JuMP to do the reformulation, or we can rewrite our model to
# have an objective function that Clarabel natively supports:

@objective(model, Max, 1.0 * t + 0.0);

# Re-printing the active bridges:

print_active_bridges(model)

# we get `* Supported objective: MOI.ScalarAffineFunction{Float64}`.

# We can manually implement some other reformulations to change our model to
# something that Clarabel more closely supports by:
#
#  * Replacing the [`MOI.VectorOfVariables`](@ref) in [`MOI.PositiveSemidefiniteConeTriangle`](@ref)
#    constraint `@variable(model, Z[1:n, 1:n], PSD)` with the
#    [`MOI.VectorAffineFunction`](@ref) in [`MOI.PositiveSemidefiniteConeTriangle`](@ref)
#    `@constraint(model, Z >= 0, PSDCone())`.
#
#  * Replacing the [`MOI.VectorOfVariables`](@ref) in [`MOI.PositiveSemidefiniteConeSquare`](@ref)
#    constraint `[s z'; z Z] >= 0, PSDCone()` with the
#    [`MOI.VectorAffineFunction`](@ref) in [`MOI.PositiveSemidefiniteConeTriangle`](@ref)
#    `@constraint(model, LinearAlgebra.Symmetric([s z'; z Z]) >= 0, PSDCone())`.
#
#  * Replacing the [`MOI.ScalarAffineFunction`](@ref) in [`MOI.GreaterThan`](@ref)
#    constraints with the vectorized equivalent of
#    [`MOI.VectorAffineFunction`](@ref) in [`MOI.Nonnegatives`](@ref)
#
#  * Replacing the [`MOI.VectorOfVariables`](@ref) in [`MOI.RootDetConeSquare`](@ref)
#    constraint with [`MOI.VectorAffineFunction`](@ref) in
#    [`MOI.RootDetConeTriangle`](@ref).

model = Model(Clarabel.Optimizer)
set_silent(model)
@variable(model, z[1:n])
@variable(model, s)
@variable(model, t)
## The former @variable(model, Z[1:n, 1:n], PSD)
@variable(model, Z[1:n, 1:n], Symmetric)
@constraint(model, Z >= 0, PSDCone())
## The former [s z'; z Z] >= 0, PSDCone()
@constraint(model, LinearAlgebra.Symmetric([s z'; z Z]) >= 0, PSDCone())
## The former constraint S[i, :]' * Z * S[i, :] - 2 * S[i, :]' * z + s <= 1
f = [1 - S[i, :]' * Z * S[i, :] + 2 * S[i, :]' * z - s for i in 1:m]
@constraint(model, f in MOI.Nonnegatives(m))
## The former constraint [t; vec(Z)] in MOI.RootDetConeSquare(n)
@constraint(model, 1 * [t; triangle_vec(Z)] .+ 0 in MOI.RootDetConeTriangle(n))
## The former @objective(model, Max, t)
@objective(model, Max, 1 * t + 0)
optimize!(model)
assert_is_solved_and_feasible(model)
Test.@test isapprox(D, value.(Z); atol = 1e-3)  #src
solve_time_1 = solve_time(model)

# This formulation gives the much smaller graph:

print_active_bridges(model)

# Note that we still need to bridge [`MOI.PositiveSemidefiniteConeTriangle`](@ref)
# constraints because Clarabel uses the [`MOI.Scaled`](@ref) PSD cone.

model = Model(Clarabel.Optimizer)
set_silent(model)
@variable(model, z[1:n])
@variable(model, s)
@variable(model, t)
@variable(model, Z[1:n, 1:n], Symmetric)
## The former @constraint(model, Z in PSDCone())
f = triangle_vec(Z)
scale_f = [1.0, sqrt(2), 1.0]
@constraint(
    model,
    scale_f .* f in MOI.Scaled(MOI.PositiveSemidefiniteConeTriangle(n)),
)
## The former LinearAlgebra.Symmetric([s z'; z Z]) >= 0, PSDCone()
g = triangle_vec(LinearAlgebra.Symmetric([s z'; z Z]))
scale_g = [1.0, sqrt(2), 1.0, sqrt(2), sqrt(2), 1.0]
@constraint(
    model,
    scale_g .* g in MOI.Scaled(MOI.PositiveSemidefiniteConeTriangle(1 + n)),
)
f = [1 - S[i, :]' * Z * S[i, :] + 2 * S[i, :]' * z - s for i in 1:m]
@constraint(model, f in MOI.Nonnegatives(m))
@constraint(model, 1 * [t; triangle_vec(Z)] .+ 0 in MOI.RootDetConeTriangle(n))
@objective(model, Max, 1 * t + 0)
optimize!(model)
assert_is_solved_and_feasible(model)
Test.@test isapprox(D, value.(Z); atol = 1e-3)  #src
solve_time_2 = solve_time(model)

# This formulation gives the much smaller graph:

print_active_bridges(model)

# Now there is only a single `Unsupported constraint` bullet, showing how JuMP
# reformulated the [`MOI.RootDetConeTriangle`](@ref) constraint by adding a mix
# of [`MOI.PositiveSemidefiniteConeTriangle`](@ref) and
# [`MOI.GeometricMeanCone`](@ref) constraints.

# Because Clarabel doesn't natively support the [`MOI.GeometricMeanCone`](@ref),
# these constraints were further bridged using a
# [`MOI.Bridges.Constraint.GeoMeanToPowerBridge`](@ref)
# to a series of [`MOI.PowerCone`](@ref) constraints.

# However, there are many other ways that a [`MOI.GeometricMeanCone`](@ref) can
# be reformulated into something that Clarabel supports. Let's see what happens
# if we use [`remove_bridge`](@ref) to remove the
# [`MOI.Bridges.Constraint.GeoMeanToPowerBridge`](@ref):

remove_bridge(model, MOI.Bridges.Constraint.GeoMeanToPowerBridge)
optimize!(model)
assert_is_solved_and_feasible(model)

# This time, the solve took:

solve_time_3 = solve_time(model)

# where previously it took

solve_time_2

# Why was the solve time different?

print_active_bridges(model)

# This time, JuMP used a [`MOI.Bridges.Constraint.GeoMeantoRelEntrBridge`](@ref)
# to reformulate the constraint into a set of [`MOI.RelativeEntropyCone`](@ref)
# constraints, which were further reformulated into a set of supported
# [`MOI.ExponentialCone`](@ref) constraints.

# Since the two models are equivalent, we can conclude that for this particular
# model, the formulations have similar performance.

# In general though, the performance of a particular reformulation is problem-
# and solver-specific. Therefore, JuMP chooses to minimize the number of bridges
# in the default reformulation, leaving you to explore alternative formulations
# using the tools and techniques shown in this tutorial.

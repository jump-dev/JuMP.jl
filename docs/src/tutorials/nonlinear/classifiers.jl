# Copyright 2023 James Foster and contributors                                  #src
# This Source Code Form is subject to the terms of the Mozilla Public License   #src
# v.2.0. If a copy of the MPL was not distributed with this file, You can       #src
# obtain one at https://mozilla.org/MPL/2.0/.                                   #src

# # Classifiers

# Classification problems deal with constructing functions, called *classifiers*,
# that can efficiently classify data into two or more distinct sets. 
# A common application is classifying previously unseen data points
# after training a classifier on known data.

# The theory and models in this tutorial come from 
# Section 9.4 "Classification Problems" of the book 
# [*Linear Programming with MATLAB*](https://doi.org/10.1137/1.9780898718775) (2007)
# by M. C. Ferris, O. L. Mangasarian, and S. J. Wright.

# ## Required packages

# This tutorial uses the following packages

using JuMP
import Ipopt
import LinearAlgebra
import Plots
import Random
import Test #src

# ## Data and visualisation

# To start, let's generate some points to test with.
# The argument ``m`` is the number of 2-dimensional points:

function generate_test_points(m; random_seed = 1)
    rng = Random.MersenneTwister(random_seed)
    return 2.0 .* rand(rng, Float64, m, 2)
end

# For the sake of the example, let's take ``m = 100``:

P = generate_test_points(100);

# Note that the points are represented row-wise in the generated array.
# Let's visualise the points using the `Plots` package:

r = 1.01 * maximum(abs.(P))
plot = Plots.scatter(
    P[:, 1],
    P[:, 2];
    xlim = (0, r),
    ylim = (0, r),
    label = nothing,
    color = :white,
    shape = :circle,
    size = (600, 600),
    legend = false,
)

# We want to split the points into two distinct sets on either side of a dividing line.
# We'll then label each point depending on which side of the line it happens to fall.
# Based on the labels of the point, we'll show how to create a classifier using a JuMP model.
# We can then test how well our classifier reproduces the original labels and the boundary between them.

# Let's make a line to divide the point into two sets by defining a gradient and constant:

w0 = [5, 3]
g0 = 8
line(v::AbstractArray; w = w0, g = g0) = LinearAlgebra.dot(w, v) - g
line(x::Real; w = w0, g = g0) = -(w[1] * x - g) / w[2];

# Julia's multiple dispatch feature allows us to define the vector and single-variable form
# of the `line` function under the same name.

# Let's add this to the plot:

Plots.plot!(
    plot,
    x -> line(x),
    0.0:0.01:2.0;
    seriestype = :line,
    linewidth = 5,
    lineopacity = 0.6,
    c = :darkcyan,
)

# Now we label the points relative to which side of the line they are.

P_pos = hcat(filter(v -> line(v) > 0, eachrow(P))...)'
P_neg = hcat(filter(v -> line(v) < 0, eachrow(P))...)'
@assert size(P_pos, 1) + size(P_neg, 1) == size(P, 1) #src
Plots.scatter!(
    plot,
    P_pos[:, 1],
    P_pos[:, 2];
    shape = :cross,
    markercolor = :blue,
    markersize = 8,
)
Plots.scatter!(
    plot,
    P_neg[:, 1],
    P_neg[:, 2];
    shape = :xcross,
    markercolor = :crimson,
    markersize = 8,
)

# The goal is to show we can reconstruct the line from *just* the points and labels.

# ## Formulation: linear support vector machine

# Firstly, we will put the point set back together row-wise as a matrix, with the labelled points group together:

A = [P_pos; P_neg]
m, n = size(A)

# To keep track of the labels, we'll use a diagonal matrix where entry ``i`` of the diagonal is the
# label for point ``i`` (row ``i`` of the matrix).

D = LinearAlgebra.Diagonal([ones(size(P_pos, 1)); -ones(size(P_neg, 1))])

# It is numerically useful to have the labels +1 for `S_pos`  and -1 for `S_neg`.

# A classifier known as the linear *support vector machine* looks for the line function data
# ``w`` and ``g`` such that the function ``L(v) = w^T v - g`` satisfies 
# ``L(p) < 0`` for a point ``p`` in `S_neg` and ``L(p) > 0`` on `S_pos`. 
# The linearly constrained quadratic program that implements this as follows:

# ```math
# \begin{aligned}
# & \min_{w \in \mathbb{R}^n, \; g \in \mathbb{R}, \; y \in \mathbb{R}^m} & \frac{1}{2} w^T w + C \; \sum_{i=1}^m y_i \\
# & \text{subject to} & D_{ii}( A_{i :} w - g ) + y_i & \geq 1, & i = 1 \ldots m \\
# & & y \ge 0.
# \end{aligned}
# ```

# We need a value for the positive penalty parameter ``C``:

C = 1e3;

# ## JuMP formulation

# Now let's build and the JuMP model. We'll get ``w`` and ``g`` from the optimal solution after the
# solve.

m, n = size(A)
model = Model(Ipopt.Optimizer)
set_silent(model)
@variable(model, w[1:n])
@variable(model, g)
@variable(model, y[1:m] >= 0)
@constraint(model, [i in 1:m], D[i, i] * (A[i, :]' * w - g) + y[i] >= 1)
@objective(model, Min, 0.5 * w' * w + C * sum(y))
optimize!(model)
Test.@test termination_status(model) == LOCALLY_SOLVED    #src
Test.@test primal_status(model) == FEASIBLE_POINT  #src
solution_summary(model)

# ## Results

# We recover the solution values

w_sol, g_sol, y_sol = value.(w), value(g), value.(y)
println("Minimum slack: ", minimum(y_sol), "\nMaximum slack: ", maximum(y_sol))

# With the solution, we can ask: was the value of the penalty constant "sufficiently large" 
# for this data set? This can be judged in part by the range of the slack variables.

# Let's add this to the plot as well and check how we did:

Plots.plot!(
    plot,
    x -> line(x; w = w_sol, g = g_sol),
    0.0:0.01:2.0;
    seriestype = :line,
    linewidth = 5,
    linestyle = :dashdotdot,
    lineopacity = 0.9,
    c = :darkblue,
)

# We find that we have recovered the dividing line from just the information of the points
# and their labels.

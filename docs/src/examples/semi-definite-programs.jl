# Copyright 2017, Iain Dunning, Joey Huchette, Miles Lubin, and contributors    #src
# This Source Code Form is subject to the terms of the Mozilla Public License   #src
# v.2.0. If a copy of the MPL was not distributed with this file, You can       #src
# obtain one at https://mozilla.org/MPL/2.0/.                                   #src

# # Semi-definite programs

# These examples use the following packages:

using JuMP
import LinearAlgebra
import Random
import SCS
import Test

# ## K-means clustering

# From "Approximating K-means-type clustering via semidefinite programming" By
# Jiming Peng and Yu Wei.
#
# Given a set of points $a_1, \ldots, a_m$  in $R_n$, allocate them to k clusters.

function example_cluster(; verbose = true)
    ## Data points
    n = 2
    m = 6
    a = Any[
        [2.0, 2.0], [2.5, 2.1], [7.0, 7.0], [2.2, 2.3], [6.8, 7.0], [7.2, 7.5]
    ]
    k = 2
    ## Weight matrix
    W = zeros(m, m)
    for i in 1:m
        for j in i + 1:m
            W[i, j] = W[j, i] = exp(-LinearAlgebra.norm(a[i] - a[j]) / 1.0)
        end
    end
    model = Model(SCS.Optimizer)
    set_silent(model)
    ## Z >= 0, PSD
    @variable(model, Z[1:m, 1:m], PSD)
    @constraint(model, Z .>= 0)
    ## min Tr(W(I-Z))
    I = Matrix(1.0 * LinearAlgebra.I, m, m)
    @objective(model, Min, LinearAlgebra.tr(W * (I - Z)))
    ## Z e = e
    @constraint(model, Z * ones(m) .== ones(m))
    ## Tr(Z) = k
    @constraint(model, LinearAlgebra.tr(Z) == k)
    optimize!(model)
    Z_val = value.(Z)
    ## A simple rounding scheme
    which_cluster = zeros(Int, m)
    num_clusters = 0
    for i in 1:m
        Z_val[i, i] <= 1e-6 && continue
        if which_cluster[i] == 0
            num_clusters += 1
            which_cluster[i] = num_clusters
            for j in i + 1:m
                if LinearAlgebra.norm(Z_val[i, j] - Z_val[i, i]) <= 1e-6
                    which_cluster[j] = num_clusters
                end
            end
        end
    end
    Test.@test which_cluster == [1, 1, 2, 1, 2, 2]  #src
    if verbose
        ## Print results
        for cluster in 1:k
            println("Cluster $cluster")
            for i in 1:m
                if which_cluster[i] == cluster
                    println(a[i])
                end
            end
        end
    end
    return
end

example_cluster()

# ## The correlation problem

# Given three random variables A, B, C and given bounds on two of the three
# correlation coefficients:
#
#     -0.2 <= ρ_AB <= -0.1
#     0.4 <= ρ_BC <=  0.5
#
# We can use the following property of the correlations to determine bounds on
# ρ_AC by solving a SDP:
#
#     |  1    ρ_AB  ρ_AC |
#     | ρ_AB   1    ρ_BC |  ≽ 0
#     | ρ_AC  ρ_BC   1   |

function example_corr_sdp()
    model = Model(SCS.Optimizer)
    set_silent(model)
    @variable(model, X[1:3, 1:3], PSD)
    ## Diagonal is 1s
    @constraint(model, X[1, 1] == 1)
    @constraint(model, X[2, 2] == 1)
    @constraint(model, X[3, 3] == 1)
    ## Bounds on the known correlations
    @constraint(model, X[1, 2] >= -0.2)
    @constraint(model, X[1, 2] <= -0.1)
    @constraint(model, X[2, 3] >=  0.4)
    @constraint(model, X[2, 3] <=  0.5)
    ## Find upper bound
    @objective(model, Max, X[1, 3])
    optimize!(model)
    println("An upper bound for X[1, 3] is $(value(X[1, 3]))")
    Test.@test value(X[1, 3]) ≈ 0.87195 atol = 1e-4  #src
    ## Find lower bound
    @objective(model, Min, X[1, 3])
    optimize!(model)
    println("A lower bound for X[1, 3] is $(value(X[1, 3]))")
    Test.@test value(X[1, 3]) ≈ -0.978 atol = 1e-3  #src
    return
end

example_corr_sdp()

# ## Robust uncertainty sets

# Computes the Value at Risk for a data-driven uncertainty set; see "Data-Driven
# Robust Optimization" (Bertsimas 2013), section 6.1 for details. Closed-form
# expressions for the optimal value are available.

function example_robust_uncertainty()
    R = 1
    d = 3
    𝛿 = 0.05
    ɛ = 0.05
    N = ceil((2 + 2 * log(2 / 𝛿))^2) + 1
    c = randn(d)
    μhat = rand(d)
    M = rand(d, d)
    Σhat = 1 / (d - 1) * (M - ones(d) * μhat')' * (M - ones(d) * μhat')
    Γ1(𝛿, N) = R / sqrt(N) * (2 + sqrt(2 * log(1 / 𝛿)))
    Γ2(𝛿, N) = 2 * R^2 / sqrt(N) * (2 + sqrt(2 * log(2 / 𝛿)))
    model = Model(SCS.Optimizer)
    set_silent(model)
    @variable(model, Σ[1:d, 1:d], PSD)
    @variable(model, u[1:d])
    @variable(model, μ[1:d])
    @constraint(model, [Γ1(𝛿 / 2, N); μ - μhat] in SecondOrderCone())
    @constraint(model, [Γ2(𝛿 / 2, N); vec(Σ - Σhat)] in SecondOrderCone())
    @SDconstraint(model, [((1 - ɛ) / ɛ) (u - μ)'; (u - μ) Σ] >= 0)
    @objective(model, Max, LinearAlgebra.dot(c, u))
    optimize!(model)
    I = Matrix(1.0 * LinearAlgebra.I, d, d)
    exact =
        LinearAlgebra.dot(μhat, c) +
        Γ1(𝛿 / 2, N) * LinearAlgebra.norm(c) +
        sqrt((1 - ɛ) / ɛ) * sqrt(LinearAlgebra.dot(c, (Σhat + Γ2(𝛿 / 2, N) * I) * c))
    Test.@test objective_value(model) ≈ exact atol = 1e-3
    return
end

example_robust_uncertainty()

# ## SDP relaxations: max-cut

# Solves a semidefinite programming relaxation of the MAXCUT graph problem:
#
#     max   0.25 * L•X
#     s.t.  diag(X) == e
#           X ≽ 0
#
# Where `L` is the weighted graph Laplacian. Uses this relaxation to generate a
# solution to the original MAXCUT problem using the method from the paper:
#
# Goemans, M. X., & Williamson, D. P. (1995). Improved approximation algorithms
# for maximum cut and satisfiability problems using semidefinite programming.
# Journal of the ACM (JACM), 42(6), 1115-1145.

function solve_max_cut_sdp(num_vertex, weights)
    ## Calculate the (weighted) Lapacian of the graph: L = D - W.
    laplacian = LinearAlgebra.diagm(0 => weights * ones(num_vertex)) - weights
    ## Solve the SDP relaxation
    model = Model(SCS.Optimizer)
    set_silent(model)
    @variable(model, X[1:num_vertex, 1:num_vertex], PSD)
    @objective(model, Max, 1 / 4 * LinearAlgebra.dot(laplacian, X))
    @constraint(model, LinearAlgebra.diag(X) .== 1)
    optimize!(model)
    ## Compute the Cholesky factorization of X, i.e., X = V^T V.
    opt_X = LinearAlgebra.Hermitian(value.(X), :U)  # Tell Julia its PSD.
    factorization = LinearAlgebra.cholesky(opt_X, Val(true); check = false)
    V = (factorization.P * factorization.L)'
    ## Normalize columns.
    for i in 1:num_vertex
        V[:, i] ./= LinearAlgebra.norm(V[:, i])
    end
    ## Generate random vector on unit sphere.
    Random.seed!(num_vertex)
    r = rand(num_vertex)
    r /= LinearAlgebra.norm(r)
    ## Iterate over vertices, and assign each vertex to a side of cut.
    cut = ones(num_vertex)
    for i in 1:num_vertex
        if LinearAlgebra.dot(r, V[:, i]) <= 0
            cut[i] = -1
        end
    end

    return cut, 0.25 * sum(laplacian .* (cut * cut'))
end

function example_max_cut_sdp()
    ##   [1] --- 5 --- [2]
    ##
    ## Solution:
    ##  (S, S′)  = ({1}, {2})
    cut, cutval = solve_max_cut_sdp(2, [0.0 5.0; 5.0 0.0])
    Test.@test cut[1] != cut[2]
    ##   [1] --- 5 --- [2]
    ##    |  \          |
    ##    |    \        |
    ##    7      6      1
    ##    |        \    |
    ##    |          \  |
    ##   [3] --- 1 --- [4]
    ##
    ## Solution:
    ##  (S, S′)  = ({1}, {2, 3, 4})
    W = [0.0 5.0 7.0 6.0;
         5.0 0.0 0.0 1.0;
         7.0 0.0 0.0 1.0;
         6.0 1.0 1.0 0.0]
    cut, cutval = solve_max_cut_sdp(4, W)
    Test.@test cut[1] != cut[2]
    Test.@test cut[2] == cut[3] == cut[4]
    ##   [1] --- 1 --- [2]
    ##    |             |
    ##    |             |
    ##    5             9
    ##    |             |
    ##    |             |
    ##   [3] --- 2 --- [4]
    ##
    ## Solution:
    ##  (S, S′)  = ({1, 4}, {2, 3})
    W = [0.0 1.0 5.0 0.0;
         1.0 0.0 0.0 9.0;
         5.0 0.0 0.0 2.0;
         0.0 9.0 2.0 0.0]
    cut, cutval = solve_max_cut_sdp(4, W)
    Test.@test cut[1] == cut[4]
    Test.@test cut[2] == cut[3]
    Test.@test cut[1] != cut[2]
    return
end

example_max_cut_sdp()

# ## The minimum distortion problem

# This example arises from computational geometry, in particular the problem of
# embedding a general finite metric space into a euclidean space.
#
# It is known that the 4-point metric space defined by the star graph:
#
#     x
#      \\
#       x — x
#      /
#     x
#
# where distances are computed by length of the shortest path between vertices,
# cannot be exactly embedded into a euclidean space of any dimension.
#
# Here we will formulate and solve an SDP to compute the best possible embedding,
# that is, the embedding f() that minimizes the distortion c such that
#
#     (1 / c) * D(a, b) ≤ ||f(a) - f(b)|| ≤ D(a, b)
#
# for all points (a, b), where D(a, b) is the distance in the metric space.
#
# Any embedding can be characterized by its Gram matrix Q, which is PSD, and
#
#     ||f(a) - f(b)||^2 = Q[a, a] + Q[b, b] - 2 * Q[a, b]
#
# We can therefore constrain
#
#     D[i, j]^2 ≤ Q[i, i] + Q[j, j] - 2 * Q[i, j] ≤ c^2 * D[i, j]^2
#
# and minimize c^2, which gives us the SDP formulation below.
#
# For more detail, see "Lectures on discrete geometry" by J. Matoušek, Springer,
# 2002, pp. 378-379.

function example_min_distortion()
    model = Model(SCS.Optimizer)
    set_silent(model)
    D = [
        0.0 1.0 1.0 1.0;
        1.0 0.0 2.0 2.0;
        1.0 2.0 0.0 2.0;
        1.0 2.0 2.0 0.0
    ]
    @variable(model, c² >= 1.0)
    @variable(model, Q[1:4, 1:4], PSD)
    for i in 1:4
        for j in (i + 1):4
            @constraint(model, D[i, j]^2 <= Q[i, i] + Q[j, j] - 2 * Q[i, j])
            @constraint(model, Q[i, i] + Q[j, j] - 2 * Q[i, j] <= c² * D[i, j]^2)
        end
    end
    @objective(model, Min, c²)
    optimize!(model)
    Test.@test termination_status(model) == MOI.OPTIMAL
    Test.@test primal_status(model) == MOI.FEASIBLE_POINT
    Test.@test objective_value(model) ≈ 4/3 atol = 1e-4
    return
end

example_min_distortion()

# ## Minimum ellipses

# This example is from the Boyd & Vandenberghe book "Convex Optimization". Given
# a set of ellipses centered on the origin
#
#     E(A) = { u | u^T inv(A) u <= 1 }
#
# find a "minimal" ellipse that contains the provided ellipses.
#
# We can formulate this as an SDP:
#
#     minimize    trace(WX)
#     subject to  X >= A_i,    i = 1,...,m
#                 X PSD
#
# where W is a PD matrix of weights to choose between different solutions.

function example_min_ellipse()
    ## We will use three ellipses: two "simple" ones, and a random one.
    As = [
        [2.0  0.0; 0.0  1.0],
        [1.0  0.0; 0.0  3.0],
        [2.86715 1.60645; 1.60645 1.12639]
    ]
    ## We change the weights to see different solutions, if they exist
    weights = [1.0 0.0; 0.0 1.0]
    model = Model(SCS.Optimizer)
    set_silent(model)
    @variable(model, X[i=1:2, j=1:2], PSD)
    @objective(model, Min, LinearAlgebra.tr(weights * X))
    for As_i in As
        @SDconstraint(model, X >= As_i)
    end
    optimize!(model)
    Test.@test termination_status(model) == MOI.OPTIMAL
    Test.@test primal_status(model) == MOI.FEASIBLE_POINT
    Test.@test objective_value(model) ≈ 6.46233 atol = 1e-5
    Test.@test value.(X) ≈ [3.1651 0.8022; 0.8022 3.2972] atol = 1e-4
    return
end

example_min_ellipse()

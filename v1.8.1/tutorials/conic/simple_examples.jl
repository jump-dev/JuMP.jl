# Copyright 2017, Iain Dunning, Joey Huchette, Miles Lubin, and contributors    #src
# This Source Code Form is subject to the terms of the Mozilla Public License   #src
# v.2.0. If a copy of the MPL was not distributed with this file, You can       #src
# obtain one at https://mozilla.org/MPL/2.0/.                                   #src

# # Simple examples

# This tutorial is a collection of examples of small conic programs. It uses
# the following packages:

using JuMP
import LinearAlgebra
import Random
import SCS
import Test

# ## Maximum cut via SDP

# The [maximum cut problem](https://en.wikipedia.org/wiki/Maximum_cut) is a
# classical example in graph theory, where we seek to partition a graph into
# two complementary sets, such that the weight of edges between the two sets is
# maximized. This problem is NP-hard, but it is possible to obtain an
# approximate solution using the semidefinite programming relaxation:
#
#     max   0.25 * L•X
#     s.t.  diag(X) == e
#           X ≽ 0
#
# where `L` is the weighted graph Laplacian. For more details, see:
#
# Goemans, M. X., & Williamson, D. P. (1995). Improved approximation algorithms
# for maximum cut and satisfiability problems using semidefinite programming.
# Journal of the ACM (JACM), 42(6), 1115-1145.

"""
    svd_cholesky(X::AbstractMatrix, rtol)

Return the matrix `U` of the Cholesky decomposition of `X` as `U' * U`.
Note that we do not use the `LinearAlgebra.cholesky` function because it
requires the matrix to be positive definite while `X` may be only
positive *semi*definite.

We use the convention `U' * U` instead of `U * U'` to be consistent with
`LinearAlgebra.cholesky`.
"""
function svd_cholesky(X::AbstractMatrix)
    F = LinearAlgebra.svd(X)
    ## We now have `X ≈ `F.U * D² * F.U'` where:
    D = LinearAlgebra.Diagonal(sqrt.(F.S))
    ## So `X ≈ U' * U` where `U` is:
    return (F.U * D)'
end

function solve_max_cut_sdp(weights)
    N = size(weights, 1)
    ## Calculate the (weighted) Laplacian of the graph: L = D - W.
    L = LinearAlgebra.diagm(0 => weights * ones(N)) - weights
    model = Model(SCS.Optimizer)
    set_silent(model)
    @variable(model, X[1:N, 1:N], PSD)
    for i in 1:N
        set_start_value(X[i, i], 1.0)
    end
    @objective(model, Max, 0.25 * LinearAlgebra.dot(L, X))
    @constraint(model, LinearAlgebra.diag(X) .== 1)
    optimize!(model)
    V = svd_cholesky(value(X))
    Random.seed!(N)
    r = rand(N)
    r /= LinearAlgebra.norm(r)
    cut = [LinearAlgebra.dot(r, V[:, i]) > 0 for i in 1:N]
    S = findall(cut)
    T = findall(.!cut)
    println("Solution:")
    println(" (S, T) = ({", join(S, ", "), "}, {", join(T, ", "), "})")
    return S, T
end

# Given the graph
# ```raw
# [1] --- 5 --- [2]
# ```
# The solution is `(S, T)  = ({1}, {2})`

S, T = solve_max_cut_sdp([0 5; 5 0])

# Given the graph
# ```raw
# [1] --- 5 --- [2]
#  |  \          |
#  |    \        |
#  7      6      1
#  |        \    |
#  |          \  |
# [3] --- 1 --- [4]
# ```
# The solution is `(S, T)  = ({1}, {2, 3, 4})`

S, T = solve_max_cut_sdp([0 5 7 6; 5 0 0 1; 7 0 0 1; 6 1 1 0])

# Given the graph
# ```raw
# [1] --- 1 --- [2]
#  |             |
#  |             |
#  5             9
#  |             |
#  |             |
# [3] --- 2 --- [4]
# ```
# The solution is `(S, T)  = ({1, 4}, {2, 3})`

S, T = solve_max_cut_sdp([0 1 5 0; 1 0 0 9; 5 0 0 2; 0 9 2 0])

# ## K-means clustering via SDP

function example_k_means_clustering()
    a = [[2.0, 2.0], [2.5, 2.1], [7.0, 7.0], [2.2, 2.3], [6.8, 7.0], [7.2, 7.5]]
    m = length(a)
    num_clusters = 2
    W = zeros(m, m)
    for i in 1:m, j in i+1:m
        W[i, j] = W[j, i] = exp(-LinearAlgebra.norm(a[i] - a[j]) / 1.0)
    end
    model = Model(SCS.Optimizer)
    set_silent(model)
    @variable(model, Z[1:m, 1:m] >= 0, PSD)
    @objective(model, Min, LinearAlgebra.tr(W * (LinearAlgebra.I - Z)))
    @constraint(model, [i = 1:m], sum(Z[i, :]) .== 1)
    @constraint(model, LinearAlgebra.tr(Z) == num_clusters)
    optimize!(model)
    Z_val = value.(Z)
    current_cluster, visited = 0, Set{Int}()
    solution = [1, 1, 2, 1, 2, 2]  #src
    for i in 1:m
        if !(i in visited)
            current_cluster += 1
            println("Cluster $current_cluster")
            for j in i:m
                if isapprox(Z_val[i, i], Z_val[i, j]; atol = 1e-3)
                    println(a[j])
                    push!(visited, j)
                    Test.@test solution[j] == current_cluster  #src
                end
            end
        end
    end
    return
end

example_k_means_clustering()

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

function example_correlation_problem()
    model = Model(SCS.Optimizer)
    set_silent(model)
    @variable(model, X[1:3, 1:3], PSD)
    S = ["A", "B", "C"]
    ρ = Containers.DenseAxisArray(X, S, S)
    @constraint(model, [i in S], ρ[i, i] == 1)
    @constraint(model, -0.2 <= ρ["A", "B"] <= -0.1)
    @constraint(model, 0.4 <= ρ["B", "C"] <= 0.5)
    @objective(model, Max, ρ["A", "C"])
    optimize!(model)
    println("An upper bound for ρ_AC is $(value(ρ["A", "C"]))")
    Test.@test value(ρ["A", "C"]) ≈ 0.87195 atol = 1e-4  #src
    @objective(model, Min, ρ["A", "C"])
    optimize!(model)
    println("A lower bound for ρ_AC is $(value(ρ["A", "C"]))")
    Test.@test value(ρ["A", "C"]) ≈ -0.978 atol = 1e-3  #src
    return
end

example_correlation_problem()

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

function example_minimum_distortion()
    model = Model(SCS.Optimizer)
    set_silent(model)
    D = [0.0 1.0 1.0 1.0; 1.0 0.0 2.0 2.0; 1.0 2.0 0.0 2.0; 1.0 2.0 2.0 0.0]
    @variable(model, c² >= 1.0)
    @variable(model, Q[1:4, 1:4], PSD)
    for i in 1:4, j in (i+1):4
        @constraint(model, D[i, j]^2 <= Q[i, i] + Q[j, j] - 2 * Q[i, j])
        @constraint(model, Q[i, i] + Q[j, j] - 2 * Q[i, j] <= c² * D[i, j]^2)
    end
    @objective(model, Min, c²)
    optimize!(model)
    Test.@test termination_status(model) == OPTIMAL
    Test.@test primal_status(model) == FEASIBLE_POINT
    Test.@test objective_value(model) ≈ 4 / 3 atol = 1e-4
    return
end

example_minimum_distortion()

# ## Robust uncertainty sets

# This example computes the Value at Risk for a data-driven uncertainty set;
# see "Data-Driven Robust Optimization" (Bertsimas 2013), Section 6.1 for
# details. Closed-form expressions for the optimal value are available.

function example_robust_uncertainty_sets()
    R, d, 𝛿, ɛ = 1, 3, 0.05, 0.05
    N = ceil((2 + 2 * log(2 / 𝛿))^2) + 1
    c, μhat, M = randn(d), rand(d), rand(d, d)
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
    @constraint(model, [((1-ɛ)/ɛ) (u - μ)'; (u-μ) Σ] >= 0, PSDCone())
    @objective(model, Max, c' * u)
    optimize!(model)
    exact =
        μhat' * c +
        Γ1(𝛿 / 2, N) * LinearAlgebra.norm(c) +
        sqrt((1 - ɛ) / ɛ) *
        sqrt(c' * (Σhat + Γ2(𝛿 / 2, N) * LinearAlgebra.I) * c)
    Test.@test objective_value(model) ≈ exact atol = 1e-2
    return
end

example_robust_uncertainty_sets()

# Copyright 2017, Iain Dunning, Joey Huchette, Miles Lubin, and contributors    #src
# This Source Code Form is subject to the terms of the Mozilla Public License   #src
# v.2.0. If a copy of the MPL was not distributed with this file, You can       #src
# obtain one at https://mozilla.org/MPL/2.0/.                                   #src

# # SDP relaxations: max-cut

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

using JuMP
import LinearAlgebra
import Random
import SCS
import Test

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
    W = [
        0.0 5.0 7.0 6.0
        5.0 0.0 0.0 1.0
        7.0 0.0 0.0 1.0
        6.0 1.0 1.0 0.0
    ]
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
    W = [
        0.0 1.0 5.0 0.0
        1.0 0.0 0.0 9.0
        5.0 0.0 0.0 2.0
        0.0 9.0 2.0 0.0
    ]
    cut, cutval = solve_max_cut_sdp(4, W)
    Test.@test cut[1] == cut[4]
    Test.@test cut[2] == cut[3]
    Test.@test cut[1] != cut[2]
    return
end

example_max_cut_sdp()

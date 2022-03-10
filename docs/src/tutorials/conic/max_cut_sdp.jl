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
import SCS
import Test

function solve_max_cut_sdp(num_vertex, weights)
    ## Calculate the (weighted) Lapacian of the graph: L = D - W.
    laplacian = LinearAlgebra.diagm(0 => weights * ones(num_vertex)) - weights
    ## Solve the SDP relaxation
    model = Model(SCS.Optimizer)
    set_silent(model)
    ## Start with X as the identity matrix to avoid numerical issues.
    @variable(
        model,
        X[i = 1:num_vertex, j = 1:num_vertex],
        PSD,
        start = (i == j ? 1.0 : 0.0),
    )
    @objective(model, Max, 1 / 4 * LinearAlgebra.dot(laplacian, X))
    @constraint(model, LinearAlgebra.diag(X) .== 1)
    optimize!(model)
    @assert termination_status(model) == MOI.OPTIMAL
    ## The paper uses some linear algebra to figure out which nodes belong to
    ## which partition of the cut. We use a simpler way, and one that is more
    ## robust to numerical error: comparing the rows, of which every element is
    ## ±1.
    opt_X = value.(X)
    cut = map(1:num_vertex) do i
        return isapprox(opt_X[1, :], opt_X[i, :]; atol=1e-4) ? 1.0 : -1.0
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

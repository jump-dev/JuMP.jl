#  Copyright 2017, Iain Dunning, Joey Huchette, Miles Lubin, and contributors
#  This Source Code Form is subject to the terms of the Mozilla Public
#  License, v. 2.0. If a copy of the MPL was not distributed with this
#  file, You can obtain one at http://mozilla.org/MPL/2.0/.
#############################################################################
# JuMP
# An algebraic modeling language for Julia
# See http://github.com/JuliaOpt/JuMP.jl
#############################################################################

using JuMP, SCS, LinearAlgebra, Test
"""
    example_cluster()

From "Approximating K-means-type clustering via semidefinite programming" By
Jiming Peng and Yu Wei.

Given a set of points a_1, ..., a_m  in R_n, allocate them to k clusters.
"""
function example_cluster(; verbose = true)
    # Data points
    n = 2
    m = 6
    a = Any[[2.0, 2.0], [2.5, 2.1], [7.0, 7.0],
            [2.2, 2.3], [6.8, 7.0], [7.2, 7.5]]
    k = 2

    # Weight matrix
    W = zeros(m, m)
    for i in 1:m
        for j in i + 1:m
            W[i, j] = W[j, i] = exp(-norm(a[i] - a[j]) / 1.0)
        end
    end

    model = Model(SCS.Optimizer)
    set_silent(model)
    # Z >= 0, PSD
    @variable(model, Z[1:m, 1:m], PSD)
    @constraint(model, Z .>= 0)
    # min Tr(W(I-Z))
    @objective(model, Min, tr(W * (Matrix(1.0I, m, m) - Z)))
    # Z e = e
    @constraint(model, Z * ones(m) .== ones(m))
    # Tr(Z) = k
    @constraint(model, tr(Z) == k)

    JuMP.optimize!(model)

    Z_val = JuMP.value.(Z)

    # A simple rounding scheme
    which_cluster = zeros(Int, m)
    num_clusters = 0
    for i in 1:m
        Z_val[i, i] <= 1e-6 && continue
        if which_cluster[i] == 0
            num_clusters += 1
            which_cluster[i] = num_clusters
            for j in i + 1:m
                if norm(Z_val[i, j] - Z_val[i, i]) <= 1e-6
                    which_cluster[j] = num_clusters
                end
            end
        end
    end

    @test which_cluster == [1, 1, 2, 1, 2, 2]

    if verbose
        # Print results
        for cluster in 1:k
            println("Cluster $cluster")
            for i in 1:m
                if which_cluster[i] == cluster
                    println(a[i])
                end
            end
        end
    end
end

example_cluster(verbose = false)

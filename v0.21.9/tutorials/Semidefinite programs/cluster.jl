# Copyright 2017, Iain Dunning, Joey Huchette, Miles Lubin, and contributors    #src
# This Source Code Form is subject to the terms of the Mozilla Public License   #src
# v.2.0. If a copy of the MPL was not distributed with this file, You can       #src
# obtain one at https://mozilla.org/MPL/2.0/.                                   #src

# # K-means clustering via SDP

# From "Approximating K-means-type clustering via semidefinite programming" By
# Jiming Peng and Yu Wei.
#
# Given a set of points $a_1, \ldots, a_m$  in $R_n$, allocate them to k clusters.

using JuMP
import LinearAlgebra
import SCS
import Test  #src

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

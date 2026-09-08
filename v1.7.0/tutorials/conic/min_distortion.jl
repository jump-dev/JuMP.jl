# Copyright 2017, Iain Dunning, Joey Huchette, Miles Lubin, and contributors    #src
# This Source Code Form is subject to the terms of the Mozilla Public License   #src
# v.2.0. If a copy of the MPL was not distributed with this file, You can       #src
# obtain one at https://mozilla.org/MPL/2.0/.                                   #src

# # The minimum distortion problem

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

using JuMP
import SCS
import Test

function example_min_distortion()
    model = Model(SCS.Optimizer)
    set_silent(model)
    D = [
        0.0 1.0 1.0 1.0
        1.0 0.0 2.0 2.0
        1.0 2.0 0.0 2.0
        1.0 2.0 2.0 0.0
    ]
    @variable(model, c² >= 1.0)
    @variable(model, Q[1:4, 1:4], PSD)
    for i in 1:4
        for j in (i+1):4
            @constraint(model, D[i, j]^2 <= Q[i, i] + Q[j, j] - 2 * Q[i, j])
            @constraint(
                model,
                Q[i, i] + Q[j, j] - 2 * Q[i, j] <= c² * D[i, j]^2
            )
        end
    end
    @objective(model, Min, c²)
    optimize!(model)
    Test.@test termination_status(model) == OPTIMAL
    Test.@test primal_status(model) == FEASIBLE_POINT
    Test.@test objective_value(model) ≈ 4 / 3 atol = 1e-4
    return
end

example_min_distortion()

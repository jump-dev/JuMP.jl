# Copyright 2017, Iain Dunning, Joey Huchette, Miles Lubin, and contributors    #src
# This Source Code Form is subject to the terms of the Mozilla Public License   #src
# v.2.0. If a copy of the MPL was not distributed with this file, You can       #src
# obtain one at https://mozilla.org/MPL/2.0/.                                   #src

# # Minimal ellipses

# This example comes from section 8.4.1 of the book *Convex Optimization* (Boyd and Vandenberghe, 2004).

# Given a set of `m` ellipses `E(A, b, c) = { x : x' A x + 2 b' x + c ≤ 0 }`, we find the ellipse of smallest area that encloses the given ellipses.

# It is convenient to parameterize as the ellipse `{ x : || Px + q || ≤ 1 }`.

# Then the optimal `P` and `q` are given by the convex program

# ````
#       maximize   log(det(P))
      
#     subject to   t[i] ≥ 0, i = 1 ... m
    
#                  [ P^2 - t[i] A[i]   P q - t[i] b[i]       0
#                    P q - t[i] b[i]    -1 - t[i] c[i]    (P q)'
#                                 0              (P q)   - P^2   ]    PSD, i = 1 ... m
# ````
# with helper variables `t[i]`. 

using JuMP
using SCS
using Plots

# Input ellipses, parameterized as
# x' A x + 2 b' x + c ≤ 0
As = [
        [1.2576 -0.3873; -0.3873 0.3467],
        [1.4125 -2.1777; -2.1777 6.7775],
        [1.7018  0.8141;  0.8141 1.7538],
        [0.9742 -0.7202; -0.7202 1.5444],
        [0.6798 -0.1424; -0.1424 0.6871],
        [0.1796 -0.1423; -0.1423 2.6181]
    ]
bs = [
        [ 0.2722,  0.1969],
        [-1.228 , -0.0521],
        [-0.4049,  1.5713],
        [ 0.0265,  0.5623],
        [-0.4301, -1.0157],
        [-0.3286,  0.557 ]
    ]
cs = [0.1831, 0.3295, 0.2077, 0.2362, 0.3284, 0.4931]

function minimal_ellipse(As, bs, cs)
    # Build the model under the change of variables
    #       Psqr = P^2
    #    q_tilde = Pq
    mdl = Model(SCS.Optimizer)

    m = length(As)
    n, _ = size(first(As))

    @variable(mdl, Psqr[1:n, 1:n], PSD)
    @variable(mdl, tau[1:m] ≥ 0)
    @variable(mdl, q_tilde[1:n])
    @variable(mdl, logdetP)
    @constraint(mdl, [logdetP; [Psqr[i, j] for i in 1:n for j in i:n]] in MOI.RootDetConeTriangle(n))

    for (A, b, c, t) in zip(As, bs, cs, tau)
        @constraint(mdl,
            -[      Psqr - t*A     q_tilde - t*b    zeros(n,n)  ;
                (q_tilde - t*b)'        -1 - t*c      q_tilde'  ;
                     zeros(n,n)          q_tilde        -Psqr    ]
            in PSDCone())
    end

    @objective(mdl, Max, logdetP)

    optimize!(mdl)

    @assert termination_status(mdl) == OPTIMAL
    @assert primal_status(mdl) == FEASIBLE_POINT

    # Restore original parameterization
    P = sqrt(value.(Psqr))
    q = P\value.(q_tilde)
    return P, q
end

P, q = minimal_ellipse(As, bs, cs)
@assert isapprox(P, [0.423694 -0.039639; -0.039639 0.316316], atol=1e-4)
@assert isapprox(q, [-0.396177, -0.021368], atol=1e-4)

# Plot results
pl = plot(size=(600,600))
thetas = range(0, 2pi+0.05, step=0.05)

for (A, b, c) in zip(As, bs, cs)
    sqrtA = sqrt(A)
    b_tilde = sqrtA\b
    alpha = b' * (A\b) - c
    rhs = hcat(sqrt(alpha) * cos.(thetas) .- b_tilde[1], sqrt(alpha) * sin.(thetas) .- b_tilde[2])
    ellipse = sqrtA\rhs'
    plot!(pl, ellipse[1, :], ellipse[2, :], label=nothing, c=:navy)
end

plot!(pl,[tuple(P\([cos(theta), sin(theta)] - q) ...) for theta in thetas], c=:crimson, label=nothing)

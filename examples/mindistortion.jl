#############################################################################
# JuMP
# An algebraic modeling langauge for Julia
# See http://github.com/JuliaOpt/JuMP.jl
#############################################################################
# mindistortion.jl
#
# This example arises from computational geometry, in particular the problem
# of embedding a general finite metric space into a euclidean space.
#
# It is known that the 4-point metric space defined by the star graph:
#   x
#    \
#     x — x
#    /
#   x
# where distances are computed by length of the shortest path between
# vertices, cannot be exactly embedded into a euclidean space of any
# dimension.
#
# Here we will formulate and solve an SDP to compute the best possible
# embedding, that is, the embedding f() that minimizes the distortion c such
# that (1/c)*D(a,b) ≤ ||f(a)-f(b)|| ≤ D(a,b) for all points (a,b), where
# D(a,b) is the distance in the metric space.
#
# Any embedding can be characterized by its Gram matrix Q, which is PSD,
# and ||f(a)-f(b)||^2 = Q[a,a] + Q[b,b] - 2Q[a,b]
# We can therefore constrain
# D[i,j]^2 ≤ Q[i,i] + Q[j,j] - 2Q[i,j] ≤ c^2*D[i,j]^2
# and minimize c^2, which gives us the SDP formulation below.
#
# For more detail, see
# "Lectures on discrete geometry" by J. Matoušek, Springer, 2002, pp. 378-379.

using JuMP, SCS


m = Model(solver=SCSSolver())

D = [0.0 1.0 1.0 1.0
     1.0 0.0 2.0 2.0
     1.0 2.0 0.0 2.0
     1.0 2.0 2.0 0.0]

@variable(m, cSq >= 1.0)

@variable(m, Q[1:4,1:4], SDP)

for i in 1:4
    for j in (i+1):4
        @constraint(m, D[i,j]^2 <= Q[i,i] + Q[j,j] - 2Q[i,j])
        @constraint(m, Q[i,i] + Q[j,j] - 2Q[i,j] <= D[i,j]^2*cSq )
    end
end

@objective(m, Min, cSq)

solve(m)

println(getvalue(cSq))

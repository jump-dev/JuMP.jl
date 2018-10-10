##############################################################################
# JuMP
# An algebraic modeling langauge for Julia
# See http://github.com/JuliaOpt/JuMP.jl
#############################################################################
# cluster.jl
#
# From "Approximating K-means-type clustering via semidefinite programming"
# By Jiming Peng and Yu Wei
#
# Given a set of points a_1, ..., a_m  in R_n, allocate them to k clusters
#############################################################################

using JuMP
using SCS

solver = SCSSolver()

# Data points
n = 2
m = 6
a = Any[[2.0, 2.0], [2.5, 2.1], [7.0, 7.0],
        [2.2, 2.3], [6.8, 7.0], [7.2, 7.5]]
k = 2

# Weight matrix
W = zeros(m,m)
for i = 1:m
    for j = i+1:m
        dist = exp(-norm(a[i] - a[j])/1.0)
        W[i,j] = dist
        W[j,i] = dist
    end
end

mod = Model(solver=SCSSolver())

# Z >= 0, PSD
@variable(mod, Z[1:m,1:m], SDP)
@constraint(mod, Z .>= 0)

# min Tr(W(I-Z))
@objective(mod, Min, trace(W * (eye(m,m) - Z)))

# Z e = e
@constraint(mod, Z*ones(m) .== ones(m))
#for i = 1:m
#    addConstraint(mod, sum([ Z[i,j] for j = 1:m]) == 1)
#end

# Tr(Z) = k
@constraint(mod, trace(Z) == k)

solve(mod)

Z_val = getvalue(Z)[:,:]
println("Raw solution")
println(round(Z_val,4))

# A simple rounding scheme
which_cluster = zeros(Int,m)
num_clusters = 0
for i = 1:m
    Z_val[i,i] <= 1e-6 && continue

    if which_cluster[i] == 0
        num_clusters += 1
        which_cluster[i] = num_clusters
        for j = i+1:m
            if norm(Z_val[i,j] - Z_val[i,i]) <= 1e-6
                which_cluster[j] = num_clusters
            end
        end
    end
end

# Print results
for cluster = 1:k
    println("Cluster $cluster")
    for i = 1:m
        if which_cluster[i] == cluster
            println(a[i])
        end
    end
end

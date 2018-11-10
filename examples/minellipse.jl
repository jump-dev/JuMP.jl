#############################################################################
# JuMP
# An algebraic modeleling langauge for Julia
# See http://github.com/JuliaOpt/JuMP.jl
#############################################################################
# minellipse.jl
#
# This example is from the Boyd & Vandenberghe book "Convex Optimization".
# If PyPlot is installed, passing the argument `plot` will plot
# the solution, e.g. julia minellipse.jl plot
#
# Given a set of ellipses centered on the origin
# E(A) = { u | u^T inv(A) u <= 1 }
# find a "minimal" ellipse that contains the provided ellipses
#
# We can formulate this as an SDP:
#     minimize  trace(WX)
#   subject to  X >= A_i,    i = 1,...,m
#               X PSD
# where W is a PD matrix of weights to choose between different solutions
#############################################################################

using JuMP
using SCS
using LinearAlgebra
#using PyPlot  # Comment out if not installed

# solver = SCSSolver(eps=1e-6)

# We will use three ellipses
m = 3
# Two "simple" ones
As = Any[ [2.0  0.0;
           0.0  1.0],
          [1.0  0.0;
           0.0  3.0]]
# and a random one
randA = rand(2,2)
push!(As, (randA' * randA) * (rand()*2+1))

# We change the weights to see different solutions, if they exist
W = [1.0 0.0
     0.0 1.0];

model = Model(with_optimizer(SCS.Optimizer))
@variable(model, X[1:2, 1:2], PSD)
@objective(model, Min, tr(W*X))
for i in 1:m
    @SDconstraint(model, X >= As[i])
end
JuMP.optimize!(model)

X_val = JuMP.value.(X)
println(X_val)

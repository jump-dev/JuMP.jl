#############################################################################
# JuMP
# An algebraic modeling langauge for Julia
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
#using PyPlot  # Comment out if not installed

solver = SCSSolver(eps=1e-6)

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

mod = Model(solver=solver)
@variable(mod, X[1:2,1:2], SDP)
@objective(mod, Min, trace(W*X))
for i = 1:m
    @SDconstraint(mod, X >= As[i])
end
solve(mod)

X_val = getvalue(X)
println(X_val)

#= This code no longer seems to be working.
ERROR: LoadError: PyError (:PyObject_Call) <type 'exceptions.ValueError'>
ValueError('The truth value of an array with more than one element is ambiguous. Use a.any() or a.all()',)

# Plot it, if desired (e.g. julia minellipse.jl plot)
if length(ARGS) > 0 && Pkg.installed("PyPlot") !== nothing
    # Setup the figure
    fig = figure()
    ax = fig[:gca]()
    ax[:set_xticks]([-4:+4])
    ax[:set_yticks]([-4:+4])

    # Draw provided ellipses
    for i = 1:m
        xs = Float64[]
        ys = Float64[]
        for angle in linspace(0, 2*pi, 100)
            u = [cos(angle),sin(angle)]
            x = As[i] * u
            push!(xs, x[1])
            push!(ys, x[2])
        end
        plot(xs, ys, "b", linewidth=2.0)
    end

    # Draw bounding ellipse
    xs = Float64[]
    ys = Float64[]
    for angle in linspace(0, 2*pi, 100)
        u = [cos(angle),sin(angle)]
        x = X_val * u
        push!(xs, x[1])
        push!(ys, x[2])
    end
    plot(xs, ys, "r", linewidth=2.0)

    fig[:canvas][:draw]()
    readline()
end
=#

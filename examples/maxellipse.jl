#############################################################################
# JuMP
# An algebraic modelling langauge for Julia
# See http://github.com/JuliaOpt/JuMP.jl
#############################################################################
# maxellipse.jl
#
# This example is from the Mosek modelling manual. 
# If PyPlot is installed, passing the argument `plot` will plot
# the solution, e.g. julia maxellipse.jl plot
#
# Given a polytope
# P = { x :  a_i x <= b_i   i = 1...m}
# we seek to find the maximal inscribed ellipse 
# E = { x :  x = Cu + d, || u || <= 1 }
#
# The ellipse is contained in P iff
# max  a_i^T (Cu + d) <= b_i  for all i = 1...m
#  u 
# which we can rewrite (e.g. using duality, like in robust optimization) as
# ||C a_i|| + a_i^T d <= b_i  for all i = 1...m
#
# The volume of this ellipse can be approximated by det(C) ^ (1/n)
#
# So we can solve the following problem 
#     maximize   det(C) ^ (1/n)
#   subject to   ||C a_i|| + a_i^T d <= b_i  for all i = 1...m
#                C is p.s.d.
# 
# Alternatively, we can solve
#     maximize   t
#   subject to   ||C a_i|| + a_i^T d <= b_i  for all i = 1...m
#                [ C    Z       ] 
#                [ Z^T  Diag(Z) ] PSD
#                t <= (Z_11 * Z_22 * ... * Z_nn) ^ 1/n
# which is an SDP
#############################################################################

using JuMP
using Mosek
using PyPlot  # Comment out if not installed

# 1  ^  2
#   /+\
#  /---\ 
#    3
#m = 3
#a = {[-1.0, +1.0],   
#     [+1.0, +1.0],
#     [ 0.0, -1.0]}
#b = [1.0, 1.0, 0.0]

# Same as above but with top of triangle chopped off
m = 4
a = {[-1.0, +1.0],   
     [+1.0, +1.0],
     [ 0.0, -1.0],
     [ 0.0,  1.0]}
b = [1.0, 1.0, 0.0, 0.7]

mod = Model()

@defVar(mod, 0 <= t <= 100)
@setObjective(mod, Max, t)

@defVar(mod, -1 <= d[1:2] <= 1)

# ||C a_i|| + a_i^T d <= b_i  for all i = 1...m
@defSDPVar(mod, C[2,2])
for i = 1:m
    addConstraint(mod, norm(C*a[i]) + dot(a[i],d) <= b[i])
end

# [ C    Z       ] 
# [ Z^T  Diag(Z) ] PSD
@defMatrixVar(mod, Z[2,2])
lhs = [C Z; Z' diagm(Z)]
addConstraint(mod, lhs >= 0)

# t <= (Z_11 * Z_22 * ... * Z_nn) ^ 1/n
# t <= (Z_11 * Z_22) ^ 1/2
# This in general requires rotated cone shenanigans
# but for the 2D case is simple enough.
addConstraint(mod, t^2 <= Z[1,1] * Z[2,2])

solve(mod)

C_val = getValue(C)[:,:]
d_val = getValue(d)[:]

println(C_val)

# Plot it, if desired (e.g. julia maxellipse.jl plot)
if length(ARGS) > 0 && Pkg.installed("PyPlot") != nothing
    # Setup the figure
    fig = figure()
    ax = fig[:gca]()
    ax[:set_xticks]([-2:+2])
    ax[:set_yticks]([-2:+2])

    # Draw polyhedron
    for i = 1:m
      x1 = -2
      y1 = (b[i] - a[i][1]*x1)/a[i][2]
      x2 = +2
      y2 = (b[i] - a[i][1]*x2)/a[i][2]
      plot((x1, x2), (y1, y2), "k", linewidth=2.0)
    end

    # Draw ellipse
    # x = Cu + d, || u || <= 1
    xs = Float64[]
    ys = Float64[]
    for angle in linspace(0, 2*pi, 100)
        u = [cos(angle),sin(angle)]
        x = C_val * u + d_val
        push!(xs, x[1])
        push!(ys, x[2])
    end
    plot(xs, ys, "b", linewidth=2.0)

    fig[:canvas][:draw]()
    readline()
end
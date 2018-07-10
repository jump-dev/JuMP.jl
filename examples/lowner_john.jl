# 
# Löwner-John ellipsoid example
#
# Compute the inner and outer ellipsoids of a polytope.
#
# Adapter from example in the MOSEK distro.
#
# The inner ellipsoidal approximation to a polytope
#
#    S = { x \in R^n | Ax < b }.
#
# maximizes the volume of the inscribed ellipsoid,
#
#    { x | x = C*u + d, || u ||_2 <= 1 }.
#
# The volume is proportional to det(C)^(1/n), so the
# problem can be solved as
#
#   maximize         t
#   subject to       t       <= det(C)^(1/n)
#               || C*ai ||_2 <= bi - ai^T * d,  i=1,...,m
#                 C is PSD
#
# which is equivalent to a mixed conic quadratic and semidefinite
# programming problem.
#
#
# The outer ellipsoidal approximation to a polytope given
# as the convex hull of a set of points
#
#     S = conv{ x1, x2, ... , xm }
#
# minimizes the volume of the enclosing ellipsoid,
#
#   { x | || P*(x-c) ||_2 <= 1 }
#
# The volume is proportional to det(P)^{-1/n}, so the problem can
# be solved as
#
#   minimize         t
#   subject to       t       >= det(P)^(-1/n)
#               || P*xi + c ||_2 <= 1,  i=1,...,m
#                 P is PSD.
#
# References:
# [1] "Lectures on Modern Optimization", Ben-Tal and Nemirovski, 2000.
#

using JuMP
import MathOptInterface
const MOI = MathOptInterface

# Models the convex set 
#
#   S = { (x, t) \in R^n x R | x >= 0, t <= (x1 * x2 * ... * xn)^(1/n) }
#
function geometric_mean(m, x, t)
    n = length(x)
    if n == 1
        @constraint(m,t-x[1] <= 0)
    else
        t2 = @variable(m)
        @constraint(m,1.0*[t2; x[n]; t] in MOI.PowerCone(1-1.0/n))
        geometric_mean(m,x[1:n-1], t2)
    end
end

# Purpose: Models the hypograph of the n-th power of the
# determinant of a positive definite matrix. See [1,2] for more details.
#
#   The convex set (a hypograph)
#
#   C = { (X, t) \in S^n_+ x R |  t <= det(X)^{1/n} },
#
#   can be modeled as the intersection of a semidefinite cone
#
#   [ X, Z; Z^T Diag(Z) ] >= 0  
#
#   and a number of rotated quadratic cones and affine hyperplanes,
#
#   t <= (Z11*Z22*...*Znn)^{1/n}  (see geometric_mean).
function det_root(m,X,t)
    n = floor(Int,sqrt(length(X)))

    @variable(m,Y[1:2*n,1:2*n],PSD)
    @constraint(m,reshape(Y[n+1:2*n,1:n] - Y[n+1:2*n,n+1:2*n],n*n) in MOI.Zeros(n*n))
    @constraint(m,reshape(X - Y[1:n,1:n], n*n) in MOI.Zeros(n*n))

    geometric_mean(m,[Y[i+n,i+n] for i in 1:n],t)
end

#  The inner ellipsoidal approximation to a polytope 
#
#     S = { x \in R^n | Ax < b }.
#
#  maximizes the volume of the inscribed ellipsoid,
#
#     { x | x = C*u + d, || u ||_2 <= 1 }.
#
#  The volume is proportional to det(C)^(1/n), so the
#  problem can be solved as 
#
#    maximize         t
#    subject to       t       <= det(C)^(1/n)
#                || C*ai ||_2 <= bi - ai^T * d,  i=1,...,m
#                C is PSD
#
#  which is equivalent to a mixed conic quadratic and semidefinite
#  programming problem.
function lownerjohn_inner(solver :: MOI.AbstractOptimizer,
                          A      :: Array{Float64,2},
                          b      :: Vector{Float64})
    m,n = size(A)

    model = Model(mode=JuMP.Direct, backend=solver)

    @variable(model, t >= 0)
    @variable(model, C[1:n,1:n], PSD)
    @variable(model, d[1:n])

    # [ b-A*d ; A*C ] in QCone

    for i in 1:m
        @constraint(model, [ b[i] - (A[i,:]' * d) ; (A[i,:]'*C)'] in MOI.SecondOrderCone(n+1))
    end

    det_root(model,C,t)

    @objective(model, Max, t)

    optimize(model)

    reshape([JuMP.resultvalue(item) for item in C],(n,n)),[JuMP.resultvalue(item) for item in d]
end

# The outer ellipsoidal approximation to a polytope given 
# as the convex hull of a set of points
#
#   S = conv{ x1, x2, ... , xm }
#
# minimizes the volume of the enclosing ellipsoid,
#
#   { x | || P*x-c ||_2 <= 1 }
#
# The volume is proportional to det(P)^{-1/n}, so the problem can
# be solved as
#
#   maximize         t
#   subject to       t       <= det(P)^(1/n)
#               || P*xi - c ||_2 <= 1,  i=1,...,m
#               P is PSD.
function lownerjohn_outer(solver :: MOI.AbstractOptimizer,
                          x      :: Array{Float64,2})
    m,n = size(x)
    model = Model(mode=JuMP.Direct, backend=solver)

    @variable(model,t >= 0)
    @variable(model,P[1:n,1:n],PSD)
    @variable(model,c[1:n])

    # (1, Px-c) in cone
    for i in 1:m
        @constraint(model,[1 ; (x[i,:]' * P)' - c] in MOI.SecondOrderCone(n+1))
    end

    # t <= det(P)^{1/n}
    det_root(model,P,t)

    @objective(model,Max,t)

    optimize(model)

    [ JuMP.resultvalue(item) for item in P ], [ JuMP.resultvalue(item) for item in c ]
end





# Probably do something else here:
import MathOptInterfaceMosek
solver = MathOptInterfaceMosek.MosekOptimizer

# Test

# Vertices of a pentagon in 2D
p = [0.0  0.0;
     1.0  3.0;
     5.5  4.5;
     7.0  4.0;
     7.0  1.0;
     3.0 -2.0]

n = size(p,1)
# The hyperplane representation of the same polytope
A = [ (-p[:,2] + [p[n,2];p[1:n-1,2]])  (p[:,1] - [p[n,1];p[1:n-1,1]]) ]
b = A[:,1] .* p[:,1] + A[:,2] .* p[:,2]

Po, co = lownerjohn_outer(solver(), p)
Ci, di = lownerjohn_inner(solver(), A, b)
Poinv = Po ^ -1


println("Solution:")
println("  Outer ellipsoid: Po = $Po, co = $co")
println("  Inner ellipsoid: Ci = $Ci, di = $di")


if Pkg.installed("Plots") != nothing && Pkg.installed("GR") != nothing
    println("Plotting solution...")
    import Plots
    Plots.gr()
    Plots.plot([p[:,1] ; p[1,1]], [p[:,2] ; p[1,2]])    

    Plots.plot!(t -> Ci[1,:]' * [ cos(t),sin(t) ] + di[1],
                t -> Ci[2,:]' * [ cos(t),sin(t) ] + di[2],
                0,2*pi)

    Plots.plot!(t -> Poinv[1,:]' * ( [cos(t),sin(t)] + c ),
                t -> Poinv[2,:]' * ( [cos(t),sin(t)] + c ),
                0,2*pi)
    
    println("Plot done. Press enter to continue.")
    readline()
else
    println("Install Plots package to plot solution")    
end

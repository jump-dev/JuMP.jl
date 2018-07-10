using JuMP
import MathOptInterface
import MathOptInterfaceMosek
const MOI = MathOptInterface


#
# In the following two examples V is the factored covariance matrix,
# i.e. Q = VV', μ is the vector of expected returns for a set of
# stocks, δ is the target expected return and S is a measure for the
# risk (specifically, the maximum variance).
#
# The interesting constraint here is the constraint t > x'Q x. We notice that
#   x' Q x = x'VVx = (x'V)^2 = || x' V ||^2
# And write
#   t > || x' V ||^2
# Or equivalently
#   (0.5, t, x' V) in RotatedSecondOrderCone
#
# Our initial total wealth is 1, thus sum(x)=1.
#

# Find the position x with minimal risk such that we get at least thge
# expected return δ.
#
#
function portfolio1(solver, V, μ, δ)
    model = Model(optimizer = solver())
    n = size(μ,1)
    m = size(V,2)

    
    @variable(model, t >= 0)
    @variable(model, x[1:n] >= 0)

    @constraint(model, μ'x .>= δ)

    @constraint(model, sum(x) == 1.0)
    println("$(size(V)), $(size(μ)), $(size(x))")
    @constraint(model, [0.5; t; V' * x] in MOI.RotatedSecondOrderCone(m+2))

    @objective(model, Min, t)

    optimize(model)

    JuMP.resultvalue(t), [ JuMP.resultvalue(item) for item in x ]
end

# Find the position x with maximum return such that the variance is at most σ.
function portfolio2(solver, V, μ, R)
    model = Model(optimizer = solver())
    n = length(μ)
    m = size(V,2)
    
    @variable(model, t >= 0)
    @variable(model, x[1:n] >= 0)

    @constraint(model, sum(x) == 1.0)
    @constraint(model, [0.5; t; V' * x] in MOI.RotatedSecondOrderCone(m+2))
    @constraint(model, t <= R)

    @objective(model, Max, μ'x)

    optimize(model)

    JuMP.resultvalue(t), [ JuMP.resultvalue(item) for item in x ]
end




solver = MathOptInterfaceMosek.MosekOptimizer
## Example data:

# Vector of expected returns
μ = [1.05, 1.3, 0.9, 1.0, 0.9, 1.5]

# lower bound on return of investment
δ = 1.2

# upper bound for standard deviation
σ = 0.018

# factor of covariance (we use the data-matrix)
V = [-0.0644  -0.0352  -0.0019   0.1015 ;
      0.1183  -0.0428  -0.0013  -0.0741 ;
     -0.0798  -0.1201   0.1040   0.0958 ;
     -0.0855  -0.1174   0.1042   0.0987 ;
     -0.0158  -0.0074   0.0110   0.0123 ;
      0.0636  -0.1341   0.0206   0.0497 ]


t1,x1 = portfolio1(solver,V,μ,δ)
t2,x2 = portfolio2(solver,V,μ,σ)

println("Minimize risk with δ = $δ:")
println("  variance = $t1")
println("  x = $x1")
println("  expected return = $(μ'x1)")
println("")
println("Maximize expected return with std.dev <= $σ:")
println("  variance = $t2")
println("  x = $x2")
println("  expected return = $(μ'x2)")

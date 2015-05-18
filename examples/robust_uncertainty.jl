#############################################################################
# JuMP
# An algebraic modelling langauge for Julia
# See http://github.com/JuliaOpt/JuMP.jl
#############################################################################
# robust_uncertainty.jl
#
# Computes the Value at Risk for a data-driven uncertainty set; see
# "Data-Driven Robust Optimization" (Bertsimas 2013), section 6.1 for
# details. Closed-form expressions for the optimal value are available.
#############################################################################


using JuMP

R = 1
d = 3
𝛿 = 0.05
ɛ = 0.05
N = ceil((2+2log(2/𝛿))^2) + 1

Γ1(𝛿,N) = (R/sqrt(N))*(2+sqrt(2*log(1/𝛿)))
Γ2(𝛿,N) = (2R^2/sqrt(N))*(2+sqrt(2*log(2/𝛿)))

μhat = rand(d)
M = rand(d,d)
# Σhat = 1/(d-1)*(M-ones(d)*μhat')'*(M-ones(d)*μhat')
Σhat = (M-ones(d)*μhat')'*(M-ones(d)*μhat')

m = Model()

@defVar(m, Σ[1:d,1:d], SDP)
@defVar(m, u[1:d])
@defVar(m, μ[1:d])

@defVar(m, t1 >= 0)
@defVar(m, L1[1:d])
@addConstraint(m, L1 .== (μ-μhat))
@addConstraint(m, sum{L1[i]^2, i=1:d} <= t1^2)
@addConstraint(m, t1 <= Γ1(𝛿/2,N))

@defVar(m, t2 >= 0)
@defVar(m, L2[1:d,1:d])
@addConstraint(m, L2 .== (Σ-Σhat))
@addConstraint(m, sum{L2[i,j]^2, i=1:d, j=1:d} <= t2^2)
@addConstraint(m, t2 <= Γ2(𝛿/2,N))

A = [(1-ɛ)/ɛ (u-μ)';
     (u-μ)     Σ   ]
@addSDPConstraint(m, A >= 0)

c = randn(d)
@setObjective(m, Max, dot(c,u))

solve(m)

object = getObjectiveValue(m)
exact = dot(μhat,c) + Γ1(𝛿/2,N)*norm(c) + sqrt((1-ɛ)/ɛ)*sqrt(dot(c,(Σhat+Γ2(𝛿/2,N)*eye(d,d))*c))

println("objective value:  $(object)")
println("error from exact: $(abs(exact-object))")
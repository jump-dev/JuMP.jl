#############################################################################
# JuMP
# An algebraic modeleling langauge for Julia
# See http://github.com/JuliaOpt/JuMP.jl
#############################################################################

using JuMP, SCS, LinearAlgebra, Test

#=
This example is from the Boyd & Vandenberghe book "Convex Optimization". Given a
set of ellipses centered on the origin
    E(A) = { u | u^T inv(A) u <= 1 }
find a "minimal" ellipse that contains the provided ellipses.

We can formulate this as an SDP:
    minimize  trace(WX)
  subject to  X >= A_i,    i = 1,...,m
              X PSD
where W is a PD matrix of weights to choose between different solutions.
=#

function example_minellipse(; verbose = true)
    # We will use three ellipses: two "simple" ones, and a random one.
    rand_A = rand(2, 2)
    As = [
        [2.0  0.0; 0.0  1.0],
        [1.0  0.0; 0.0  3.0],
        (rand_A' * rand_A) * (2 * rand() + 1)
    ]
    # We change the weights to see different solutions, if they exist
    weights = [1.0 0.0; 0.0 1.0]
    model = Model(with_optimizer(SCS.Optimizer))
    @variable(model, X[1:2, 1:2], PSD)
    @objective(model, Min, tr(weights * X))
    for As_i in As
        @SDconstraint(model, X >= As_i)
    end
    JuMP.optimize!(model)

    X_val = JuMP.value.(X)
    if verbose
        println(X_val)
    end
    @test JuMP.termination_status(model) == MOI.OPTIMAL
    @test JuMP.primal_status(model) == MOI.FEASIBLE_POINT
    @test JuMP.objective_value(model) ≈ 5.0 atol = 1e-6
    @test JuMP.value.(X) ≈ [2.0 0.0; 0.0 3.0] atol = 1e-6
end

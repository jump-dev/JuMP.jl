#############################################################################
# JuMP
# An algebraic modeleling langauge for Julia
# See http://github.com/JuliaOpt/JuMP.jl
#############################################################################

using JuMP, SCS, LinearAlgebra, Test

"""
    example_minellipse()

This example is from the Boyd & Vandenberghe book "Convex Optimization". Given a
set of ellipses centered on the origin
    E(A) = { u | u^T inv(A) u <= 1 }
find a "minimal" ellipse that contains the provided ellipses.

We can formulate this as an SDP:
    minimize  trace(WX)
  subject to  X >= A_i,    i = 1,...,m
              X PSD
where W is a PD matrix of weights to choose between different solutions.
"""
function example_minellipse()
    # We will use three ellipses: two "simple" ones, and a random one.
    As = [
        [2.0  0.0; 0.0  1.0],
        [1.0  0.0; 0.0  3.0]  #,
        # [2.86715 1.60645; 1.60645 1.12639]
    ]
    # We change the weights to see different solutions, if they exist
    weights = [1.0 0.0; 0.0 1.0]
    model = Model(with_optimizer(SCS.Optimizer, verbose = 0))
    @variable(model, X[i=1:2, j=1:2], PSD, start = [2.0 0.0; 0.0 3.0][i, j])
    @objective(model, Min, tr(weights * X))
    for As_i in As
        @SDconstraint(model, X >= As_i)
    end
    JuMP.optimize!(model)

    @test JuMP.termination_status(model) == MOI.OPTIMAL
    @test JuMP.primal_status(model) == MOI.FEASIBLE_POINT
    @test JuMP.objective_value(model) ≈ 5.0 atol = 1e-4
    @test JuMP.value.(X) ≈ [2.0 0.0; 0.0 3.0] atol = 1e-4
end

example_minellipse()

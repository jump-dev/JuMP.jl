#  Copyright 2017, Iain Dunning, Joey Huchette, Miles Lubin, and contributors
#  This Source Code Form is subject to the terms of the Mozilla Public
#  License, v. 2.0. If a copy of the MPL was not distributed with this
#  file, You can obtain one at http://mozilla.org/MPL/2.0/.
#############################################################################
# JuMP
# An algebraic modeling language for Julia
# See http://github.com/JuliaOpt/JuMP.jl
#############################################################################

using JuMP, SCS, LinearAlgebra, Test
const MOI = JuMP.MathOptInterface

"""
    example_min_ellipse()

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
function example_min_ellipse()
    # We will use three ellipses: two "simple" ones, and a random one.
    As = [
        [2.0  0.0; 0.0  1.0],
        [1.0  0.0; 0.0  3.0],
        [2.86715 1.60645; 1.60645 1.12639]
    ]
    # We change the weights to see different solutions, if they exist
    weights = [1.0 0.0; 0.0 1.0]
    model = Model(SCS.Optimizer)
    set_silent(model)
    @variable(model, X[i=1:2, j=1:2], PSD, start = [2.0 0.0; 0.0 3.0][i, j])
    @objective(model, Min, tr(weights * X))
    for As_i in As
        @SDconstraint(model, X >= As_i)
    end
    JuMP.optimize!(model)

    @test JuMP.termination_status(model) == MOI.OPTIMAL
    @test JuMP.primal_status(model) == MOI.FEASIBLE_POINT
    @test JuMP.objective_value(model) ≈ 6.46233 atol = 1e-5
    @test JuMP.value.(X) ≈ [3.1651 0.8022; 0.8022 3.2972] atol = 1e-4
end

example_min_ellipse()

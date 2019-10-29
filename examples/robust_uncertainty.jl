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

"""
    example_robust_uncertainty()

Computes the Value at Risk for a data-driven uncertainty set; see "Data-Driven
Robust Optimization" (Bertsimas 2013), section 6.1 for details. Closed-form
expressions for the optimal value are available.
"""
function example_robust_uncertainty()
    R = 1
    d = 3
    𝛿 = 0.05
    ɛ = 0.05
    N = ceil((2 + 2 * log(2 / 𝛿))^2) + 1

    c = randn(d)

    μhat = rand(d)
    M = rand(d, d)
    Σhat = 1 / (d - 1) * (M - ones(d) * μhat')' * (M - ones(d) * μhat')

    Γ1(𝛿, N) = R / sqrt(N) * (2 + sqrt(2 * log(1 / 𝛿)))
    Γ2(𝛿, N) = 2 * R^2 / sqrt(N) * (2 + sqrt(2 * log(2 / 𝛿)))

    model = Model(SCS.Optimizer)
    set_silent(model)

    @variable(model, Σ[1:d, 1:d], PSD)
    @variable(model, u[1:d])
    @variable(model, μ[1:d])
    @constraint(model, [Γ1(𝛿 / 2, N); μ - μhat] in SecondOrderCone())
    @constraint(model, [Γ2(𝛿 / 2, N); vec(Σ - Σhat)] in SecondOrderCone())
    @SDconstraint(model, [((1 - ɛ) / ɛ) (u - μ)'; (u - μ) Σ] >= 0)
    @objective(model, Max, dot(c, u))

    JuMP.optimize!(model)

    exact = dot(μhat, c) + Γ1(𝛿 / 2, N) * norm(c) + sqrt((1 - ɛ) / ɛ) *
        sqrt(dot(c, (Σhat + Γ2(𝛿 / 2, N) * Matrix(1.0I, d, d)) * c))
    @test JuMP.objective_value(model) ≈ exact atol = 1e-3
end

example_robust_uncertainty()

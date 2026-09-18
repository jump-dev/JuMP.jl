# Copyright 2017, Iain Dunning, Joey Huchette, Miles Lubin, and contributors    #src
# This Source Code Form is subject to the terms of the Mozilla Public License   #src
# v.2.0. If a copy of the MPL was not distributed with this file, You can       #src
# obtain one at https://mozilla.org/MPL/2.0/.                                   #src

# # Robust uncertainty sets

# Computes the Value at Risk for a data-driven uncertainty set; see "Data-Driven
# Robust Optimization" (Bertsimas 2013), section 6.1 for details. Closed-form
# expressions for the optimal value are available.

using JuMP
import SCS
import LinearAlgebra
import Test

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
    @constraint(model, [((1-ɛ)/ɛ) (u - μ)'; (u-μ) Σ] in PSDCone())
    @objective(model, Max, LinearAlgebra.dot(c, u))
    optimize!(model)
    I = Matrix(1.0 * LinearAlgebra.I, d, d)
    exact =
        LinearAlgebra.dot(μhat, c) +
        Γ1(𝛿 / 2, N) * LinearAlgebra.norm(c) +
        sqrt((1 - ɛ) / ɛ) *
        sqrt(LinearAlgebra.dot(c, (Σhat + Γ2(𝛿 / 2, N) * I) * c))
    Test.@test objective_value(model) ≈ exact atol = 1e-2
    return
end

example_robust_uncertainty()

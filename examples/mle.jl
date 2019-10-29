#  Copyright 2017, Iain Dunning, Joey Huchette, Miles Lubin, and contributors
#  This Source Code Form is subject to the terms of the Mozilla Public
#  License, v. 2.0. If a copy of the MPL was not distributed with this
#  file, You can obtain one at http://mozilla.org/MPL/2.0/.
#############################################################################
# JuMP
# An algebraic modeling language for Julia
# See http://github.com/JuliaOpt/JuMP.jl
#############################################################################

using JuMP, Ipopt, Random, Statistics, Test

"""
    example_mle()

Use nonlinear optimization to compute the maximum likelihood estimate (MLE) of
the parameters of a normal distribution aka the sample mean and variance
"""
function example_mle(; verbose = true)
    n = 1_000
    Random.seed!(1234)
    data = randn(n)
    model = Model(Ipopt.Optimizer)
    set_silent(model)
    @variable(model, μ, start = 0.0)
    @variable(model, σ >= 0.0, start = 1.0)
    @NLobjective(model, Max, n / 2 * log(1 / (2 * π * σ^2)) -
        sum((data[i] - μ)^2 for i in 1:n) / (2 * σ^2)
    )
    JuMP.optimize!(model)
    if verbose
        println("μ = ", JuMP.value(μ))
        println("mean(data) = ", mean(data))
        println("σ^2 = ", JuMP.value(σ)^2)
        println("var(data) = ", var(data))
        println("MLE objective: ", JuMP.objective_value(model))
    end
    @test JuMP.value(μ) ≈ mean(data) atol = 1e-3
    @test JuMP.value(σ)^2 ≈ var(data) atol = 1e-2

    # constrained MLE?
    @NLconstraint(model, μ == σ^2)

    JuMP.optimize!(model)
    if verbose
        println("\nWith constraint μ == σ^2:")
        println("μ = ", JuMP.value(μ))
        println("σ^2 = ", JuMP.value(σ)^2)
        println("Constrained MLE objective: ", JuMP.objective_value(model))
    end
    @test JuMP.value(μ) ≈ JuMP.value(σ)^2
end

example_mle(verbose = false)

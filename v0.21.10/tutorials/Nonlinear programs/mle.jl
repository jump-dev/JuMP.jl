# Copyright 2017, Iain Dunning, Joey Huchette, Miles Lubin, and contributors    #src
# This Source Code Form is subject to the terms of the Mozilla Public License   #src
# v.2.0. If a copy of the MPL was not distributed with this file, You can       #src
# obtain one at https://mozilla.org/MPL/2.0/.                                   #src

# # Maximum likelihood estimation

# Use nonlinear optimization to compute the maximum likelihood estimate (MLE) of
# the parameters of a normal distribution, a.k.a., the sample mean and variance.

using JuMP
import Ipopt
import Random
import Statistics
import Test

function example_mle(; verbose = true)
    n = 1_000
    Random.seed!(1234)
    data = randn(n)
    model = Model(Ipopt.Optimizer)
    set_silent(model)
    @variable(model, μ, start = 0.0)
    @variable(model, σ >= 0.0, start = 1.0)
    @NLobjective(
        model,
        Max,
        n / 2 * log(1 / (2 * π * σ^2)) -
        sum((data[i] - μ)^2 for i in 1:n) / (2 * σ^2)
    )
    optimize!(model)
    if verbose
        println("μ             = ", value(μ))
        println("mean(data)    = ", Statistics.mean(data))
        println("σ^2           = ", value(σ)^2)
        println("var(data)     = ", Statistics.var(data))
        println("MLE objective = ", objective_value(model))
    end
    Test.@test value(μ) ≈ Statistics.mean(data) atol = 1e-3
    Test.@test value(σ)^2 ≈ Statistics.var(data) atol = 1e-2
    ## You can even do constrained MLE!
    @NLconstraint(model, μ == σ^2)
    optimize!(model)
    Test.@test value(μ) ≈ value(σ)^2
    if verbose
        println()
        println("With constraint μ == σ^2:")
        println("μ                         = ", value(μ))
        println("σ^2                       = ", value(σ)^2)
        println("Constrained MLE objective = ", objective_value(model))
    end
    return
end

example_mle()

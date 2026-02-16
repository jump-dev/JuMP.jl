# Copyright 2017, Iain Dunning, Joey Huchette, Miles Lubin, and contributors    #src
# This Source Code Form is subject to the terms of the Mozilla Public License   #src
# v.2.0. If a copy of the MPL was not distributed with this file, You can       #src
# obtain one at https://mozilla.org/MPL/2.0/.                                   #src

# # Progressive Hedging

# The purpose of this tutorial is to demonstrate the Progressive Hedging
# algorithm. It may be helpful to read [Two-stage stochastic programs](@ref)
# first.

# ## Required packages

# This tutorial requires the following packages:

using JuMP
import Distributions
import Ipopt
import Printf

# ## Background

# Progressive Hedging (PH) is a popular decomposition algorithm for stochastic
# programming. It decomposes a stochastic problem into scenario subproblems that
# are solved iteratively, with penalty terms driving solutions toward consensus.

# In Progressive Hedging, each scenario subproblem includes a quadratic penalty
# term:
# ```math
# \min\limits_{x_s}: f_s(x_s) + \frac{\rho}{2} * ||x_s - \bar{x}||^2 + w_s^\top x_s
# ```
# where:
# - $x_s$ is the primal variable in scenario $s$
# - $f_s(x)$ is the original scenario objective
# - $\rho$ is the penalty parameter
# - $\bar{x}$ is the current consensus (average) solution
# - $w_s$ is the dual price (Lagrangian multiplier) in scenario $s$

# Progressive Hedging is an iterative algorithm. In each iteration, it solves
# all the penalized scenario subproblems, then it applies two updates:
#
# 1. $\bar{x} = \mathbb{E}_s[x_s]$
# 2. $w_s \pluseq \rho (x_s - \bar{x})$
#
# The algorithm terminates if $\bar{x} \approxeq x_s$ for all scenarios (the
# primal residual), and $\bar{x}$ has not changed by much between iterations
# (the dual residual).

# $\rho$ can be optionally updated between iterations. How to do so is an open
# question. There is a large literature on different updates strategies.

# In this tutorial we use parameters for $\rho$, $w$, and $\bar{x}$ to
# efficiently modify each scenario's subproblem between PH iterations.

# ## Building a single scenario

# The building block of Progressive Hedging is a separate JuMP model for each
# scenario. Here's an example, using the problem from
# [Two-stage stochastic programs](@ref):

function build_subproblem(; demand::Float64)
    model = Model(Ipopt.Optimizer)
    set_silent(model)
    @variable(model, x >= 0)
    @variable(model, 0 <= y <= demand)
    @constraint(model, y <= x)
    @variable(model, ρ in Parameter(1))
    @variable(model, x̄ in Parameter(0))
    @variable(model, w in Parameter(0))
    @expression(model, f_s, 2 * x - 5 * y + 0.1 * (x - y))
    @objective(model, Min, f_s + ρ / 2 * (x - x̄)^2 + w * x)
    return model
end

# Using the `build_subproblem` function, we can create one JuMP model for each
# scenario:

N = 10
demands = rand(Distributions.TriangularDist(150.0, 250.0, 200.0), N);
subproblems = map(demands) do demand
    return (; model = build_subproblem(; demand), probability = 1 / N)
end;

# ## The Progressive Hedging loop

# We're almost ready for our optimization loop, but first, here's a helpful
# function for logging:

function print_iteration(iter, args...)
    if mod(iter, 10) == 0
        f(x) = Printf.@sprintf("%15.4e", x)
        println(lpad(iter, 9), " ", join(f.(args), " "))
    end
    return
end

# Now we can implement our algorithm:

function solve_progressive_hedging(
    subproblems;
    iteration_limit::Int = 400,
    atol::Float64 = 1e-4,
    ρ::Float64 = 1.0,
)
    x̄_old, x̄ = 0.0, 0.0
    x, w = zeros(length(subproblems)), zeros(length(subproblems))
    println("iteration primal_residual   dual_residual")
    ## For each iteration...
    for iter in 1:iteration_limit
        ## For each subproblem...
        for (i, data) in enumerate(subproblems)
            ## Update the parameters
            set_parameter_value(data.model[:ρ], ρ)
            set_parameter_value(data.model[:x̄], x̄)
            set_parameter_value(data.model[:w], w[i])
            ## Solve the subproblem
            optimize!(data.model)
            assert_is_solved_and_feasible(data.model)
            ## Store the primal solution
            x[i] = value(data.model[:x])
        end
        ## Compute the consensus solution for the first-stage variables
        x̄ = sum(s.probability * x_s for (s, x_s) in zip(subproblems, x))
        ## Compute the primal and dual residuals
        primal_residual = maximum(abs, x_s - x̄ for x_s in x)
        dual_residual = ρ * abs(x̄ - x̄_old)
        print_iteration(iter, primal_residual, dual_residual)
        ## Check for convergence
        if primal_residual < atol && dual_residual < atol
            break
        end
        ## Update
        x̄_old = x̄
        w .+= ρ .* (x .- x̄)
    end
    return x̄
end

x̄ = solve_progressive_hedging(subproblems)

# ## Progressive Hedging with an adaptive penalty parameter

# You can also make the penalty parameter $\rho$ adaptive. How to do so is an
# open question. There is a large literature on different updates strategies.
# One approach is to increase $\rho$ if the primal residual is much larger than
# the dual residual, and to decrease $\rho$ if the dual residual is much larger
# than the primal residual.

function solve_adaptive_progressive_hedging(
    subproblems;
    iteration_limit::Int = 400,
    atol::Float64 = 1e-4,
    ρ::Float64 = 1.0,
    τ::Float64 = 1.3,
    μ::Float64 = 15.0,
)
    x̄_old, x̄ = 0.0, 0.0
    x, w = zeros(length(subproblems)), zeros(length(subproblems))
    println("iteration primal_residual   dual_residual")
    for iter in 1:iteration_limit
        for (i, data) in enumerate(subproblems)
            set_parameter_value(data.model[:ρ], ρ)
            set_parameter_value(data.model[:x̄], x̄)
            set_parameter_value(data.model[:w], w[i])
            optimize!(data.model)
            assert_is_solved_and_feasible(data.model)
            x[i] = value(data.model[:x])
        end
        x̄ = sum(s.probability * x_s for (s, x_s) in zip(subproblems, x))
        primal_residual = maximum(abs, x_s - x̄ for x_s in x)
        dual_residual = ρ * abs(x̄ - x̄_old)
        print_iteration(iter, primal_residual, dual_residual)
        if primal_residual < atol && dual_residual < atol
            break
        end
        w .+= ρ .* (x .- x̄)
        x̄_old = x̄
        ## Adaptive ρ update
        if primal_residual > μ * dual_residual
            ρ *= τ
        elseif dual_residual > μ * primal_residual
            ρ /= τ
        end
    end
    return x̄
end

x̄ = solve_adaptive_progressive_hedging(subproblems)

# Try tuning the values of `τ` and `μ`. Can you get the algorithm to converge
# in fewer iterations?

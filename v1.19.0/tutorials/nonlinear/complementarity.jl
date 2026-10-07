# Copyright 2017, Iain Dunning, Joey Huchette, Miles Lubin, and contributors    #src
# This Source Code Form is subject to the terms of the Mozilla Public License   #src
# v.2.0. If a copy of the MPL was not distributed with this file, You can       #src
# obtain one at https://mozilla.org/MPL/2.0/.                                   #src

# # Mixed complementarity problems

# This tutorial is a collection of mixed complementarity programs.

# This tutorial uses the following packages:

using JuMP
import PATHSolver
import Test  #src

# ## Background

# A mixed complementarity problem has the form:
# ```math
# \begin{align}
#     F_i(x) \perp x_i & i = 1 \ldots n \\
#     l_i \le x_i \le u_i & i = 1 \ldots n.
# \end{align}
# ```
# where the ``\perp`` constraint enforces the following relations:
#
#  - If ``l_i < x_i < u_i``, then ``F_i(x) = 0``
#  - If ``l_i = x_i``, then ``F_i(x) \ge 0``
#  - If ``x_i = u_i``, then ``F_i(x) \le 0``

# You may have seen a complementarity problem written as
# ``0 \le F(x) \perp x \ge 0``. This is a special case of a mixed
# complementarity problem in which ``l_i = 0`` and ``u_i = \infty``.

# Importantly, a mixed complementarity problem does not have an objective, and
# no other constraint types are present.

# ## Linear complementarity

# Form a mixed complementarity problem using the perp symbol `⟂` (type
# `\perp<tab>` in the REPL).

M = [0 0 -1 -1; 0 0 1 -2; 1 -1 2 -2; 1 2 -2 4]
q = [2, 2, -2, -6]
model = Model(PATHSolver.Optimizer)
set_silent(model)
@variable(model, 0 <= x[1:4] <= 10, start = 0)
@constraint(model, M * x + q ⟂ x)
optimize!(model)
Test.@test value.(x) ≈ [2.8, 0.0, 0.8, 1.2]  #src
value.(x)

# ## Other ways of writing linear complementarity problems

# You do not need to use a single vector of variables, and the complementarity
# constraints can be given in any order. In addition, you can use the perp
# symbol, the `complements(F, x)` syntax, or the [`MOI.Complements`](@ref) set.

model = Model(PATHSolver.Optimizer)
set_silent(model)
@variable(model, 0 <= w <= 10, start = 0)
@variable(model, 0 <= x <= 10, start = 0)
@variable(model, 0 <= y <= 10, start = 0)
@variable(model, 0 <= z <= 10, start = 0)
@constraint(model, complements(y - 2z + 2, x))
@constraint(model, [-y - z + 2, w] in MOI.Complements(2))
@constraint(model, w + 2x - 2y + 4z - 6 ⟂ z)
@constraint(model, w - x + 2y - 2z - 2 ⟂ y)
optimize!(model)
Test.@test value.([w, x, y, z]) ≈ [2.8, 0.0, 0.8, 1.2]  #src
value.([w, x, y, z])

# ## Transportation

# This example is a reformulation of the transportation problem from Chapter
# 3.3 of Dantzig, G.B. (1963). _Linear Programming and Extensions_. Princeton
# University Press, Princeton, New Jersey. It is based on the GAMS model
# [`gamslib_transmcp`](https://www.gams.com/latest/gamslib_ml/libhtml/gamslib_transmcp.html).

capacity = Dict("seattle" => 350, "san-diego" => 600)
demand = Dict("new-york" => 325, "chicago" => 300, "topeka" => 275)
cost = Dict(
    ("seattle" => "new-york") => 90 * 2.5 / 1_000,
    ("seattle" => "chicago") => 90 * 1.7 / 1_000,
    ("seattle" => "topeka") => 90 * 1.8 / 1_000,
    ("san-diego" => "new-york") => 90 * 2.5 / 1_000,
    ("san-diego" => "chicago") => 90 * 1.8 / 1_000,
    ("san-diego" => "topeka") => 90 * 1.4 / 1_000,
)
plants, markets = keys(capacity), keys(demand)
model = Model(PATHSolver.Optimizer)
set_silent(model)
@variable(model, w[i in plants] >= 0)
@variable(model, p[j in markets] >= 0)
@variable(model, x[i in plants, j in markets] >= 0)
@constraints(
    model,
    begin
        [i in plants, j in markets], w[i] + cost[i=>j] - p[j] ⟂ x[i, j]
        [i in plants], capacity[i] - sum(x[i, :]) ⟂ w[i]
        [j in markets], sum(x[:, j]) - demand[j] ⟂ p[j]
    end
)
optimize!(model)
Test.@test isapprox(value(p["new-york"]), 0.225; atol = 1e-3)  #src
value.(p)

# ## Expected utility of insurance

# This example is taken from a lecture of the course AAE706, given by Thomas F.
# Rutherford at the University of Wisconsin, Madison. It models the expected
# coverage of insurance `K` that a rational actor would obtain to insure a risk
# that occurs with probability `pi` and results in a loss of `L`.

pi = 0.01  # Probability of a bad outcome
L = 0.5    # Loss with a bad outcome
γ = 0.02   # Premium for coverage
σ = 0.5    # Elasticity
ρ = -1     # Risk exponent
U(C) = C^ρ / ρ
MU(C) = C^(ρ - 1)
model = Model(PATHSolver.Optimizer)
set_silent(model)
@variable(model, EU, start = 1)   # Expected utilitiy
@variable(model, EV, start = 1)   # Equivalent variation in income
@variable(model, C_G, start = 1)  # Consumption on a good day
@variable(model, C_B, start = 1)  # Consumption on a bad day
@variable(model, K, start = 1)    # Coverage
@constraints(
    model,
    begin
        (1 - pi) * U(C_G) + pi * U(C_B) - EU ⟂ EU
        100 * (((1 - pi) * C_G^ρ + pi * C_B^ρ)^(1 / ρ) - 1) - EV ⟂ EV
        1 - γ * K - C_G ⟂ C_G
        1 - L + (1 - γ) * K - C_B ⟂ C_B
        γ * ((1 - pi) * MU(C_G) + pi * MU(C_B)) - pi * MU(C_B) ⟂ K
    end
)
optimize!(model)
Test.@test isapprox(value(C_G), 0.996; atol = 1e-3)  #src
value(K)

# ## Electricity consumption

# This example is a mixed complementarity formulation of Example 3.3.1 from
# [DAertrycke2017](@cite).

# This example models a risk neutral competitive equilibrium between a producer
# and a consumer of electricity.

# In our example, we assume a producer is looking to invest in a new power
# plant with capacity ``x`` [MW]. This plant has an annualized capital cost of
# ``I`` [€/MW] and an operating cost of ``C`` [€/MWh]. There are 8760 hours in a
# year.
#
# After making the capital investment, there are five possible consumption
# scenarios, ``\omega``, which occur with probability ``\theta_\omega``. In each
# scenario , the producer makes ``Y_ω`` MW of electricity.
#
# There is one consumer in the model, who has a quadratic utility function,
# ``U(Q_ω) = A_ω Q_ω + \frac{B_ω Q_ω^2}{2}``.
#
# We now build and solve the mixed complementarity problem with a few brief
# comments. The economic justification for the model would require a larger
# tutorial than the space available here. Consult [DAertrycke2017](@cite) for
# details.

I = 90_000                     # Annualized capital cost
C = 60                         # Operation cost per MWh
τ = 8_760                      # Hours per year
θ = [0.2, 0.2, 0.2, 0.2, 0.2]  # Scenario probabilities
A = [300, 350, 400, 450, 500]  # Utility function coefficients
B = 1                          # Utility function coefficients
model = Model(PATHSolver.Optimizer)
set_silent(model)
@variable(model, x >= 0, start = 1)           # Installed capacity
@variable(model, Q[ω = 1:5] >= 0, start = 1)  # Consumption
@variable(model, Y[ω = 1:5] >= 0, start = 1)  # Production
@variable(model, P[ω = 1:5], start = 1)       # Electricity price
@variable(model, μ[ω = 1:5] >= 0, start = 1)  # Capital scarcity margin
## Unit investment cost equals annualized scarcity profit or investment is 0
@constraint(model, I - τ * θ' * μ ⟂ x)
## Difference between price and scarcity margin is equal to operation cost
@constraint(model, [ω = 1:5], C - (P[ω] - μ[ω]) ⟂ Y[ω])
## Price is equal to consumer's marginal utility
@constraint(model, [ω = 1:5], P[ω] - (A[ω] - B * Q[ω]) ⟂ Q[ω])
## Production is equal to consumption
@constraint(model, [ω = 1:5], Y[ω] - Q[ω] ⟂ P[ω])
## Production does not exceed capacity
@constraint(model, [ω = 1:5], x - Y[ω] ⟂ μ[ω])
optimize!(model)
solution_summary(model)

# An equilibrium solution is to build 389 MW:

Test.@test isapprox(value(x), 389; atol = 1)  #src
value(x)

# The production in each scenario is:

Test.@test isapprox(value.(Q), [240, 290, 340, 389, 389]; atol = 1)  #src
value.(Q)

# The price in each scenario is:

Test.@test isapprox(value.(P), [60, 60, 60, 61, 111]; atol = 1)  #src
value.(P)

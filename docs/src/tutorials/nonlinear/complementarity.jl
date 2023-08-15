# Copyright 2017, Iain Dunning, Joey Huchette, Miles Lubin, and contributors    #src
# This Source Code Form is subject to the terms of the Mozilla Public License   #src
# v.2.0. If a copy of the MPL was not distributed with this file, You can       #src
# obtain one at https://mozilla.org/MPL/2.0/.                                   #src

# # Mixed complementarity problems

# This tutorial is a collection of examples of small mixed-complementarity
# programs. See [Complementarity constraints](@ref) for the definition of a
# complementarity constraint.

# This tutorial uses the following packages:

using JuMP
import PATHSolver
import Test  #src

# ## Vector-valued linear complementarity

M = [0 0 -1 -1; 0 0 1 -2; 1 -1 2 -2; 1 2 -2 4]
q = [2, 2, -2, -6]
model = Model(PATHSolver.Optimizer)
@variable(model, 0 <= x[1:4] <= 10, start = 0)
@constraint(model, M * x + q ⟂ x)
optimize!(model)
Test.@test value.(x) ≈ [2.8, 0.0, 0.8, 1.2]  #src
value.(x)

# ## Scalar-valued linear complementarity

# You do not need to use a single vector of variables, and the complementarity
# constraints can be given in any order. In addition, you can either use the
#  perp symbol `⟂` (type `\perp<tab>` in the REPL), or you can use the
# [`MOI.Complements`](@ref) set.

model = Model(PATHSolver.Optimizer)
@variable(model, 0 <= w <= 10, start = 0)
@variable(model, 0 <= x <= 10, start = 0)
@variable(model, 0 <= y <= 10, start = 0)
@variable(model, 0 <= z <= 10, start = 0)
@constraint(model, [y - 2z + 2, x] in MOI.Complements(2))
@constraint(model, -y - z + 2 ⟂ w)
@constraint(model, w + 2x - 2y + 4z - 6 ⟂ z)
@constraint(model, w - x + 2y - 2z - 2 ⟂ y)
optimize!(model)
Test.@test value.([w, x, y, z]) ≈ [2.8, 0.0, 0.8, 1.2]  #src
value.([w, x, y, z])

# ## Transportation

# This is example is a reformulation of the transportation problem from Chapter
# 3.3 of Dantzig, G.B. (1963). _Linear Programming and Extensions_. Princeton
# University Press, Princeton, New Jersey. It is based on the GAMS model
# [gamslib_transmcp](https://www.gams.com/latest/gamslib_ml/libhtml/gamslib_transmcp.html).

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

# ## Expected utility of insurannce

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
@variable(model, EU, start = 1)   # Expected utilitiy
@variable(model, EV, start = 1)   # Equivalent variation in income
@variable(model, C_G, start = 1)  # Consumption on a good day
@variable(model, C_B, start = 1)  # Consumptio on a bad day
@variable(model, K, start = 1)    # Coverage
@constraints(
    model,
    begin
        EU - ((1 - pi) * U(C_G) + pi * U(C_B)) ⟂ EU
        EV - 100 * (((1 - pi) * C_G^ρ + pi * C_B^ρ)^(1 / ρ) - 1) ⟂ EV
        C_G - (1 - γ * K) ⟂ C_G
        C_B - (1 - L + (1 - γ) * K) ⟂ C_B
        γ * ((1 - pi) * MU(C_G) + pi * MU(C_B)) - pi * MU(C_B) ⟂ K
    end
)
optimize!(model)
Test.@test isapprox(value(C_G), 0.996; atol = 1e-3)  #src
value(K)

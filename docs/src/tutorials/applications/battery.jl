# Copyright 2017, Iain Dunning, Joey Huchette, Miles Lubin, and contributors    #src
# This Source Code Form is subject to the terms of the Mozilla Public License   #src
# v.2.0. If a copy of the MPL was not distributed with this file, You can       #src
# obtain one at https://mozilla.org/MPL/2.0/.                                   #src

# # Battery scheduling with complementarity constraints

# This tutorial schedules a battery to profit from changes in electricity
# prices. We write one model with complementarity constraints, then solve it
# using HiGHS, Ipopt, and CCOpt.
#
# **Learning intentions:**
# * Prevent simultaneous charging and discharging using complementarity.
# * Bridge complementarity to a mixed-integer formulation for HiGHS.
# * Select a Fischer–Burmeister reformulation using a bridge attribute.
# * Pass complementarity constraints to the specialized solver CCOpt.

# ## Required packages

# Install the packages with `import Pkg; Pkg.add(["JuMP", "HiGHS", "Ipopt",
# "MathOptComplements", "NLPModelsJuMP", "CCOpt"])`.
# Use MathOptComplements v0.1.2 or later for the bridge attributes below.

using JuMP
import CCOpt
import HiGHS
import Ipopt
import MathOptComplements
import NLPModelsJuMP
import Test

# Importing NLPModelsJuMP and MathOptComplements activates CCOpt's JuMP
# interface.

# ## The battery problem

# Our model is a deterministic simplification of the storage model in Section 3
# of [Chancelier2026](@cite).
# Their model includes uncertain prices, bid–ask spreads, and risk preferences.
# Here, we assume known prices and maximize trading revenue over six one-hour
# periods. The data are illustrative, rather than observations from a market.

prices = [-30.0, -10.0, 20.0, 80.0, 30.0, 60.0]  # EUR/MWh
capacity = 1.0          # MWh
power_limit = 0.6       # MW
duration = 1.0          # Hours per period
efficiency = 0.9        # Charging and discharging efficiency
initial_energy = 0.5    # MWh

# Let `charge[t]` and `discharge[t]` be power drawn from and delivered to the
# grid, respectively. The stored energy satisfies
#
# `energy[t] = energy[t-1] + duration *
#              (efficiency * charge[t] - discharge[t] / efficiency)`.
#
# We bound the stored energy and power, and require the battery to finish with
# its initial energy. This prevents earning money by depleting the initial
# stock. Revenue is `sum(duration * prices[t] * (discharge[t] - charge[t]))`.
# Negative prices mean that the operator is paid to consume electricity.
#
# The constraint `0 ≤ charge[t] ⟂ discharge[t] ≥ 0` prevents charging and
# discharging at the same time: both quantities are nonnegative and their
# product is zero.

function build_battery_model(optimizer)
    model = Model(optimizer)
    MathOptComplements.Bridges.add_all_bridges(model)
    T = length(prices)
    @variable(model, 0 <= energy[0:T] <= capacity, start = initial_energy)
    @variable(model, 0 <= charge[1:T] <= power_limit, start = 0.1)
    @variable(model, 0 <= discharge[1:T] <= power_limit, start = 0.1)
    @constraint(model, energy[0] == initial_energy)
    @constraint(model, energy[T] == initial_energy)
    @constraint(
        model,
        [t in 1:T],
        energy[t] ==
        energy[t-1] +
        duration * (efficiency * charge[t] - discharge[t] / efficiency),
    )
    @constraint(
        model,
        exclusive[t in 1:T],
        [charge[t], discharge[t]] in
        MathOptComplements.ComplementsWithSetType{MOI.Nonnegatives}(2),
    )
    @objective(
        model,
        Max,
        sum(duration * prices[t] * (discharge[t] - charge[t]) for t in 1:T),
    )
    return model
end

# We explicitly select nonnegative complementarity with
# `ComplementsWithSetType{MOI.Nonnegatives}`. Writing
# `charge[t] ⟂ discharge[t]` instead would use the variable bounds on
# `discharge[t]` to define a *mixed* complementarity condition. Its finite upper
# bound has a different meaning from the power limit we want here.

# ## HiGHS: a mixed-integer reformulation

# MathOptComplements bridges each nonnegative complementarity pair to an SOS1
# constraint. An SOS1 constraint allows at most one entry to be nonzero.
# HiGHS does not support SOS1 directly, so MathOptInterface bridges SOS1 to
# binary variables and linear constraints using the finite power bounds.
# Conceptually, the resulting formulation uses a binary mode `z[t]`:
#
# * `charge[t] ≤ power_limit * z[t]`
# * `discharge[t] ≤ power_limit * (1 - z[t])`
#
# There is no need to write these binary variables ourselves.

highs_model = build_battery_model(HiGHS.Optimizer)
set_silent(highs_model)
optimize!(highs_model)
assert_is_solved_and_feasible(highs_model)
Test.@test termination_status(highs_model) == MOI.OPTIMAL
objective_value(highs_model)

# Inspect the active bridges to see the SOS1-to-MIP conversion:

print_active_bridges(highs_model)

# Inspect the schedule. Positive net power means discharging to the grid.

[
    (
        period = t,
        price = prices[t],
        charge = value(highs_model[:charge][t]),
        discharge = value(highs_model[:discharge][t]),
        energy = value(highs_model[:energy][t]),
    ) for t in eachindex(prices)
]

# ## Ipopt: choose a reformulation with a bridge attribute

# Ipopt requires smooth nonlinear constraints. Wrap it in
# `MathOptComplements.Optimizer` to configure the complementarity bridge using
# the `DefaultComplementarityReformulation()` optimizer attribute.

ipopt_model =
    build_battery_model(() -> MathOptComplements.Optimizer(Ipopt.Optimizer()))
set_attribute(
    ipopt_model,
    MathOptComplements.DefaultComplementarityReformulation(),
    MathOptComplements.FischerBurmeisterRelaxation(1e-8),
)
set_attribute(ipopt_model, "mu_strategy", "adaptive")
set_attribute(ipopt_model, "bound_relax_factor", 0.0)
set_attribute(ipopt_model, "tol", 1e-9)
set_silent(ipopt_model)
optimize!(ipopt_model)
assert_is_solved_and_feasible(ipopt_model)
objective_value(ipopt_model)

# For nonnegative `a` and `b`, the bridge imposes
# `a + b - sqrt(a^2 + b^2 + ε) ≤ 0`.
# With `ε = 0`, this is equivalent to `a * b = 0`, but it is nonsmooth at the
# origin. A positive `ε` smooths the expression and relaxes complementarity to
# `a * b ≤ ε / 2`. Consequently, we check the original complementarity residual
# below, rather than relying only on the solver's termination status.
#
# Ipopt finds a local solution of this nonlinear relaxation. It does not provide
# the global optimality certificate that HiGHS provides for the MIP.

# ## CCOpt: solve the complementarity problem directly

# CCOpt handles mathematical programs with complementarity constraints using
# specialized relaxation and penalty algorithms built on MadNLP. See
# [Pozharskiy2026](@cite) for the algorithms and implementation.
#
# Its JuMP interface accepts our nonnegative complementarity pairs directly.
# Do not wrap CCOpt in `MathOptComplements.Optimizer`: we want CCOpt to receive
# the complementarity structure, rather than a generic nonlinear reformulation.
# We request a stopping tolerance of `1e-6` and check the physical residuals
# against the same tolerance below.

ccopt_model = build_battery_model(CCOpt.Optimizer)
set_attribute(ccopt_model, "tol", 1e-6)
set_silent(ccopt_model)
optimize!(ccopt_model)
assert_is_solved_and_feasible(ccopt_model)
objective_value(ccopt_model)

# CCOpt is also a local solver. The HiGHS optimum provides a reference for
# assessing the revenues obtained by both local solvers.

# ## Check the solutions

# Compare revenue, status, and the maximum product of charging and discharging
# power. The product should be zero, up to the chosen numerical tolerance.

function battery_result(model)
    residual = maximum(value.(model[:charge]) .* value.(model[:discharge]))
    return (
        revenue = objective_value(model),
        complementarity_residual = residual,
        status = termination_status(model),
    )
end

results = (
    HiGHS = battery_result(highs_model),
    Ipopt = battery_result(ipopt_model),
    CCOpt = battery_result(ccopt_model),
)

# Check physical feasibility in addition to the bridged model's solver status.

for model in (highs_model, ipopt_model, ccopt_model)
    e = value.(model[:energy])
    c = value.(model[:charge])
    d = value.(model[:discharge])
    Test.@test all(-1e-6 <= e[t] <= capacity + 1e-6 for t in 0:length(prices))
    Test.@test all(-1e-6 <= x <= power_limit + 1e-6 for x in [c; d])
    Test.@test isapprox(e[0], initial_energy; atol = 1e-6)
    Test.@test isapprox(e[length(prices)], initial_energy; atol = 1e-6)
    Test.@test all(
        isapprox(
            e[t],
            e[t-1] + duration * (efficiency * c[t] - d[t] / efficiency);
            atol = 1e-6,
        ) for t in eachindex(prices)
    )
    Test.@test maximum(abs.(c .* d)) <= 1e-6
end

# The local solutions cannot improve on the global MIP optimum, apart from
# small numerical differences allowed by the smoothing and solver tolerances.

Test.@test results.Ipopt.revenue <= results.HiGHS.revenue + 1e-3
Test.@test results.CCOpt.revenue <= results.HiGHS.revenue + 1e-3
results

# ## When can we omit complementarity?

# Proposition 2 of [Chancelier2026](@cite) gives a sufficient condition for their
# convex relaxation to be exact: `s⁻ / c⁻ ≤ s⁺ / c⁺`, where `s⁻` and `s⁺` are
# bid and ask prices, and `c⁻` and `c⁺` are discharge and charge coefficients.
# In our model, these are `s⁻ = s⁺ = price`, `c⁺ = efficiency`, and
# `c⁻ = 1 / efficiency`. Negative prices violate this condition. Without
# complementarity, simultaneous charging and discharging could consume energy
# through losses and collect payments at negative prices.
#
# When the condition holds, their construction recovers a complementary policy
# from the convex relaxation without sacrificing revenue. For the stochastic
# problem, preserving convexity allows them to use SDDP.

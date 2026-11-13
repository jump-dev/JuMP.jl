# Copyright 2017, Iain Dunning, Joey Huchette, Miles Lubin, and contributors    #src
# This Source Code Form is subject to the terms of the Mozilla Public License   #src
# v.2.0. If a copy of the MPL was not distributed with this file, You can       #src
# obtain one at https://mozilla.org/MPL/2.0/.                                   #src

# # Two-stage stochastic programs

# The purpose of this tutorial is to demonstrate how to model and solve a
# two-stage stochastic program.

# !!! info
#     The JuMP extension [InfiniteOpt.jl](../../packages/InfiniteOpt.md) can also be
#     used to model and solve two-stage stochastic programs.
#     The JuMP extension [SDDP.jl](../../packages/SDDP.md) can be
#     used to model and solve multi-stage stochastic programs.

# This tutorial uses the following packages

using JuMP
import Distributions
import HiGHS
import Plots
import StatsPlots
import Statistics

# ## Background

# During the week, you are a busy practitioner of Operations Research. To escape
# the drudgery of mathematics, you decide to open a side business selling creamy
# mushroom pies with puff pastry. After a few weeks, it quickly becomes apparent
# that operating a food business is not so easy.

# The pies must be prepared in the morning, _before_ you open for the day and
# can gauge the level of demand. If you bake too many, the unsold pies at the
# end of the day must be discarded and you have wasted time and money on their
# production. But if you bake too few, then there may be un-served customers and
# you could have made more money by baking more pies.

# After a few weeks of poor decision making, you decide to put your knowledge of
# Operations Research to good use, starting with some data collection.

# Each pie costs you \$2 to make, and you sell them at \$5 each. Disposal of an
# unsold pie costs \$0.10. Based on three weeks of data collected, in which you
# made 200 pies each week, you sold 150, 190, and 200 pies. Thus, as a guess,
# you assume a triangular distribution of demand with a minimum of 150, a median
# of 200, and a maximum of 250.

# We can model this problem by a two-stage stochastic program. In the first
# stage, we decide a quantity of pies to make ``x``. We make this decision
# before we observe the demand ``d_\omega``. In the second stage, we sell
# ``y_\omega`` pies, and incur any costs for unsold pies.

# We can formulate this problem as follows:
# ```math
# \begin{aligned}
# \max\limits_{x,y_\omega} \;\; & -2x + \mathbb{E}_\omega[5y_\omega - 0.1(x - y_\omega)] \\
#   & y_\omega \le x              & \quad \forall \omega \in \Omega \\
#   & 0 \le y_\omega \le d_\omega & \quad \forall \omega \in \Omega \\
#   & x \ge 0.
# \end{aligned}
# ```

# ## Sample Average approximation

# If the distribution of demand is continuous, then our problem has an infinite
# number of variables and constraints. To form a computationally tractable
# problem, we instead use a finite set of samples drawn from the distribution.
# This is called sample average approximation (SAA).

D = Distributions.TriangularDist(150.0, 250.0, 200.0)
N = 100
d = sort!(rand(D, N));
Ω = 1:N
P = fill(1 / N, N);
StatsPlots.histogram(d; bins = 20, label = "", xlabel = "Demand")

# ## JuMP model

# The implementation of our two-stage stochastic program in JuMP is:

model = Model(HiGHS.Optimizer)
set_silent(model)
@variable(model, x >= 0)
@variable(model, 0 <= y[ω in Ω] <= d[ω])
@constraint(model, [ω in Ω], y[ω] <= x)
@expression(model, z[ω in Ω], 5y[ω] - 0.1 * (x - y[ω]))
@objective(model, Max, -2x + sum(P[ω] * z[ω] for ω in Ω))
optimize!(model)
@assert is_solved_and_feasible(model)
solution_summary(model)

# The optimal number of pies to make is:

value(x)

# The distribution of total profit is:

total_profit = [-2 * value(x) + value(z[ω]) for ω in Ω]

# Let's plot it:

"""
    bin_distribution(x::Vector{Float64}, N::Int)

A helper function that discretizes `x` into bins of width `N`.
"""
bin_distribution(x, N) = N * (floor(minimum(x) / N):ceil(maximum(x) / N))

plot = StatsPlots.histogram(
    total_profit;
    bins = bin_distribution(total_profit, 25),
    label = "",
    xlabel = "Profit [\$]",
    ylabel = "Number of outcomes",
)
μ = Statistics.mean(total_profit)
Plots.vline!(
    plot,
    [μ];
    label = "Expected profit (\$$(round(Int, μ)))",
    linewidth = 3,
)
plot

# ## Risk measures

# A risk measure is a function which maps a random variable to a real number.
# Common risk measures include the mean (expectation), median, mode, and
# maximum. We need a risk measure to convert the distribution of second stage
# costs into a single number that can be optimized.

# Our model currently uses the expectation risk measure, but others are possible
# too. One popular risk measure is the conditional value at risk (CVaR).

# CVaR has a parameter $\gamma$, and it computes the expectation of the worst
# $\gamma$ fraction of outcomes.

# If we are maximizing, so that small outcomes are bad, the definition of CVaR
# is:
# ```math
# CVaR_{\gamma}[Z] = \max\limits_{\xi} \;\; \xi - \frac{1}{\gamma}\mathbb{E}_\omega\left[(\xi - Z)_+\right]
# ```
# which can be formulated as the linear program:
# ```math
# \begin{aligned}
# CVaR_{\gamma}[Z] = \max\limits_{\xi, z_\omega} \;\; & \xi - \frac{1}{\gamma}\sum P_\omega z_\omega\\
#  & z_\omega \ge \xi - Z_\omega & \quad \forall \omega \\
#  & z_\omega \ge 0 & \quad \forall \omega.
# \end{aligned}
# ```

function CVaR(Z::Vector{Float64}, P::Vector{Float64}; γ::Float64)
    @assert 0 < γ <= 1
    N = length(Z)
    model = Model(HiGHS.Optimizer)
    set_silent(model)
    @variable(model, ξ)
    @variable(model, z[1:N] >= 0)
    @constraint(model, [i in 1:N], z[i] >= ξ - Z[i])
    @objective(model, Max, ξ - 1 / γ * sum(P[i] * z[i] for i in 1:N))
    optimize!(model)
    @assert is_solved_and_feasible(model)
    return objective_value(model)
end

# When `γ` is `1.0`, we compute the mean of the profit:

cvar_10 = CVaR(total_profit, P; γ = 1.0)

#-

Statistics.mean(total_profit)

# As `γ` approaches `0.0`, we compute the worst-case (minimum) profit:

cvar_00 = CVaR(total_profit, P; γ = 0.0001)

#-

minimum(total_profit)

# By varying `γ` between `0` and `1` we can compute some trade-off of these two
# extremes:

cvar_05 = CVaR(total_profit, P; γ = 0.5)

# Let's plot these outcomes on our distribution:

plot = StatsPlots.histogram(
    total_profit;
    bins = bin_distribution(total_profit, 25),
    label = "",
    xlabel = "Profit [\$]",
    ylabel = "Number of outcomes",
)
Plots.vline!(
    plot,
    [cvar_10 cvar_05 cvar_00];
    label = ["γ = 1.0" "γ = 0.5" "γ = 0.0"],
    linewidth = 3,
)
plot

# ## Risk averse sample average approximation

# Because CVaR can be formulated as a linear program, we can form a risk averse
# sample average approximation model by combining the two formulations:

γ = 0.4
model = Model(HiGHS.Optimizer)
set_silent(model)
@variable(model, x >= 0)
@variable(model, 0 <= y[ω in Ω] <= d[ω])
@constraint(model, [ω in Ω], y[ω] <= x)
@expression(model, Z[ω in Ω], 5 * y[ω] - 0.1(x - y[ω]))
@variable(model, ξ)
@variable(model, z[ω in Ω] >= 0)
@constraint(model, [ω in Ω], z[ω] >= ξ - Z[ω])
@objective(model, Max, -2x + ξ - 1 / γ * sum(P[ω] * z[ω] for ω in Ω))
optimize!(model)
@assert is_solved_and_feasible(model)

# When ``\gamma = 0.4``, the optimal number of pies to bake is:

value(x)

# The distribution of total profit is:

risk_averse_total_profit = [value(-2x + Z[ω]) for ω in Ω]
bins = bin_distribution([total_profit; risk_averse_total_profit], 25)
plot = StatsPlots.histogram(total_profit; label = "Expectation", bins = bins)
StatsPlots.histogram!(
    plot,
    risk_averse_total_profit;
    label = "CV@R",
    bins = bins,
    alpha = 0.5,
)
plot

# ## Next steps
#
#  * Try solving this problem for different numbers of samples and different
#    distributions.
#  * Refactor the example to avoid hard-coding the costs. What happens to the
#    solution if the cost of disposing unsold pies increases?
#  * Plot the optimal number of pies to make for different values of the risk
#    aversion parameter $\gamma$. What is the relationship?

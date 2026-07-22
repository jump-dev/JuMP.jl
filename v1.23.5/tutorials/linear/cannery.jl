# Copyright 2017, Iain Dunning, Joey Huchette, Miles Lubin, and contributors    #src
# This Source Code Form is subject to the terms of the Mozilla Public License   #src
# v.2.0. If a copy of the MPL was not distributed with this file, You can       #src
# obtain one at https://mozilla.org/MPL/2.0/.                                   #src

# # The cannery problem

# **This tutorial was originally contributed by Louis Luangkesorn.**

# This tutorial solves the cannery problem from Dantzig,
# _Linear Programming and Extensions_, Princeton University Press, Princeton,
# NJ, 1963. This class of problem is known as a transshipment problem.

# The purpose of this tutorial is to demonstrate how to use JSON data in the
# formulation of a JuMP model.

# ## Required packages

# This tutorial requires the following packages:

using JuMP
import HiGHS
import JSON
import Test

# ## Formulation

# The cannery problem assumes we are optimizing the shipment of cases of
# cans from production plants ``p \in P`` to markets ``m \in M``.

# Each production plant ``p`` has a capacity ``c_p``, and each market ``m``
# has a demand ``d_m``. The shipping cost per case of cans from plant ``p``
# to market ``m`` is ``d_{p,m}``.

# We wish to find the distribution plan ``x_{p,m}``, the number of cases of cans
# to ship from plant ``p`` to market ``m``, for ``p \in P`` and ``m \in M``
# that minimizes the shipping costs. We can formulate our problem as the
# following linear program:
# ```math
# \begin{aligned}
# \min & \sum\limits_{p \in P}\sum\limits_{m \in M} d_{p,m} x_{p,m} \\
# \text{s.t.} & \sum\limits_{m \in M} x_{p,m} \le c_p, && \forall p \in P \\
#             & \sum\limits_{p \in P} x_{p,m} \ge d_m, && \forall m \in M \\
#             & x_{p,m} \ge 0, && \forall p \in P, m \in M
# \end{aligned}
# ```

# ## Data

# A key feature of the tutorial is to demonstrate how to load data from JSON.

# For simplicity, we've hard-coded it below. But if the data was available as a
# `.json` file, we could use `data = JSON.parsefile(filename)` to read in the
# data.

data = JSON.parse("""
{
    "plants": {
        "Seattle": {"capacity": 350},
        "San-Diego": {"capacity": 600}
    },
    "markets": {
        "New-York": {"demand": 300},
        "Chicago": {"demand": 300},
        "Topeka": {"demand": 300}
    },
    "distances": {
        "Seattle => New-York": 2.5,
        "Seattle => Chicago": 1.7,
        "Seattle => Topeka": 1.8,
        "San-Diego => New-York": 2.5,
        "San-Diego => Chicago": 1.8,
        "San-Diego => Topeka": 1.4
    }
}
""")

# Create the set of plants:

P = keys(data["plants"])

# Create the set of markets:

M = keys(data["markets"])

# We also need a function to compute the distance from plant to market:

distance(p::String, m::String) = data["distances"]["$(p) => $(m)"]

# ## JuMP formulation

# Now we're ready to convert our mathematical formulation into a JuMP model.

# First, create a new JuMP model. Since we have a linear program, we'll use
# HiGHS as our optimizer:

model = Model(HiGHS.Optimizer)

# Our decision variables are indexed over the set of plants and markets:

@variable(model, x[P, M] >= 0)

# We need a constraint that each plant can ship no more than its capacity:

@constraint(model, [p in P], sum(x[p, :]) <= data["plants"][p]["capacity"])

# and each market must receive at least its demand:

@constraint(model, [m in M], sum(x[:, m]) >= data["markets"][m]["demand"])

# Finally, our objective is to minimize the transportation distance:

@objective(model, Min, sum(distance(p, m) * x[p, m] for p in P, m in M));

# ## Solution

# Let's optimize and look at the solution:

optimize!(model)
solution_summary(model)

# What's the optimal shipment?

Test.@test is_solved_and_feasible(model)
Test.@test isapprox(objective_value(model), 1_680.0, atol = 1e-6)  #src
for p in P, m in M
    println(p, " => ", m, ": ", value(x[p, m]))
end

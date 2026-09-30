# Copyright 2017, Iain Dunning, Joey Huchette, Miles Lubin, and contributors    #src
# This Source Code Form is subject to the terms of the Mozilla Public License   #src
# v.2.0. If a copy of the MPL was not distributed with this file, You can       #src
# obtain one at https://mozilla.org/MPL/2.0/.                                   #src

# # The diet problem

# This tutorial solves the classic "diet problem", also known as the
# [Stigler diet](https://en.wikipedia.org/wiki/Stigler_diet).

# ## Required packages

using JuMP
import DataFrames
import HiGHS
import Test  #hide

# ## Formulation

# Suppose we wish to cook a nutritionally balanced meal by choosing the quantity
# of each food ``f`` to eat from a set of foods ``F`` in our kitchen.

# Each food ``f`` has a cost, ``c_f``, as well as a macronutrient profile
# ``a_{m,f}`` for each macronutrient ``m \in M``.

# Because we care about a nutritionally balanced meal, we set some minimum and
# maximum limits for each nutrient, which we denote ``l_m`` and ``u_m``
# respectively.

# Furthermore, because we are optimizers, we seek the minimum cost solution.

# With a little effort, we can formulate our dinner problem as the following
# linear program:
# ```math
# \begin{aligned}
# \min & \sum\limits_{f \in F} c_f x_f \\
# \text{s.t.}\ \ & l_m \le \sum\limits_{f \in F} a_{m,f} x_f \le u_m, && \forall m \in M \\
# & x_f \ge 0, && \forall f \in F
# \end{aligned}
# ```

# In the rest of this tutorial, we will create and solve this problem in JuMP,
# and learn what we should cook for dinner.

# ## Data

# First, we need some data for the problem:

foods = DataFrames.DataFrame(
    [
        "hamburger" 2.49 410 24 26 730
        "chicken" 2.89 420 32 10 1190
        "hot dog" 1.50 560 20 32 1800
        "fries" 1.89 380 4 19 270
        "macaroni" 2.09 320 12 10 930
        "pizza" 1.99 320 15 12 820
        "salad" 2.49 320 31 12 1230
        "milk" 0.89 100 8 2.5 125
        "ice cream" 1.59 330 8 10 180
    ],
    ["name", "cost", "calories", "protein", "fat", "sodium"],
)

# Here, ``F`` is `foods.name` and ``c_f`` is `foods.cost`. (We're also playing
# a bit loose the term "macronutrient" by including calories and sodium.)

# !!! tip
#     Although we hard-coded the data here, you could also read it in from a
#     file. See [Getting started with data and plotting](@ref) for details.

# We also need our minimum and maximum limits:

limits = DataFrames.DataFrame(
    [
        "calories" 1800 2200
        "protein" 91 Inf
        "fat" 0 65
        "sodium" 0 1779
    ],
    ["name", "min", "max"],
)

# ## JuMP formulation

# Now we're ready to convert our mathematical formulation into a JuMP model.

# First, create a new JuMP model. Since we have a linear program, we'll use
# HiGHS as our optimizer:

model = Model(HiGHS.Optimizer)

# Next, we create a set of decision variables `x`, indexed over the foods in the
# `data` DataFrame. Each `x` has a lower bound of `0`.

@variable(model, x[foods.name] >= 0)

# Our objective is to minimize the total cost of purchasing food. We can write
# that as a sum over the rows in `data`.

@objective(
    model,
    Min,
    sum(food["cost"] * x[food["name"]] for food in eachrow(foods)),
)

# For the next component, we need to add a constraint that our total intake of
# each component is within the limits contained in the `limits` DataFrame.
# To make this more readable, we introduce a JuMP `@expression`

for limit in eachrow(limits)
    intake = @expression(
        model,
        sum(food[limit["name"]] * x[food["name"]] for food in eachrow(foods)),
    )
    @constraint(model, limit.min <= intake <= limit.max)
end

# What does our model look like?

print(model)

# ## Solution

optimize!(model)

Test.@test primal_status(model) == FEASIBLE_POINT   #hide
Test.@test objective_value(model) ≈ 11.8288 atol = 1e-4 #hide

solution_summary(model)

# Success! We found an optimal solution. Let's see what the optimal solution is:

for food in foods.name
    println(food, " = ", value(x[food]))
end

# That's a lot of milk and ice cream! And sadly, we only get `0.6` of a
# hamburger.

# ## Problem modification

# JuMP makes it easy to take an existing model and modify it by adding extra
# constraints. Let's see what happens if we add a constraint that we can buy at
# most 6 units of milk or ice cream combined.

@constraint(model, x["milk"] + x["ice cream"] <= 6)

#-

optimize!(model)

Test.@test termination_status(model) == INFEASIBLE  #hide
Test.@test primal_status(model) == NO_SOLUTION      #hide

solution_summary(model)

# Uh oh! There exists no feasible solution to our problem. Looks like we're stuck
# eating ice cream for dinner.

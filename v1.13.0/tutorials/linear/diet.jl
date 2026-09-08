# Copyright 2017, Iain Dunning, Joey Huchette, Miles Lubin, and contributors    #src
# This Source Code Form is subject to the terms of the Mozilla Public License   #src
# v.2.0. If a copy of the MPL was not distributed with this file, You can       #src
# obtain one at https://mozilla.org/MPL/2.0/.                                   #src

# # The diet problem

# This tutorial solves the classic "diet problem," also known as the
# [Stigler diet](https://en.wikipedia.org/wiki/Stigler_diet).

# ## Required packages

# This tutorial requires the following packages:

using JuMP
import DataFrames
import HiGHS
import Test  #hide

# ## Formulation

# Suppose we wish to cook a nutritionally balanced meal by choosing the quantity
# of each food ``f`` to eat from a set of foods ``F`` in our kitchen.

# Each food ``f`` has a cost, ``c_f``, as well as a macro-nutrient profile
# ``a_{m,f}`` for each macro-nutrient ``m \in M``.

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
# a bit loose the term "macro-nutrient" by including calories and sodium.)

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
    ["nutrient", "min", "max"],
)

# ## JuMP formulation

# Now we're ready to convert our mathematical formulation into a JuMP model.

# First, create a new JuMP model. Since we have a linear program, we'll use
# HiGHS as our optimizer:

model = Model(HiGHS.Optimizer)
set_silent(model)

# Next, we create a set of decision variables `x`, with one element for each row
# in the DataFrame, and each `x` has a lower bound of `0`:

@variable(model, x[foods.name] >= 0)

# To simplify things later on, we store the vector as a new column `x` in the
# DataFrame `foods`:

foods.x = Array(x)

# Our objective is to minimize the total cost of purchasing food:

@objective(model, Min, sum(foods.cost .* foods.x));

# For the next component, we need to add a constraint that our total intake of
# each component is within the limits contained in the `limits` DataFrame:

@constraint(
    model,
    [row in eachrow(limits)],
    row.min <= sum(foods[!, row.nutrient] .* foods.x) <= row.max,
);

# What does our model look like?

print(model)

# ## Solution

# Let's optimize and take a look at the solution:

optimize!(model)
Test.@test primal_status(model) == FEASIBLE_POINT        #hide
Test.@test objective_value(model) ≈ 11.8288 atol = 1e-4  #hide
solution_summary(model)

# We found an optimal solution. Let's see what the optimal solution is:

for row in eachrow(foods)
    println(row.name, " = ", value(row.x))
end

# That's a lot of milk and ice cream, and sadly, we only get `0.6` of a
# hamburger.

# We can also use the function [`Containers.rowtable`](@ref) to easily convert
# the result into a DataFrame:

table = Containers.rowtable(value, x; header = [:food, :quantity])
solution = DataFrames.DataFrame(table)

# This makes it easy to perform analyses our solution:

filter!(row -> row.quantity > 0.0, solution)

# ## Problem modification

# JuMP makes it easy to take an existing model and modify it by adding extra
# constraints. Let's see what happens if we add a constraint that we can buy at
# most 6 units of milk or ice cream combined.

dairy_foods = ["milk", "ice cream"]
is_dairy = map(name -> name in dairy_foods, foods.name)
@constraint(model, sum(foods[is_dairy, :x]) <= 6)
optimize!(model)
Test.@test termination_status(model) == INFEASIBLE  #hide
Test.@test primal_status(model) == NO_SOLUTION      #hide
solution_summary(model)

# There exists no feasible solution to our problem. Looks like we're stuck
# eating ice cream for dinner.

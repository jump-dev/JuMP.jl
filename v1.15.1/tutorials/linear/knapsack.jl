# Copyright 2017, Iain Dunning, Joey Huchette, Miles Lubin, and contributors    #src
# This Source Code Form is subject to the terms of the Mozilla Public License   #src
# v.2.0. If a copy of the MPL was not distributed with this file, You can       #src
# obtain one at https://mozilla.org/MPL/2.0/.                                   #src

# # The knapsack problem example

# The purpose of this tutorial is to demonstrate how to formulate and solve a
# simple optimization problem.

# ## Required packages

# This tutorial requires the following packages:

using JuMP
import HiGHS

# ## Formulation

# The [knapsack problem](https://en.wikipedia.org/wiki/Knapsack_problem)
# is a classical optimization problem: given a set of items and a container with
# a fixed capacity, choose a subset of items having the greatest combined
# value that will fit within the container without exceeding the capacity.

# The name of the problem suggests its analogy to packing for a trip,
# where the baggage weight limit is the capacity and the goal is to pack the
# most profitable combination of belongings.

# We can formulate the knapsack problem as the integer linear program:
# ```math
# \begin{aligned}
# \max \; & \sum_{i=1}^n c_i x_i      \\
# s.t. \; & \sum_{i=1}^n w_i x_i \le C, \\
#         & x_i \in \{0,1\},\quad \forall i=1,\ldots,n,
# \end{aligned}
# ```
# where ``C`` is the capacity, and there is a choice between ``n`` items, with
# item ``i`` having weight ``w_i``, profit ``c_i``. Decision variable ``x_i`` is
# equal to 1 if the item is chosen and 0 if not.

# This formulation can be written more compactly as:
# ```math
# \begin{aligned}
# \max \; & c^\top x       \\
# s.t. \; & w^\top x \le C \\
#         & x \text{ binary }.
# \end{aligned}
# ```

# ## Data

# The data for the problem consists of two vectors (one for the profits and one
# for the weights) along with a capacity.

# There are five objects:

n = 5;

# For our example, we use a capacity of 10 units:

capacity = 10.0;

# and the profit and cost data:

profit = [5.0, 3.0, 2.0, 7.0, 4.0];
weight = [2.0, 8.0, 4.0, 2.0, 5.0];

# ## JuMP formulation

# Let's begin constructing the JuMP model for our knapsack problem.

# First, we'll create a `Model` object for holding model elements as we
# construct each part. We'll also set the solver that will ultimately be called
# to solve the model, once it's constructed.

model = Model(HiGHS.Optimizer)

# Next we need the decision variables representing which items are chosen:

@variable(model, x[1:n], Bin)

# We now want to constrain those variables so that their combined
# weight is less than or equal to the given capacity:

@constraint(model, sum(weight[i] * x[i] for i in 1:n) <= capacity)

# Finally, our objective is to maximize the combined profit of the chosen items:

@objective(model, Max, sum(profit[i] * x[i] for i in 1:n))

# Let's print a human-readable description of the model and check that the model
# looks as expected:

print(model)

# We can now solve the optimization problem and inspect the results.

optimize!(model)
solution_summary(model)

# The items chosen are

items_chosen = [i for i in 1:n if value(x[i]) > 0.5]

# ## Writing a function

# After working interactively, it is good practice to implement your model in a
# function.

# The function can be used to ensure that the model is given well-defined input
# data with validation checks, and that the solution process went as expected.

function solve_knapsack_problem(;
    profit::Vector{Float64},
    weight::Vector{Float64},
    capacity::Float64,
)
    n = length(weight)
    ## The profit and weight vectors must be of equal length.
    @assert length(profit) == n
    model = Model(HiGHS.Optimizer)
    set_silent(model)
    @variable(model, x[1:n], Bin)
    @objective(model, Max, profit' * x)
    @constraint(model, weight' * x <= capacity)
    optimize!(model)
    @assert termination_status(model) == OPTIMAL
    @assert primal_status(model) == FEASIBLE_POINT
    println("Objective is: ", objective_value(model))
    println("Solution is:")
    for i in 1:n
        print("x[$i] = ", round(Int, value(x[i])))
        println(", c[$i] / w[$i] = ", profit[i] / weight[i])
    end
    chosen_items = [i for i in 1:n if value(x[i]) > 0.5]
    return return chosen_items
end

solve_knapsack_problem(; profit = profit, weight = weight, capacity = capacity)

# We observe that the chosen items (1, 4, and 5) have the best
# profit to weight ratio in this particular example.

# ## Next steps

# Here are some things to try next:
#
# * Call the function with different data. What happens as the capacity
#   increases?
# * What happens if the profit and weight vectors are different lengths?
# * Instead of creating a binary variable with `Bin`, we could have written
#   `@variable(model, 0 <= x[1:n] <= 1, Int)`. Verify that this formulation
#   finds the same solution. What happens if we are allowed to take more than
#   one of each item?

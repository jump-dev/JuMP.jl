# Copyright 2017, Iain Dunning, Joey Huchette, Miles Lubin, and contributors    #src
# This Source Code Form is subject to the terms of the Mozilla Public License   #src
# v.2.0. If a copy of the MPL was not distributed with this file, You can       #src
# obtain one at https://mozilla.org/MPL/2.0/.                                   #src

# # The knapsack problem example

# The purpose of this example is demonstrate how to formulate and solve a 
# simple optimization problem.
# We use the knapsack problem for this purpose.

# ## Required packages

# This tutorial requires the following packages:

using JuMP
import HiGHS
import Test

# ## Formulation

# The simple [knapsack problem](https://en.wikipedia.org/wiki/Knapsack_problem) 
# is a well-known type of optimization problem: given a set of items and
# a container with a fixed capacity, choose a subset of items having the greatest combined
# value that will fit within the container without exceeding the capacity.
# The name of the problem suggests its analogy to packing for a trip,
# where the baggage weight limit is the capacity and the goal is to pack the
# most profitable combination of belongings.

# We can formulate the problem as an integer linear program

# ```math
# \begin{aligned}
# \max \; & \sum_{i=1}^n c_i x_i      \\
# s.t. \; & \sum_{i=1}^n w_i x_i \le C, \\
#         & x_i \in \{0,1\},\quad \forall i=1,\ldots,n
# \end{aligned}
# ```

# or more compactly as

# ```math
# \begin{aligned}
# \max \; & c^\top x       \\
# s.t. \; & w^\top x \le C \\
#         & x \text{ binary },
# \end{aligned}
# ```

# where there is a choice between ``n`` items, with item ``i`` having weight ``w_i``,
# profit ``c_i`` and a decision variable ``x_i`` equal to 1 if the item is chosen
# and 0 if not. 
# The capacity is a single real number ``C`` of the same number type as the 
# individual weights.

# ## Data

# The data for the problem consists of two vectors (one for the profits
# and one for the weights) along with a capacity.
# For our example, we use a capacity of 10 units
capacity = 10;
# and the vector data
profit = [5, 3, 2, 7, 4];
weight = [2, 8, 4, 2, 5];

# ## JuMP formulation

# Let's begin constructing the JuMP model for our knapsack problem.
# First, we'll create a `Model` object for holding model elements as
# we construct each part. We'll also set the solver that will
# ultimately be called to solve the model, once it's constructed.
model = Model(HiGHS.Optimizer)

# Next we need the decision variables representing which items are chosen.
@variable(model, x[1:5], Bin)

# We now want to constrain those variables so that their combined
# weight is less than or equal to the given capacity.
@constraint(model, sum(weight[i] * x[i] for i in 1:5) <= capacity)

# Finally, we set an objective: maximise the combined profit of the 
# selected items.
@objective(model, Max, sum(profit[i] * x[i] for i in 1:5))

# Let's print a human-readable description of the model and 
# check that the model looks as expected:
print(model)

# We can now solve the optimization problem and inspect the results.
optimize!(model)
solution_summary(model)

# The items chosen are

items_chosen = [i for i in 1:5 if value(x[i]) > 0.5]

# After working interactively, it is good practice to implement 
# your model in a function.
# A function can be used to ensure that the model is given
# well-defined input data by validation checks, and that the
# solution process went as expected.

function solve_knapsack_problem(;
    profit,
    weight::Vector{T},
    capacity::T,
) where {T<:Real}
    N = length(weight)

    ## The profit and weight vectors must be of equal length.
    @assert length(profit) == N

    model = Model(HiGHS.Optimizer)
    set_silent(model)

    ## Declare the decision variables as binary (0 or 1).
    @variable(model, x[1:N], Bin)

    ## Objective: maximize profit.
    @objective(model, Max, profit' * x)

    ## Constraint: can carry all items.
    @constraint(model, weight' * x <= capacity)

    ## Solve problem using MIP solver
    optimize!(model)
    println("Objective is: ", objective_value(model))
    println("Solution is:")
    for i in 1:N
        print("x[$i] = ", Int(value(x[i])))
        println(", c[$i]/w[$i] = ", profit[i] / weight[i])
    end
    Test.@test termination_status(model) == OPTIMAL
    Test.@test primal_status(model) == FEASIBLE_POINT
    Test.@test objective_value(model) == 16.0
    return
end

solve_knapsack_problem(; profit = profit, weight = weight, capacity = capacity)

# We observe that the chosen items (1, 4 and 5) have the best 
# profit to weight ratio in this particular example.

# Copyright 2017, Iain Dunning, Joey Huchette, Miles Lubin, and contributors    #src
# This Source Code Form is subject to the terms of the Mozilla Public License   #src
# v.2.0. If a copy of the MPL was not distributed with this file, You can       #src
# obtain one at https://mozilla.org/MPL/2.0/.                                   #src

# # Multi-objective knapsack

# This tutorial explains how to create and solve a multi-objective linear
# program. In addition, it demonstrates how to work with solvers which return
# multiple solutions.

# ## Required packages

# This tutorial requires the following packages:

using JuMP
import HiGHS
import MultiObjectiveAlgorithms as MOA
import Plots
import Test  #hide

# [`MultiObjectiveAlgorithms.jl`](https://github.com/jump-dev/MultiObjectiveAlgorithms.jl)
# is a package which implements a variety of algorithms for solving
# multi-objective optimization problems. Because it is a long package name, we
# import it instead as `MOA`.

# ## Formulation

# The [knapsack problem](https://en.wikipedia.org/wiki/Knapsack_problem) is a
# classic problem in mixed-integer programming. Given a collection of items
# ``i \in I``, each of which has an associated weight, ``w_i``, and profit,
# ``p_i``, the knapsack problem determines which profit-maximizing subset of
# items to pack into a knapsack such that the total weight is less than a
# capacity ``c``. The mathematical formulation is:

# ```math
# \begin{aligned}
# \max & \sum\limits_{i \in I} p_i x_i \\
# \text{s.t.}\ \ & \sum\limits_{i \in I} w_i x_i \le c\\
# & x_i \in \{0, 1\} && \forall i \in I
# \end{aligned}
# ```
# where ``x_i`` is ``1`` if we pack item ``i`` into the knapsack and ``0``
# otherwise.

# For this tutorial, we extend the single-objective knapsack problem by adding
# another objective: given a desirability rating, ``r_i``, we wish to maximize
# the total desirability of the items in our knapsack. Thus, our mathematical
# formulation is now:

# ```math
# \begin{aligned}
# \max & \sum\limits_{i \in I} p_i x_i \\
#      & \sum\limits_{i \in I} r_i x_i \\
# \text{s.t.}\ \ & \sum\limits_{i \in I} w_i x_i \le c\\
# & x_i \in \{0, 1\} && \forall i \in I
# \end{aligned}
# ```

# ## Data

# The data for this example was taken from [vOptGeneric](https://github.com/vOptSolver/vOptGeneric.jl),
# and the original author was [@xgandibleux](https://github.com/xgandibleux).

profit = [77, 94, 71, 63, 96, 82, 85, 75, 72, 91, 99, 63, 84, 87, 79, 94, 90]
desire = [65, 90, 90, 77, 95, 84, 70, 94, 66, 92, 74, 97, 60, 60, 65, 97, 93]
weight = [80, 87, 68, 72, 66, 77, 99, 85, 70, 93, 98, 72, 100, 89, 67, 86, 91]
capacity = 900
N = length(profit)

# Comparing the capacity to the total weight of all the items:

capacity / sum(weight)

# shows that we can take approximately 64% of the items.

# Plotting the items, we see that there are a range of items with different
# profits and desirability. Some items have a high profit and a high
# desirability, others have a low profit and a high desirability (and vice
# versa).

Plots.scatter(
    profit,
    desire;
    xlabel = "Profit",
    ylabel = "Desire",
    legend = false,
)

# The goal of the bi-objective knapsack problem is to choose a subset which
# maximizes both objectives.

# ## JuMP formulation

# Our JuMP formulation is a direct translation of the mathematical formulation:

model = Model()
@variable(model, x[1:N], Bin)
@constraint(model, sum(weight[i] * x[i] for i in 1:N) <= capacity)
@expression(model, profit_expr, sum(profit[i] * x[i] for i in 1:N))
@expression(model, desire_expr, sum(desire[i] * x[i] for i in 1:N))
@objective(model, Max, [profit_expr, desire_expr])

# Note how we form a multi-objective program by passing a vector of scalar
# objective functions.

# ## Solution

# To solve our model, we need an optimizer which supports multi-objective linear
# programs. One option is to use the [`MultiObjectiveAlgorithms.jl`](https://github.com/jump-dev/MultiObjectiveAlgorithms.jl)
# package.

set_optimizer(model, () -> MOA.Optimizer(HiGHS.Optimizer))
set_silent(model)

# `MultiObjectiveAlgorithms` supports many different algorithms for solving
# multiobjective optimization problems. One option is the epsilon-constraint
# method:

set_attribute(model, MOA.Algorithm(), MOA.EpsilonConstraint())

# Let's solve the problem and see the solution

optimize!(model)
solution_summary(model)

# There are 9 solutions available. We can also use [`result_count`](@ref) to see
# how many solutions are available:

result_count(model)

# ## Accessing multiple solutions

# Access the nine different solutions in the model using the `result` keyword to
# [`solution_summary`](@ref), [`value`](@ref), and [`objective_value`](@ref):

solution_summary(model; result = 5)

#-

objective_value(model; result = 5)

# Note that because we set a vector of two objective functions, the objective
# value is a vector with two elements. We can also query the value of each
# objective separately:

value(profit_expr; result = 5)

# ## Visualizing objective space

# Unlike single-objective optimization problems, multi-objective optimization
# problems do not have a single optimal solution. Instead, the solutions
# returned represent possible trade-offs that the decision maker can choose
# between the two objectives. A common way to visualize this is by plotting
# the objective values of each of the solutions:

plot = Plots.scatter(
    [value(profit_expr; result = i) for i in 1:result_count(model)],
    [value(desire_expr; result = i) for i in 1:result_count(model)];
    xlabel = "Profit",
    ylabel = "Desire",
    title = "Objective space",
    label = "",
    xlims = (915, 960),
)
for i in 1:result_count(model)
    y = objective_value(model; result = i)
    Plots.annotate!(y[1] - 1, y[2], (i, 10))
end
ideal_point = objective_bound(model)
Plots.scatter!([ideal_point[1]], [ideal_point[2]]; label = "Ideal point")

# Visualizing the objective space lets the decision maker choose a solution that
# suits their personal preferences. For example, result `#7` is close to the
# maximum value of profit, but offers significantly higher desirability compared
# with solutions `#8` and `#9`.

# The set of items that are chosen in solution `#7` are:

items_chosen = [i for i in 1:N if value(x[i]; result = 7) > 0.9]

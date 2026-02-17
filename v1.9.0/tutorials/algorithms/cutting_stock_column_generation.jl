# Copyright 2017, Iain Dunning, Joey Huchette, Miles Lubin, and contributors    #src
# This Source Code Form is subject to the terms of the Mozilla Public License   #src
# v.2.0. If a copy of the MPL was not distributed with this file, You can       #src
# obtain one at https://mozilla.org/MPL/2.0/.                                   #src

# # Column generation

# This tutorial describes how to implement the [Cutting stock problem](https://en.wikipedia.org/wiki/Cutting_stock_problem)
# in JuMP using column generation. It uses the following packages:

using JuMP
import GLPK
import SparseArrays
import Test  #src

# ## Mathematical formulation

# The cutting stock problem is about cutting large rolls of paper into smaller
# pieces. There is a demand different sizes of pieces to meet, and all large
# rolls have the same width. The goal is to meet the demand while maximizing the
# total profit.

# We denote the set of possible sized pieces that a roll can be cut into by
# ``i\in 1,\ldots,I``. Each piece ``i`` has a width, ``w_i``, and a demand,
# ``d_i``. The width of the large roll is ``W``.

# Here's the data that we are going to use in this tutorial:

struct Piece
    w::Float64
    d::Int
end

struct Data
    pieces::Vector{Piece}
    W::Float64
end

function Base.show(io::IO, d::Data)
    println(io, "Data for the cutting stock problem:")
    println(io, "  W = $(d.W)")
    println(io, "with pieces:")
    println(io, "   i   w_i d_i")
    println(io, "  ------------")
    for (i, p) in enumerate(d.pieces)
        println(io, lpad(i, 4), " ", lpad(p.w, 5), " ", lpad(p.d, 3))
    end
    return
end

function get_data()
    data = [
        75.0 38
        75.0 44
        75.0 30
        75.0 41
        75.0 36
        53.8 33
        53.0 36
        51.0 41
        50.2 35
        32.2 37
        30.8 44
        29.8 49
        20.1 37
        16.2 36
        14.5 42
        11.0 33
        8.6 47
        8.2 35
        6.6 49
        5.1 42
    ]
    return Data([Piece(data[i, 1], data[i, 2]) for i in axes(data, 1)], 100.0)
end

data = get_data()

# To formulate the cutting stock problem as a mixed-integer linear program, we
# assume that there is a set of large rolls ``j=1,\ldots,J`` to use. Then, we
# introduce two classes of decision variables:
# * ``x_{ij} \ge 0, \text{integer}, \forall i=1,\ldots,I, j=1,\ldots,J``
# * ``y_{j} \in \{0, 1\} \forall j=1,\ldots,J.``
# ``y_j`` is a binary variable that indicates if we use roll ``j``, and
# ``x_{ij}`` counts how many pieces of size ``i`` that we cut from roll ``j``.

# Our mixed-integer linear program is therefore:
# ```math
# \begin{align}
#     \min & \sum\limits_{j=1}^J y_j \\
#     \;\;\text{s.t.} & \sum\limits_{i=1}^N w_i x_{ij} \le W y_j & \forall j=1,\ldots,J \\
#          & \sum\limits_{j=1}^J x_{ij} \ge d_i & \forall i=1,\ldots,I \\
#          & x_{ij} \ge 0 & \forall i=1,\ldots,N, j=1,\ldots,J \\
#          & x_{ij} \in \mathbb{Z} & \forall i=1,\ldots,I, j=1,\ldots,J \\
#          & y_{j} \in \{0, 1\} & \forall j=1,\ldots,J \\
# \end{align}
# ```
# The objective is to minimize the number of rolls that we use, and the two
# constraints ensure that we respect the total width of each large roll and that
# we satisfy demand exactly.

I = length(data.pieces)
J = 1000  # Some large number
model = Model(GLPK.Optimizer)
@variable(model, x[1:I, 1:J] >= 0, Int)
@variable(model, y[1:J], Bin)
@constraint(
    model,
    [j in 1:J],
    sum(data.pieces[i].w * x[i, j] for i in 1:I) <= data.W * y[j],
)
@constraint(model, [i in 1:I], sum(x[i, j] for j in 1:J) >= data.pieces[i].d)
@objective(model, Min, sum(y[j] for j in 1:J))
model

# Unfortunately, we won't attempt to solve this formulation because it takes a
# very long time to solve. (Try it and see.)

## optimize!(model)

# However, there is a formulation that solves much faster, and that is to use a
# column generation scheme.

# ## Column generation theory

# The key insight for column generation is to recognize that the ``x`` variables
# above encode _cutting patterns_. For example, if we look only at the roll
# ``j=1``, then feasible solutions are:
# * ``x_{1,1} = 1``, ``x_{13,1} = 1`` and all the rest ``0``, which is 1 roll of
#   piece \#1 and 1 roll of piece \#13
# * ``x_{1,20} = 19`` and all the rest ``0``, which is 19 rolls of piece \#20.

# Cutting patterns like ``x_{1,1} = 1`` and ``x_{2,1} = 1`` are infeasible
# because the combined length is greater than ``W``.

# Since there are a finite number of ways that we could cut a roll into a
# valid cutting pattern, we can create a set of all possible cutting patterns
# ``p = 1,\ldots,P``, with data ``a_{i,p}`` indicating how many pieces of size
# ``i`` we cut in pattern ``p``. Then, we can formulate our mixed-integer linear
# program as:
# ```math
# \begin{align}
#     \min & \sum\limits_{p=1}^P x_p \\
#     \;\;\text{s.t.} & \sum\limits_{p=1}^P a_{ip} x_p \ge d_i & \forall i=1,\ldots,I \\
#          & x_{p} \ge 0 & \forall p=1,\ldots,P \\
#          & x_{p} \in \mathbb{Z} & \forall p=1,\ldots,P
# \end{align}
# ```

# Unfortunately, there will be a very large number of these patterns, so it is
# often intractable to enumerate all columns ``p=1,\ldots,P``.

# Column generation is an iterative algorithm that starts with a small set of
# initial patterns, and then cleverly chooses new columns to add to the main
# MILP so that we find the optimal solution without having to enumerate every
# column.

# ## Choosing new columns

# For now we assume that we have our mixed-integer linear program with a subset
# of the columns. If we have all of the columns that appear in an optimal
# solution then we are done. Otherwise, how do we choose a new column that leads
# to an improved solution?

# Column generation chooses a new column by relaxing the integrality constraint
# on ``x`` and looking at the dual variable ``\pi_i`` associated with demand
# constraint ``i``.

# Using the economic interpretation of the dual variable, we can say that a one
# unit increase in demand for piece ``i`` will cost an extra ``\pi_i`` rolls.
# Alternatively, we can say that a one unit increase in the left-hand side
# (for example, due to a new cutting pattern) will _save_ us ``\pi_i`` rolls.
# Therefore, we want a new column that maximizes the savings associated with
# the dual variables, while respecting the total width of the roll:
# ```math
# \begin{align}
#     \max & \sum\limits_{i=1}^I \pi_i y_i \\
#     \;\;\text{s.t.} & \sum\limits_{i=1}^I w_i y_{i} \le W \\
#          & y_{i} \ge 0 & \forall i=1,\ldots,I \\
#          & y_{i} \in \mathbb{Z} & \forall i=1,\ldots,I \\
# \end{align}
# ```
# If this problem, called the _pricing problem_, has an objective value greater
# than ``1``, then we estimate than adding `y` as the coefficients of a new
# column will decrease the objective by more than the cost of an extra roll.

# Here is code to solve the pricing problem:

function solve_pricing(data::Data, π::Vector{Float64})
    I = length(π)
    model = Model(GLPK.Optimizer)
    set_silent(model)
    @variable(model, y[1:I] >= 0, Int)
    @constraint(model, sum(data.pieces[i].w * y[i] for i in 1:I) <= data.W)
    @objective(model, Max, sum(π[i] * y[i] for i in 1:I))
    optimize!(model)
    if objective_value(model) > 1
        return round.(Int, value.(y))
    end
    return nothing
end

# ## Choosing the initial set of patterns

# For the initial set of patterns, we create a trivial cutting pattern which
# cuts as many pieces of size ``i`` as will fit, or the amount demanded,
# whichever is smaller.

patterns = Vector{Int}[]
for i in 1:I
    pattern = zeros(Int, I)
    pattern[i] = floor(Int, min(data.W / data.pieces[i].w, data.pieces[i].d))
    push!(patterns, pattern)
end
P = length(patterns)

# We can visualize the patterns by looking at the sparse matrix of the
# non-zeros:

SparseArrays.sparse(hcat(patterns...))

# ## Solving the problem

# First, we create our initial linear program:

model = Model(GLPK.Optimizer)
set_silent(model)
@variable(model, x[1:P] >= 0)
@objective(model, Min, sum(x))
@constraint(model, demand[i = 1:I], patterns[i]' * x == data.pieces[i].d)
model

# Then, we run the iterative column generation scheme:

while true
    ## Solve the linear relaxation
    optimize!(model)
    ## Obtain a new dual vector
    π = dual.(demand)
    ## Solve the pricing problem
    new_pattern = solve_pricing(data, π)
    ## Stop iterating if there is no new pattern
    if new_pattern === nothing
        break
    end
    push!(patterns, new_pattern)
    ## Create a new column
    push!(x, @variable(model, lower_bound = 0))
    ## Update the objective coefficients
    set_objective_coefficient(model, x[end], 1.0)
    ## Update the non-zeros in the coefficient matrix
    for i in 1:I
        if new_pattern[i] > 0
            set_normalized_coefficient(demand[i], x[end], new_pattern[i])
        end
    end
end

# Let's have a look at the patterns now:

SparseArrays.sparse(hcat(patterns...))

# Nice! We found over 20 new patterns.

# Here's pattern 21:

for i in 1:I
    if patterns[21][i] > 0
        println(patterns[21][i], " unit(s) of piece $i")
    end
end

# ## Looking at the solution

# Since we only solved a linear relaxation, some of our columns have fractional
# solutions. We can create a integer feasible solution by rounding up the orders:

for p in 1:length(x)
    v = ceil(Int, value(x[p]))
    if v > 0
        println(lpad(v, 2), " roll(s) of pattern $p")
    end
end

# This requires 343 rolls:

Test.@test sum(ceil.(Int, value.(x))) == 343  #src
sum(ceil.(Int, value.(x)))

# Alternatively, we can re-introduce the integrality constraints and resolve the
# problem:

set_integer.(x)
optimize!(model)
for p in 1:length(x)
    v = round(Int, value(x[p]))
    if v > 0
        println(lpad(v, 2), " roll(s) of pattern $p, each roll of which makes:")
        for i in 1:I
            if patterns[p][i] > 0
                println("  ", patterns[p][i], " unit(s) of piece $i")
            end
        end
    end
end

# This now requires 334 rolls:

Test.@test sum(ceil.(Int, value.(x))) == 334  #src
total_rolls = sum(ceil.(Int, value.(x)))

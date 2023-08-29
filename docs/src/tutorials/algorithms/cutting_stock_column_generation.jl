# Copyright 2017, Iain Dunning, Joey Huchette, Miles Lubin, and contributors    #src
# This Source Code Form is subject to the terms of the Mozilla Public License   #src
# v.2.0. If a copy of the MPL was not distributed with this file, You can       #src
# obtain one at https://mozilla.org/MPL/2.0/.                                   #src

# # Column generation

# This tutorial describes how to implement the [Cutting stock problem](https://en.wikipedia.org/wiki/Cutting_stock_problem)
# in JuMP using column generation. It uses the following packages:

using JuMP
import DataFrames
import HiGHS
import Plots
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
J = 1_000  # Some large number
model = Model(HiGHS.Optimizer)
set_silent(model)
@variable(model, x[1:I, 1:J] >= 0, Int)
@variable(model, y[1:J], Bin)
@constraint(
    model,
    [j in 1:J],
    sum(data.pieces[i].w * x[i, j] for i in 1:I) <= data.W * y[j],
)
@constraint(model, [i in 1:I], sum(x[i, :]) >= data.pieces[i].d)
@objective(model, Min, sum(y))
model

# Unfortunately, we cant't solve this formulation because it takes a very long
# time to solve. (Try removing the time limit.)

set_time_limit_sec(model, 5.0)
optimize!(model)
solution_summary(model)

# However, there is a formulation that solves much faster, and that is to use a
# column generation scheme.

# ## Column generation theory

# The key insight for column generation is to recognize that the ``x`` variables
# above encode _cutting patterns_. For example, if we look only at the roll
# ``j=1``, then feasible solutions are:
# * ``x_{1,1} = 1``, ``x_{13,1} = 1`` and all the rest ``0``, which is 1 roll of
#   piece \#1 and 1 roll of piece \#13
# * ``x_{20,1} = 19`` and all the rest ``0``, which is 19 rolls of piece \#20.

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
    model = Model(HiGHS.Optimizer)
    set_silent(model)
    @variable(model, y[1:I] >= 0, Int)
    @constraint(model, sum(data.pieces[i].w * y[i] for i in 1:I) <= data.W)
    @objective(model, Max, sum(π[i] * y[i] for i in 1:I))
    optimize!(model)
    if objective_value(model) <= 1
        return nothing  # Benefit of pattern is less than the cost of a new roll
    end
    return SparseArrays.sparse(round.(Int, value.(y)))
end

# If we solve the pricing problem with an artifical dual vector:

π = [1.0 / i for i in 1:I]
solve_pricing(data, π)

# the solution is a roll with 1 unit of size 1, 1 unit of size 17, and 3 units
# of size 20.

# If we solve the pricing problem with a dual vector of zeros, then the benefit
# of the new pattern is less than the cost of a roll, and so the function
# returns `nothing`:

solve_pricing(data, zeros(I))

# ## Choosing the initial set of patterns

# For the initial set of patterns, we create a trivial cutting pattern which
# cuts as many pieces of size ``i`` as will fit, or the amount demanded,
# whichever is smaller.

patterns = SparseArrays.SparseVector{Int,Int}[
    SparseArrays.sparsevec(
        [i],
        floor(Int, min(data.W / data.pieces[i].w, data.pieces[i].d)),
        I,
    ) for i in 1:I
]
P = length(patterns)

# We can visualize the patterns by looking at the sparse matrix of the
# non-zeros:

SparseArrays.sparse(hcat(patterns...))

# ## Solving the problem

# First, we create our initial linear program:

model = Model(HiGHS.Optimizer)
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
    rows, counts = SparseArrays.findnz(new_pattern)
    for (row, count) in zip(rows, counts)
        set_normalized_coefficient(demand[row], x[end], count)
    end
end

# Let's have a look at the patterns now:

SparseArrays.sparse(hcat(patterns...))

# We found over 20 new patterns.

# Here's pattern 21:

patterns[21]

# Here's another way to visualize the patterns:

function cutting_locations(data::Data, pattern::SparseArrays.SparseVector)
    locations = Float64[]
    offset = 0.0
    for (i, c) in zip(SparseArrays.findnz(pattern)...)
        for _ in 1:c
            offset += data.pieces[i].w
            push!(locations, offset)
        end
    end
    return locations
end

plot = Plots.bar(;
    xlims = (0, length(patterns) + 1),
    ylims = (0, data.W),
    xlabel = "Pattern",
    ylabel = "Roll length",
)
for (i, p) in enumerate(patterns)
    locations = cutting_locations(data, p)
    Plots.bar!(
        plot,
        fill(i, length(locations)),
        reverse(locations);
        bar_width = 1.0,
        label = false,
    )
end
plot

# ## Looking at the solution

# Since we only solved a linear relaxation, some of our columns have fractional
# solutions. We can create a integer feasible solution by rounding up the orders:

solution = DataFrames.DataFrame([
    (pattern = p, rolls = ceil(Int, value(x_p))) for (p, x_p) in enumerate(x)
])
filter!(row -> row.rolls > 0, solution)

# This requires 343 rolls:

Test.@test sum(solution.rolls) == 343  #src
sum(solution.rolls)

# Alternatively, we can re-introduce the integrality constraints and resolve the
# problem:

set_integer.(x)
optimize!(model)
solution = DataFrames.DataFrame([
    (pattern = p, rolls = round(Int, value(x_p))) for (p, x_p) in enumerate(x)
])
filter!(row -> row.rolls > 0, solution)

# This now requires 334 rolls:

Test.@test sum(solution.rolls) == 334  #src
sum(solution.rolls)

# Note that this may not be the global minimum because we are not adding new
# columns during the solution of the mixed-integer problem `model` (an algorithm
# known as [branch and price](https://en.wikipedia.org/wiki/Branch_and_price)).
# Nevertheless, the column generation algorithm typically finds good integer
# feasible solutions to an otherwise intractable optimization problem.

# ## Next steps

# * Our objective function is to minimize the total number of rolls. What is the
#   total length of waste? How does that compare to the total demand?
# * Writing the optimization algorithm is only part of the challenge. Can you
#   develop a better way to communicate the solution to stakeholders?

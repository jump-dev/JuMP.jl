# Copyright 2017, Iain Dunning, Joey Huchette, Miles Lubin, and contributors    #src
# This Source Code Form is subject to the terms of the Mozilla Public License   #src
# v.2.0. If a copy of the MPL was not distributed with this file, You can       #src
# obtain one at https://mozilla.org/MPL/2.0/.                                   #src

# # Constraint programming

# JuMP supports a range of constraint-programming type constraints via the
# corresponding sets in MathOptInterface. For most constraints, there are
# reformulations built-in that convert the constraint programming constraint
# into a mixed-integer programming equivalent.

# Because of this reformulation, all variables must be integer, and they must
# typically have finite bounds. An error will be thrown if the reformulation
# requires finiteness and you have a variable with non-finite bounds.

# This tutorial uses the following packages:

using JuMP
import HiGHS

# ## AllDifferent

# The [`MOI.AllDifferent`](@ref) set ensures that every element in a list takes
# a different integer value.

model = Model(HiGHS.Optimizer)
set_silent(model)
@variable(model, 1 <= x[1:4] <= 4, Int)
@constraint(model, x in MOI.AllDifferent(4))
optimize!(model)
value.(x)

# No two elements in `x` should have the same value!

# ## BinPacking

# The [`MOI.BinPacking`](@ref) set can be used to divide up a set of items into
# different groups, such that the sum of their weights does not exceed the
# capacity of a bin.

weights, capacity = Float64[1, 1, 2, 2, 3], 3.0;
number_of_bins = 3
model = Model(HiGHS.Optimizer)
set_silent(model)
@variable(model, 1 <= x[1:length(weights)] <= number_of_bins, Int)
@constraint(model, x in MOI.BinPacking(capacity, weights))
optimize!(model)
value.(x)

# Here, the value of `x[i]` is the bin that item `i` was placed into.

# ## Circuit

# The [`MOI.Circuit`](@ref) set is used to construct a tour of a list of `N`
# variables. They will each be assigned an integer from `1` to `N`, that
# describes the successor to each variable in the list:

model = Model(HiGHS.Optimizer)
set_silent(model)
@variable(model, x[1:4], Int)
@constraint(model, x in MOI.Circuit(4))
optimize!(model)

# Let's see what tour was found, starting at node number `1`:
y = round.(Int, value.(x))
tour = Int[1]
while length(tour) < length(y)
    push!(tour, y[tour[end]])
end
tour

# ## CountAtLeast

# The [`MOI.CountAtLeast`](@ref) set is used to ensure that at least `n`
# elements in a set of variables belong to a set of values.

# For example, here is a model with three variables, constrained between 0 and
# 5:

model = Model(HiGHS.Optimizer)
set_silent(model)
@variable(model, 0 <= x[1:3] <= 5, Int)

# If we want to ensure that at least one element of each set `{x[1], x[2]}` and
# `{x[2], x[3]}` is in the set `{3}`, then we create a list of variables by
# concatenating the sets together:

variables = [x[1], x[2], x[2], x[3]]

# Then we need a partition list that contains the number of elements in each
# set of variables:

partitions = [2, 2]

# Finally, we need a set of values that the elements must be a part of:

values = Set([3])

# And the number of elements that must be part of the set `values`:

n = 1

# The constraint is:

@constraint(model, variables in MOI.CountAtLeast(n, partitions, values))

# To ensure the uniqueness of the solution, we'll add a constraint that `x[2]`
# must be `<= 2`. This ensures that the only feasible solution is for `x[1]` and
# `x[3]` to be `3`:

@constraint(model, x[2] <= 2)

# Let's check that we found a valid solution:

optimize!(model)
value.(x)

# ## CountBelongs

# The [`MOI.CountBelongs`](@ref) set is used to count how many elements in a set
# of variables belong to a set of values.

# For example, to count how many elements in a set of 4 variables belong to the
# set `{2, 3}`, do:

model = Model(HiGHS.Optimizer)
set_silent(model)
@variable(model, 0 <= x[i = 1:4] <= i, Int)
@variable(model, n, Int)
@objective(model, Max, sum(x))
set = Set([2, 3])
@constraint(model, [n; x] in MOI.CountBelongs(1 + length(x), set))
optimize!(model)
value(n), value.(x)

# ## CountDistinct

# The [`MOI.CountDistinct`](@ref) set is used to count the number of distinct
# elements in a set of variables.

model = Model(HiGHS.Optimizer)
set_silent(model)
@variable(model, 0 <= x[i = 1:4] <= i, Int)
@variable(model, n, Int)
@objective(model, Max, sum(x))
@constraint(model, [n; x] in MOI.CountDistinct(1 + length(x)))
optimize!(model)
value(n), value.(x)

# ## CountGreaterThan

# The [`MOI.CountGreaterThan`](@ref) set is used to strictly upper-bound the
# number of distinct elements in a set of variables that have a value equal to
# another variable.

# For example, to count the number `n` of times that `y` appears in the vector
# `x`, use:

model = Model(HiGHS.Optimizer)
set_silent(model)
@variable(model, 0 <= x[i = 1:4] <= i, Int)
@variable(model, n, Int)
@variable(model, 3 <= y <= 4, Int)
@objective(model, Max, sum(x))
@constraint(model, [n; y; x] in MOI.CountGreaterThan(1 + 1 + length(x)))
optimize!(model)
value(n), value(y), value.(x)

# Here `n` is strictly greater than the count, and there is no limit on how
# large `n` could be. For example, `n = 100` is also a feasible solution. The
# only constraint is that `n` cannot be equal to or smaller than the number of
# times that `y` appears.

# ## Table

# The [`MOI.Table`](@ref) set is used to select a single row from a matrix of
# values.

# For example, given a matrix:

table = Float64[1 1 0; 0 1 1; 1 0 1; 1 1 1]

# we can constraint a 3-element vector `x` to equal one of the rows in `table`
# via:

model = Model(HiGHS.Optimizer)
set_silent(model)
@variable(model, x[i = 1:3], Int)
@constraint(model, x in MOI.Table(table))
optimize!(model)
value.(x)

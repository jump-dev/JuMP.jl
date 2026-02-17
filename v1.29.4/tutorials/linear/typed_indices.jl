# Copyright 2017, Iain Dunning, Joey Huchette, Miles Lubin, and contributors    #src
# This Source Code Form is subject to the terms of the Mozilla Public License   #src
# v.2.0. If a copy of the MPL was not distributed with this file, You can       #src
# obtain one at https://mozilla.org/MPL/2.0/.                                   #src

# # Strategies for dealing with many indices

# Common modeling patterns may lead to variables that are defined over many
# indices. For example, a transshipment model may have `x[p, f, m, t]` to
# represent the flow of product `p` from factory `f` to market `m` in time `t`.
# If each set is represented by a range of integers, it can be easy to make
# little errors like `x[t, p, f, m]` that are hard to catch and debug. The
# purpose of this tutorial is to explain how to used typed indices to avoid such
# bugs.

# ## Required packages

# This tutorial uses the following packages:

using JuMP

# ## The problem

# Consider a transshipment model of a single product from factories to markets.
# One way to model the flow variable is:

F, M = 3, 3
model = Model()
@variable(model, x[1:F, 1:M] >= 0)

# Now the question is: what does `x[1, 2]` represent? Was it `x[f, m]` or
# `x[m, f]`? In this simple example, it's probably easy to remember the flow
# _from_ a factory _to_ a market, but it gets harder if there are more sets. Was
# time the first index or the last?

# Mis-typing the order of two sets is a common error, and it can result in bugs
# that are surprisingly hard to find and debug. If you didn't already know,
# would you be able to spot the difference in these two constraints, especially
# if the constraint expression was large or you were looking through many
# constraints trying to find out why the solution wasn't what you expected?

@constraint(model, [f in 1:F], sum(x[f, m] for m in 1:M) == 1)

#-

@constraint(model, [f in 1:F], sum(x[m, f] for m in 1:M) == 1)

# There are two approaches you can use to avoid bugs like `x[m, f]`: keyword
# indexing and typed indices.

# ## Keyword indexing

# The first approach is to use [keyword indexing](@ref dense_keyword_indexing).
# Because keyword indexing does not work for `Array` containers, we must
# explicitly choose the `DenseAxisArray` container.

F, M = 3, 3
model = Model()
@variable(model, x[f in 1:F, m in 1:M] >= 0, container = DenseAxisArray)

# Accessing `x` in the correct order works:

x[f=1, m=2]

# But the wrong order errors:

try                         #hide
    x[m=2, f=1]
catch err                   #hide
    showerror(stderr, err)  #hide
end                         #hide

# Keyword indexing means we can write constraints like this:

@constraint(model, [f in 1:F], sum(x[f=f, m=m] for m in 1:M) == 1)

# ## Typed indices

# A second approach is to define new types for each index:

struct Factory
    x::Int
end

struct Market
    x::Int
end

# Then, we define each axis as a vector of these types:

factories = Factory.(1:F)
markets = Market.(1:M)
model = Model()
@variable(model, x[factories, markets] >= 0)

# Accessing `x` in the correct order works:

x[Factory(1), Market(2)]

# But the wrong order errors:

try                         #hide
    x[Market(2), Factory(1)]
catch err                   #hide
    showerror(stderr, err)  #hide
end                         #hide

# Typed indices means we can write constraints like this:

@constraint(model, [f in 1:F], sum(x[Factory(f), Market(m)] for m in 1:M) == 1)

# or like this:

@constraint(model, [f in factories], sum(x[f, m] for m in markets) == 1)

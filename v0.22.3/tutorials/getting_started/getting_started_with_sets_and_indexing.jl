# # Getting started with sets and indexing

# Most introductory courses to linear programming will teach you to identify
# sets over which the decision variables and constraints are indexed. Therefore,
# it is common to write variables such as ``x_i`` for all ``i \in I``.

# A common stumbling block for new users to JuMP is that _JuMP does not provide
# specialized syntax for constructing and manipulating these sets_.

# We made this decision because Julia already provides a wealth of data
# structures for working with sets.

# In contrast, because tools like AMPL are stand-alone software packages, they
# had to define their own syntax for set construction and manipulation. Indeed,
# the [AMPL Book](https://ampl.com/resources/the-ampl-book/chapter-downloads/)
# has two entire chapters devoted to sets and indexing (V: Simple Sets and
# Indexing, and VI: Compound Sets and Indexing).

# The purpose of this tutorial is to demonstrate a variety of ways in which you
# can construct and manipulate sets for optimization models.

# If you haven't already, you should first read
# [Getting started with JuMP](@ref).

using JuMP

# ## Unordered sets

# Unordered sets are useful to describe non-numeric indices, such as the names
# of cities or types of products.

# The most common way to construct a set is by creating a vector:

animals = ["dog", "cat", "chicken", "cow", "pig"]
model = Model()
@variable(model, x[animals])

# We can also use things like the `keys` of a dictionary:

weight_of_animals = Dict(
    "dog" => 20.0,
    "cat" => 5.0,
    "chicken" => 2.0,
    "cow" => 720.0,
    "pig" => 150.0,
)
animal_keys = keys(weight_of_animals)
model = Model()
@variable(model, x[animal_keys])

# A third option is to use Julia's `Set` object.

animal_set = Set()
for animal in keys(weight_of_animals)
    push!(animal_set, animal)
end
animal_set

# The nice thing about `Set`s is that they automatically remove duplicates:

push!(animal_set, "dog")
animal_set

# Note how `dog` does not appear twice.

model = Model()
@variable(model, x[animal_set])

# ## Sets of numbers

# Sets of numbers are useful to decribe sets that are ordered, such as years or
# elements in a vector. The easiest way to create sets of numbers is to use
# Julia's `range` syntax.

# These can start at `1`:

model = Model()
@variable(model, x[1:4])

# but they don't have to:

model = Model()
@variable(model, x[2012:2021])

# Ranges also have a `start:step:stop` syntax. So the Olympic years are:

model = Model()
@variable(model, x[1896:4:2020])

# ## Sets of other things

# An important observation is that you can have _any_ Julia type as the element
# of a set. It doesn't have to be a `String` or a `Number`. For example, you
# can have tuples:

sources = ["A", "B", "C"]
sinks = ["D", "E"]
S = [(source, sink) for source in sources, sink in sinks]
model = Model()
@variable(model, x[S])

#-

x[("A", "D")]

# For multi-dimensional sets, you can use JuMP's syntax for constructing
# [Containers](@ref):

model = Model()
@variable(model, x[sources, sinks])

#-

x["A", "D"]

# !!! info
#     Note how we indexed `x["A", "D"]` instead of `x[("A", "D")]` as above.

# ## Set operations

# Julia has built-in support for set operations such as `union`, `intersect`,
# and  `setdiff`.

# Therefore, to create a set of all years in which the summer Olympics were
# held, we can use:

baseline = 1896:4:2020
cancelled = [1916, 1940, 1944, 2020]
off_year = [2021]
olympic_years = union(setdiff(baseline, cancelled), off_year)

# You can also find the number of elements (i.e., the cardinality) in a set
# using `length`:

length(olympic_years)

# ## Set membership operations

# To compute membership of sets, use the `in` function.

2000 in olympic_years

#-

2001 in olympic_years

# ## Indexing expressions

# Use Julia's generator syntax to compute new sets, such as the list of
# Olympic years that are divisible by 3:

olympic_3_years = [year for year in olympic_years if mod(year, 3) == 0]
model = Model()
@variable(model, x[olympic_3_years])

# Alternatively, use JuMP's syntax for constructing [Containers](@ref):

model = Model()
@variable(model, x[year in olympic_years; mod(year, 3) == 0])

# ## Compound sets

# Consider a transportation problem in which we need to ship goods between
# cities. We have been provided a list of cities:

cities = ["Auckland", "Wellington", "Christchurch", "Dunedin"]

# and a distance matrix which records the shipping distance between pairs of
# cities. If we can't ship between two cities, the distance is `0`.

distances = [0 643 1071 1426; 0 0 436 790; 0 0 0 360; 1426 0 0 0]

# Let's have a look at ways we could write a model with an objective function to
# minimize the total shipping cost. For simplicity, we'll ignore all
# constraints.

# ### Fix unused variables

# One approach is to fix all variables that we can't use to zero. Most solvers
# are smart-enough to remove these during a presolve phase, so it has a very
# small impact on performance:

N = length(cities)
model = Model()
@variable(model, x[1:N, 1:N] >= 0)
for i in 1:N, j in 1:N
    if distances[i, j] == 0
        fix(x[i, j], 0.0; force = true)
    end
end
@objective(model, Min, sum(distances[i, j] * x[i, j] for i in 1:N, j in 1:N))

# ### Filtered summation

# Another approach is to define filters whenever we want to sum over our
# decision variables:

N = length(cities)
model = Model()
@variable(model, x[1:N, 1:N] >= 0)
@objective(
    model,
    Min,
    sum(
        distances[i, j] * x[i, j] for i in 1:N, j in 1:N if distances[i, j] > 0
    ),
)

# ### Filtered indexing

# We could also use JuMP's support for [Containers](@ref):

N = length(cities)
model = Model()
@variable(model, x[i = 1:N, j = 1:N; distances[i, j] > 0])
@objective(model, Min, sum(distances[i...] * x[i] for i in eachindex(x)))

# !!! note
#     The `i...` is called  a "splat". It converts a tuple like `(1, 2)` into
#     two indices like `distances[1, 2]`.

# ### Converting to a different data structure

# Another approach, and one that is often the most readable, is to convert the
# data you have into something that is easier to work with. Originally, we had
# a vector of strings and a matrix of distances. What we really need is
# something that maps usable origin-destination pairs to distances. A dictionary
# is an obvious choice:

routes = Dict(
    (a, b) => distances[i, j] for
    (i, a) in enumerate(cities), (j, b) in enumerate(cities) if
    distances[i, j] > 0
)

# Then, we can create our model like so:

model = Model()
@variable(model, x[keys(routes)])
@objective(model, Min, sum(v * x[k] for (k, v) in routes))

# This has a number of benefits over the other approaches, including a compacter
# algebraic model and variables that are named in a more meaningful way.

# !!! tip
#     If you're struggling to formulate a problem using the available syntax in
#     JuMP, it's probably a sign that you should convert your data into a
#     different form.

# ## Next steps

# The purpose of this tutorial was to show how JuMP does not have specialized
# syntax for set creation and manipulation. Instead, you should use the tools
# provided by Julia itself.

# This is both an opportunity and a challenge, because you are free to pick the
# syntax and data structures that best suit your problem, but for new users it
# can be daunting to decide which structure to use.

# Read through some of the other JuMP tutorials to get inspiration and ideas for
# how you can use Julia's syntax and data structures to your advantage.

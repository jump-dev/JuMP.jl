# Copyright 2017, Iain Dunning, Joey Huchette, Miles Lubin, and contributors    #src
# This Source Code Form is subject to the terms of the Mozilla Public License   #src
# v.2.0. If a copy of the MPL was not distributed with this file, You can       #src
# obtain one at https://mozilla.org/MPL/2.0/.                                   #src

# # The multi-commodity flow problem

# **This tutorial was originally contributed by Louis Luangkesorn.**

# This tutorial is a JuMP implementation of the multi-commodity transportation
# model described in
# [_AMPL: A Modeling Language for Mathematical Programming_](https://ampl.com/resources/the-ampl-book/),
# by R. Fourer, D.M. Gay and B.W. Kernighan.

# The purpose of this tutorial is to demonstrate creating a JuMP model from an
# SQLite database.

# ## Required packages

# This tutorial uses the following packages

using JuMP
import DataFrames
import HiGHS
import SQLite
import SQLite: DBInterface
import Tables
import Test


# ## Formulation

# The multi-commondity flow problem is a simple extension of
# [The transportation problem](@ref) to multiple types of products. Briefly, we
# start with the formulation of the transportation problem:
# ```math
# \begin{aligned}
# \min  && \sum_{i \in O, j \in D} c_{i,j} x_{i,j} \\
# s.t.  && \sum_{j \in D} x_{i, j} \le s_i && \forall i \in O \\
#       && \sum_{i \in O} x_{i, j} = d_j && \forall j \in D \\
#       && x_{i, j} \ge 0 && \forall i \in O, j \in D
# \end{aligned}
# ```
# but introduce a set of products ``P``, resulting in:
# ```math
# \begin{aligned}
# \min  && \sum_{i \in O, j \in D, k \in P} c_{i,j,k} x_{i,j,k} \\
# s.t.  && \sum_{j \in D} x_{i, j, k} \le s_{i,k} && \forall i \in O, k \in P \\
#       && \sum_{i \in O} x_{i, j, k} = d_{j,k} && \forall j \in D, k \in P \\
#       && x_{i, j,k} \ge 0 && \forall i \in O, j \in D, k \in P \\
#       && \sum_{k \in P} x_{i, j, k} \le u_{i,j} && \forall i \in O, j \in D
# \end{aligned}
# ```
# Note that the last constraint is new; it says that there is a maximum quantity
# of goods (of any type) that can be transported from origin ``i`` to
# destination ``j``.

# ## Data

# For the purpose of this tutorial, the JuMP repository contains an example
# database called `multi.sqlite`.

filename = joinpath(@__DIR__, "multi.sqlite");

# To run locally, download [`multi.sqlite`](multi.sqlite) and update `filename`
# appropriately.

# Load the database using `SQLite.DB`:

db = SQLite.DB(filename)

# A quick way to see the schema of the database is via `SQLite.tables`:

SQLite.tables(db)

# We interact with the database by executing queries, and then piping the
# results to an appropriate table. One example is a `DataFrame`:

DBInterface.execute(db, "SELECT * FROM locations") |> DataFrames.DataFrame

# But other table types are supported, such as `Tables.rowtable`:

DBInterface.execute(db, "SELECT * FROM locations") |> Tables.rowtable

# A `rowtable` is a `Vector` of `NamedTuple`s.

# You can construct more complicated SQL queries:

origins =
    DBInterface.execute(
        db,
        "SELECT location FROM locations WHERE type = \"origin\"",
    ) |> Tables.rowtable

# But for our purpose, we just want the list of strings:

origins = map(y -> y.location, origins)

# We can compose these two operations to get a list of destinations:

destinations =
    DBInterface.execute(
        db,
        "SELECT location FROM locations WHERE type = \"destination\"",
    ) |>
    Tables.rowtable |>
    x -> map(y -> y.location, x)

# And a list of products from our `products` table:

products =
    DBInterface.execute(db, "SELECT product FROM products") |>
    Tables.rowtable |>
    x -> map(y -> y.product, x)

# ## JuMP formulation

# We start by creating a model and our decision variables:

model = Model(HiGHS.Optimizer)
set_silent(model)
@variable(
    model,
    x[origin in origins, destination in destinations, product in products] >= 0,
    container = DataFrames.DataFrame,
)

# One approach when working with databases is to extract all of the data into a
# Julia datastructure. For example, let's pull the cost table into a DataFrame:

cost = DBInterface.execute(db, "SELECT * FROM cost") |> DataFrames.DataFrame

# and then join the decision variables:

function natural_join(left, right)
    on_names = intersect(names(left), names(right))
    return DataFrames.innerjoin(left, right; on = on_names)
end

cost_x = natural_join(cost, x)

# We've defined a new function, `natural_join`, to simplify the process of
# joining two DataFrames. This fuction acts like the `NATURAL JOIN` statment in
# SQL.

# Our objective is the inner product of two columns:

@objective(model, Max, cost_x.cost' * cost_x.value);

# The supply constraint is more complicated. A useful utility is a function that
# sums the `.value` column after grouping on a set of columns:

function sum_value_by(df, cols)
    gdf = DataFrames.groupby(df, cols)
    return DataFrames.combine(gdf, :value => sum => :value)
end

# Here is it in action:

sum_value_by(x, [:origin, :product])

# The constraint that the supply must be less than or equal to a capacity can
# now be written as:

supply = natural_join(
    DBInterface.execute(db, "SELECT * FROM supply") |> DataFrames.DataFrame,
    sum_value_by(x, [:origin, :product]),
)
@constraint(model, supply.value .<= supply.supply);

# The demand constraint ca be written similarly:

demand = natural_join(
    DBInterface.execute(db, "SELECT * FROM demand") |> DataFrames.DataFrame,
    sum_value_by(x, [:destination, :product]),
)
@constraint(model, demand.value .== demand.demand);

# The SQLite queries can be arbitrarily complex. For example, here's a query
# which builds every possible origin-destination pair:

od_pairs = DBInterface.execute(
    db,
    """
    SELECT a.location as 'origin',
           b.location as 'destination'
    FROM locations a
    INNER JOIN locations b
    ON a.type = 'origin' AND b.type = 'destination'
    """,
) |> DataFrames.DataFrame

# With a constraint that we cannot send more than 625 units between each pair:

od = natural_join(od_pairs, sum_value_by(x, [:origin, :destination]))
@constraint(model, od.value .<= 625);

# ## Solution

# Finally, we can optimize the model:

optimize!(model)
Test.@test is_solved_and_feasible(model)
Test.@test objective_value(model) == 225_700.0      #src
solution_summary(model)

# and obtain the solution:

x.value = value.(x.value)
x[x.value .> 0, :]

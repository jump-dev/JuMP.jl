# Copyright 2017, Iain Dunning, Joey Huchette, Miles Lubin, and contributors    #src
# This Source Code Form is subject to the terms of the Mozilla Public License   #src
# v.2.0. If a copy of the MPL was not distributed with this file, You can       #src
# obtain one at https://mozilla.org/MPL/2.0/.                                   #src

# # The multi-commodity flow problem

# **Originally contributed by:** Louis Luangkesorn

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
import Tables
import Test  #src

const DBInterface = SQLite.DBInterface

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
# database called `multi.sqlite`:

db = SQLite.DB(joinpath(@__DIR__, "multi.sqlite"))

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
@variable(model, x[origins, destinations, products] >= 0)

# One approach when working with databases is to extract all of the data into a
# Julia datastructure. For example, let's pull the cost table into a DataFrame
# and then construct our objective by iterating over the rows of the DataFrame:

cost = DBInterface.execute(db, "SELECT * FROM cost") |> DataFrames.DataFrame
@objective(
    model,
    Max,
    sum(r.cost * x[r.origin, r.destination, r.product] for r in eachrow(cost)),
);

# If we don't want to use a DataFrame, we can use a `Tables.rowtable` instead:

supply = DBInterface.execute(db, "SELECT * FROM supply") |> Tables.rowtable
for r in supply
    @constraint(model, sum(x[r.origin, :, r.product]) <= r.supply)
end

# Another approach is to execute the query, and then to iterate through the rows
# of the query using `Tables.rows`:

demand = DBInterface.execute(db, "SELECT * FROM demand")
for r in Tables.rows(demand)
    @constraint(model, sum(x[:, r.destination, r.product]) == r.demand)
end

# !!! warning
#     Iterating through the rows of a query result works by incrementing a
#     cursor inside the database. As a consequence, you cannot call
#     `Tables.rows` twice on the same query result.

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
)

# With a constraint that we cannot send more than 625 units between each pair:

for r in Tables.rows(od_pairs)
    @constraint(model, sum(x[r.origin, r.destination, :]) <= 625)
end

# ## Solution

# Finally, we can optimize the model:

optimize!(model)
Test.@test termination_status(model) == OPTIMAL     #src
Test.@test primal_status(model) == FEASIBLE_POINT   #src
Test.@test objective_value(model) == 225_700.0      #src
solution_summary(model)

# and print the solution:

begin
    println("         ", join(products, ' '))
    for o in origins, d in destinations
        v = lpad.([round(Int, value(x[o, d, p])) for p in products], 5)
        println(o, " ", d, " ", join(replace.(v, "   0" => "  . "), " "))
    end
end

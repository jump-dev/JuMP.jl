# Copyright 2017, Iain Dunning, Joey Huchette, Miles Lubin, and contributors    #src
# This Source Code Form is subject to the terms of the Mozilla Public License   #src
# v.2.0. If a copy of the MPL was not distributed with this file, You can       #src
# obtain one at https://mozilla.org/MPL/2.0/.                                   #src

# # The network multi-commodity flow problem

# This tutorial is a variation of [The multi-commodity flow problem](@ref) where
# the graph is a network instead of a bipartite graph.

# The purpose of this tutorial is to demonstrate a style of modeling that
# uses relational algebra.

# ## Required packages

# This tutorial uses the following packages:

using JuMP
import DataFrames
import HiGHS
import SQLite
import SQLite.DBInterface
import Test

# ## Formulation

# The network multi-commondity flow problem is an extension of the
# [The multi-commodity flow problem](@ref), where instead of having a bipartite
# graph of supply and demand nodes, the graph can contains a set of nodes,
# $i \in \mathcal{N}$ , which each have a (potentially zero) supply capacity,
# $u^s_{i,p}$, and (potentially zero) a demand, $d_{i,p}$ for each commodity
# $p \in P$. The nodes are connected by a set of edges $(i, j) \in \mathcal{E}$,
# which have a shipment cost $c^x_{i,j,p}$ and a total flow capacity of
# $u^x_{i,j}$.

# Our take is to choose an optimal supply for each node $s_{i,p}$, as well as
# the optimal transshipment $x_{i,j,p}$ that minimizes the total cost.

# The mathematical formulation is:
# ```math
# \begin{aligned}
# \min \;\; & \sum_{(i,j)\in\mathcal{E}, p \in P} c^x_{i,j,p} x_{i,j,p} + \sum_{i\in\mathcal{N}, p \in P} c^s_{i,p} s_{i,p} \\
# s.t. \;\; & s_{i,p} + \sum_{(j, i) \in \mathcal{E}} x_{j,i,p} - \sum_{(i,j) \in \mathcal{E}} x_{i,j,p} = d_{i,p} & \forall i \in \mathcal{N}, p \in P \\
#           & x_{i,j,p} \ge 0           & \forall (i, j) \in \mathcal{E}, p \in P \\
#           & \sum_{p \in P} x_{i,j,p} \le u^x_{i,j}           & \forall (i, j) \in \mathcal{E} \\
#           & 0 \le s_{i,p} \le u^s_{i,p} & \forall i \in \mathcal{N}, p \in P.
# \end{aligned}
# ```
#
# The purpose of this tutorial is to demonstrate how this model can be built
# using relational algebra instead of a direct math-to-code translation of the
# summations.

# ## Data

# For the purpose of this tutorial, the JuMP repository contains an example
# database called `commodity_nz.db`:

filename = joinpath(@__DIR__, "commodity_nz.db");

# To run locally, download [`commodity_nz.db`](commodity_nz.db) and update
# `filename` appropriately.

# Load the database using `SQLite.DB`:

db = SQLite.DB(filename)

# A quick way to see the schema of the database is via `SQLite.tables`:

SQLite.tables(db)

# We interact with the database by executing queries and then loading the
# results into a DataFrame:

function get_table(db, table)
    query = DBInterface.execute(db, "SELECT * FROM $table")
    return DataFrames.DataFrame(query)
end

# The `shipping` table contains the set of arcs and their distances:

df_shipping = get_table(db, "shipping")

# The `products` table contains the shipping cost per kilometer of each product:

df_products = get_table(db, "products")

# The `supply` table contains the supply capacity of each node, as well as the
# cost:

df_supply = get_table(db, "supply")

# The `demand` table contains the demand of each node:

df_demand = get_table(db, "demand")

# ## JuMP formulation

# We start by creating a model and our decision variables:

model = Model(HiGHS.Optimizer)
set_silent(model)

# For the shipping decisions, we create a new column in `df_shipping` called
# `x_flow`, which has one non-negative decision variable for each row:

df_shipping.x_flow = @variable(model, x[1:size(df_shipping, 1)] >= 0)
df_shipping

# For the supply, we add a variable to each row, and then set the upper bound to
# the capacity of each supply node:

df_supply.x_supply = @variable(model, s[1:size(df_supply, 1)] >= 0)
set_upper_bound.(df_supply.x_supply, df_supply.capacity)
df_supply

# Our objective is to minimize the shipping cost plus the supply cost. To
# compute the flow cost, we need to join the shipping table, which contains
# `distance_km` with the products table, which contains `cost_per_km`:

df_cost = DataFrames.leftjoin(df_shipping, df_products; on = [:product])
df_cost.flow_cost = df_cost.cost_per_km .* df_cost.distance_km
df_cost

# Then we can use linear algebra to compute the inner product between two
# columns:

@objective(
    model,
    Min,
    df_cost.flow_cost' * df_shipping.x_flow +
    df_supply.cost' * df_supply.x_supply
);

# For the flow capacities on each arc, we use `DataFrames.groupby` to partition
# the flow variables based on `:origin` and `:destination`, and then we
# constrain their sum to be less than a fixed capacity.

capacity = 30
for df in DataFrames.groupby(df_shipping, [:origin, :destination])
    @constraint(model, sum(df.x_flow) <= capacity)
end

# For each node in the graph, we need to compute a mass balance constraint
# which says that for each product, the supply, plus the flow into the node, and
# less the flow out of the node is equal to the demand.

# We can compute an expression for the flow out of each node using
# `DataFrames.groupby` on the `origin` and `product` columns of the
# `df_shipping` table:

df_flow_out = DataFrames.DataFrame(
    (node = i.origin, product = i.product, x_flow_out = sum(df.x_flow)) for
    (i, df) in pairs(DataFrames.groupby(df_shipping, [:origin, :product]))
)

# We can compute an expression for the flow into each node using
# `DataFrames.groupby` on the `destination` and `product` columns of the
# `df_shipping` table:

df_flow_in = DataFrames.DataFrame(
    (node = i.destination, product = i.product, x_flow_in = sum(df.x_flow))
    for (i, df) in
    pairs(DataFrames.groupby(df_shipping, [:destination, :product]))
)

# We can join the two tables together using `DataFrames.outerjoin`. We need to
# use `outerjoin` here because there might be missing rows.

df = DataFrames.outerjoin(df_flow_in, df_flow_out; on = [:node, :product])

# Next, we need to join the supply column:

df = DataFrames.leftjoin(
    df,
    DataFrames.select(df_supply, [:origin, :product, :x_supply]);
    on = [:node => :origin, :product],
)

# and then the demand column

df = DataFrames.leftjoin(
    df,
    DataFrames.select(df_demand, [:destination, :product, :demand]);
    on = [:node => :destination, :product],
)

# Now we're ready to add our mass balance constraint. Because some rows contain
# `missing` values, we need to use `coalesce` to convert any `missing` into a
# numeric value:

@constraint(
    model,
    [r in eachrow(df)],
    coalesce(r.x_supply, 0.0) + coalesce(r.x_flow_in, 0.0) -
    coalesce(r.x_flow_out, 0.0) == coalesce(r.demand, 0.0),
);

# ## Solution

# Finally, we can optimize the model:

optimize!(model)
Test.@test is_solved_and_feasible(model)
solution_summary(model)

# update the solution in the DataFrames:

df_shipping.x_flow = value.(df_shipping.x_flow)
df_supply.x_supply = value.(df_supply.x_supply);

# and display the optimal solution for flows:

DataFrames.select(
    filter!(row -> row.x_flow > 0.0, df_shipping),
    [:origin, :destination, :product, :x_flow],
)

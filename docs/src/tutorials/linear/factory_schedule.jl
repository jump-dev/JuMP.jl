# Copyright 2017, Iain Dunning, Joey Huchette, Miles Lubin, and contributors    #src
# This Source Code Form is subject to the terms of the Mozilla Public License   #src
# v.2.0. If a copy of the MPL was not distributed with this file, You can       #src
# obtain one at https://mozilla.org/MPL/2.0/.                                   #src

# # The factory schedule example

# This tutorial was originally contributed by `@Crghilardi`.

# This tutorial is a Julia translation of Part 5 from [Introduction to Linear
# Programming with Python](https://github.com/benalexkeen/Introduction-to-linear-programming).

# The purpose of this tutorial is to demonstrate how to use DataFrames and
# delimited files, and to structure your code that is robust to infeasibilities
# and permits running with different datasets.

# ## Required packages

using JuMP

import CSV
import DataFrames
import HiGHS
import StatsPlots

import Test  #src

# ## Formulation

# The Factory Scheduling Problem assumes we are optimizing the production of a
# good from factories ``f \in F`` over the course of 12 months ``m \in M``.

# If a factory ``f`` runs during a month ``m``, a fixed cost of ``a_f`` is
# incurred, the factory must produce ``x_{m,f}`` units that is within some
# minimum and maximum production levels ``l_f`` and ``u_f`` respectively, and
# each unit of production incurs a variable cost ``c_f``. Otherwise, the factory
# can be shut for the month with zero production and no fixed-cost is incurred.
# We denote the run/not-run decision by ``z_{m,f} \in \{0, 1\}``, where
# ``z_{m,f}`` is ``1`` if factory ``f`` runs in month ``m``. The factory must
# produce enough units to satisfy demand ``d_m``.

# With a little effort, we can formulate our problem as the following linear
# program:
# ```math
# \begin{aligned}
# \min & \sum\limits_{f \in F, m \in M} a_f z_{m,f} + c_f x_{m,f} \\
# \text{s.t.} & x_{m,f} \le u_f z_{m,f} && \forall f \in F, m \in M \\
#             & x_{m,f} \ge l_f z_{m,f} && \forall f \in F, m \in M \\
#             & \sum\limits_{f\in F} x_{m,f} = d_m && \forall f \in F, m \in M \\
#             & z_{m,f} \in \{0, 1\} && \forall f \in F, m \in M.
# \end{aligned}
# ```

# However, this formulation has a problem: if demand is too high, we may be
# unable to satisfy the demand constraint, and the problem will be infeasible.

# !!! tip
#     When modeling, consider ways to formulate your model such that it always
#     has a feasible solution. This greatly simplifies debugging data errors
#     that would otherwise result in an infeasible solution. In practice, most
#     practical decisions have a feasible solution. In our case, we could
#     satisfy demand (at a high cost) by buying replacement items for the buyer,
#     or running the factories in overtime to make up the difference.

# We can improve our model by adding a new variable, ``\delta_m``, which
# represents the quantity of unmet demand in each month ``m``. We penalize
# ``\delta_m`` by an arbitrarily large value of \$10,000/unit in the objective.

# ```math
# \begin{aligned}
# \min & \sum\limits_{f \in F, m \in M} a_f z_{m,f} + c_f x_{m,f} + \sum\limits_{m \in M}10000 \delta_m \\
# \text{s.t.} & x_{m,f} \le u_f z_{m,f} && \forall f \in F, m \in M \\
#             & x_{m,f} \ge l_f z_{m,f} && \forall f \in F, m \in M \\
#             & \sum\limits_{f\in F} x_{m,f} - \delta_m = d_m && \forall f \in F, m \in M \\
#             & z_{m,f} \in \{0, 1\} && \forall f \in F, m \in M \\
#             & \delta_m \ge 0 && \forall m \in M.
# \end{aligned}
# ```

# ## Data

# The JuMP GitHub repository contains two text files with the data we need for
# this tutorial.

# The first file contains a dataset of our factories, `A` and `B`, with their
# production and cost levels for each month. For the documentation, the file is
# located at:

factories_filename = joinpath(@__DIR__, "factory_schedule_factories.txt")

# and it has the following contents:

print(read(factories_filename, String))

# You can reproduce this locally by saving the contents to a new file and
# updating `factories_filename` appropriately.

# We use the `CSV` and `DataFrames` packages to read it into Julia:

factory_df = CSV.read(
    factories_filename,
    DataFrames.DataFrame;
    delim = ' ',
    ignorerepeated = true,
)

# The second file contains the demand data by month:

demand_filename = joinpath(@__DIR__, "factory_schedule_demand.txt")

print(read(demand_filename, String))

#-

demand_df = CSV.read(
    demand_filename,
    DataFrames.DataFrame;
    delim = ' ',
    ignorerepeated = true,
)

# Before moving on, it's always good practice to validate the data you read from
# external sources. The more effort you spend here, the fewer issues you will
# have later. The following function contains a few simple checks, but we could
# add others. For example, you might want to check that none of the values are
# too large (or too small), which might indicate a typo or a unit conversion
# issue (perhaps the variable costs are in \$/1000 units instead of \$/unit).

function valiate_data(
    demand_df::DataFrames.DataFrame,
    factory_df::DataFrames.DataFrame,
)
    ## Minimum production must not exceed maximum production.
    @assert all(factory_df.min_production .<= factory_df.max_production)
    ## Demand, minimum production, fixed costs, and variable costs must all be
    ## non-negative.
    @assert all(demand_df.demand .>= 0)
    @assert all(factory_df.min_production .>= 0)
    @assert all(factory_df.fixed_cost .>= 0)
    @assert all(factory_df.variable_cost .>= 0)
    return
end

valiate_data(demand_df, factory_df)

# ## JuMP formulation

# Next, we need to code our JuMP formulation. As shown in
# [Design patterns for larger models](@ref), it's always good practice to code
# your model in a function that accepts well-defined input and returns
# well-defined output.

function solve_factory_scheduling(
    demand_df::DataFrames.DataFrame,
    factory_df::DataFrames.DataFrame,
)
    ## Even though we validated the data above, it's good practice to do it here
    ## too.
    valiate_data(demand_df, factory_df)
    months, factories = unique(factory_df.month), unique(factory_df.factory)
    model = Model(HiGHS.Optimizer)
    set_silent(model)
    @variable(model, status[months, factories], Bin)
    @variable(model, production[months, factories], Int)
    @variable(model, unmet_demand[months] >= 0)
    ## We use `eachrow` to loop through the rows of the dataframe and add the
    ## relevant constraints.
    for row in eachrow(factory_df)
        m, f = row.month, row.factory
        @constraints(model, begin
            production[m, f] <= row.max_production * status[m, f]
            production[m, f] >= row.min_production * status[m, f]
        end)
    end
    @constraint(
        model,
        [row in eachrow(demand_df)],
        sum(production[row.month, :]) + unmet_demand[row.month] == row.demand,
    )
    @objective(
        model,
        Min,
        10_000 * sum(unmet_demand) +
        sum(
            row.fixed_cost * status[row.month, row.factory] +
            row.variable_cost * production[row.month, row.factory] for
            row in eachrow(factory_df)
        )
    )
    optimize!(model)
    schedules = Dict{Symbol,Vector{Float64}}(
        Symbol(f) => value.(production[:, f]) for f in factories
    )
    schedules[:unmet_demand] = value.(unmet_demand)
    return (
        termination_status = termination_status(model),
        cost = objective_value(model),
        ## This `select` statement re-orders the columns in the DataFrame.
        schedules = DataFrames.select(
            DataFrames.DataFrame(schedules),
            [:unmet_demand, :A, :B],
        ),
    )
end

# ## Solution

# Now we can call our `solve_factory_scheduling` function using the data we read
# in above.

solution = solve_factory_scheduling(demand_df, factory_df);

Test.@test solution.termination_status == OPTIMAL  #src
Test.@test solution.cost == 12_906_400.0  #src

# Let's see what `solution` contains:

solution.termination_status

#-

solution.cost

#-

solution.schedules

# These schedules will be easier to visualize as a graph:

StatsPlots.groupedbar(
    Matrix(solution.schedules);
    bar_position = :stack,
    labels = ["unmet demand" "A" "B"],
    xlabel = "Month",
    ylabel = "Production",
    legend = :topleft,
    color = ["#20326c" "#4063d8" "#a0b1ec"],
)

# Note that we don't have any unmet demand!

# ## What happens if demand increases?

# Let's run an experiment by increasing the demannd by 50% in all time periods:

demand_df.demand .*= 1.5

# Now we resolve the problem:

high_demand_solution = solve_factory_scheduling(demand_df, factory_df);

# and visualize the solution:

StatsPlots.groupedbar(
    Matrix(high_demand_solution.schedules);
    bar_position = :stack,
    labels = ["unmet demand" "A" "B"],
    xlabel = "Month",
    ylabel = "Production",
    legend = :topleft,
    color = ["#20326c" "#4063d8" "#a0b1ec"],
)

# Uh oh! We can't satisfy all of the demand.

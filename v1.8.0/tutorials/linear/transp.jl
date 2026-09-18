# Copyright 2017, Iain Dunning, Joey Huchette, Miles Lubin, and contributors    #src
# This Source Code Form is subject to the terms of the Mozilla Public License   #src
# v.2.0. If a copy of the MPL was not distributed with this file, You can       #src
# obtain one at https://mozilla.org/MPL/2.0/.                                   #src

# # The transportation problem

# **Originally contributed by:** Louis Luangkesorn

# This tutorial is an adaptation of the transportation problem described in
# [_AMPL: A Modeling Language for Mathematical Programming_](https://ampl.com/resources/the-ampl-book/),
# by R. Fourer, D.M. Gay and B.W. Kernighan.

# The purpose of this tutorial is to demonstrate how to create a JuMP model from
# an ad-hoc structured text file.

# ## Required packages

# This tutorial uses the following packages:

using JuMP
import DelimitedFiles
import HiGHS

# ## Formulation

# Suppose that we have a set of factories that produce [pogo sticks](https://en.wikipedia.org/wiki/Pogo_stick),
# and a set of retail stores in which to sell them. Each factory has a maximum
# number of pogo sticks that it can produce, and each retail store has a demand
# of pogo sticks that it can sell.

# In the transportation problem, we want to choose the number of pogo sticks to
# make and ship from each factory to each retail store that minimizes the total
# shipping cost.

# Mathematically, we represent our set of factories by a set of origins
# ``i \in O`` and our retail stores by a set of destinations ``j \in D``. The
# maximum supply at each factory is ``s_i`` and the demand from each retail
# store is ``d_j``. The cost of shipping one pogo stick from ``i`` to ``j`` is
# ``c_{i,j}``.

# With a little effort, we can model the transportation problem as the following
# linear program:
# ```math
# \begin{aligned}
# \min  && \sum_{i \in O, j \in D} c_{i,j} x_{i,j} \\
# s.t.  && \sum_{j \in D} x_{i, j} \le s_i && \forall i \in O \\
#       && \sum_{i \in O} x_{i, j} = d_j && \forall j \in D \\
#       && x_{i, j} \ge 0 && \forall i \in O, j \in D
# \end{aligned}
# ```

# ## Data

# We assume our data is in the form of a text file that has the following form.
# In practice, we would obtain this text file from the user as input, but for
# the purpose of this tutorial we're going to create it from Julia.

open(joinpath(@__DIR__, "transp.txt"), "w") do io
    print(
        io,
        """
             . FRA  DET LAN WIN  STL  FRE  LAF SUPPLY
          GARY  39   14  11  14   16   82    8   1400
          CLEV  27    .  12   .   26   95   17   2600
          PITT  24   14  17  13   28   99   20   2900
        DEMAND 900 1200 600 400 1700 1100 1000      0
        """,
    )
    return
end

# Here the rows are the origins, the columns are the destinations, and the
# values are the cost of shipping one pogo stick from the origin to the
# destination. If pogo stick cannot be transported from a source to a
# destination, then the value is `.`. The final row and final column are the
# demand and supply of each location respectively.

# We didn't account for arcs which do not exist in our formulation, but we can
# make a small change and fix ``x_{i,j} = 0`` if ``c_{i,j} = "."``.

# Our first step is to convert this text format into an appropriate Julia
# datastructure that we can work with. Since our data is tabular with named
# rows and columns, one option is JuMP's [`Containers.DenseAxisArray`](@ref)
# object:

function read_data(filename::String)
    data = DelimitedFiles.readdlm(filename)
    rows, columns = data[2:end, 1], data[1, 2:end]
    return Containers.DenseAxisArray(data[2:end, 2:end], rows, columns)
end

data = read_data(joinpath(@__DIR__, "transp.txt"))

# ## JuMP formulation

# Following [Design patterns for larger models](@ref), we code our JuMP model
# as a function which takes in an input. In this example, we print the output to
# `stdout`:

function solve_transportation_problem(data::Containers.DenseAxisArray)
    ## Get the set of supplies and demands
    O, D = axes(data)
    ## Drop the SUPPLY and DEMAND nodes from our sets
    O, D = setdiff(O, ["DEMAND"]), setdiff(D, ["SUPPLY"])
    model = Model(HiGHS.Optimizer)
    set_silent(model)
    @variable(model, x[o in O, d in D] >= 0)
    ## Remove arcs with "." cost by fixing them to 0.0.
    for o in O, d in D
        if data[o, d] == "."
            fix(x[o, d], 0.0; force = true)
        end
    end
    @objective(
        model,
        Min,
        sum(data[o, d] * x[o, d] for o in O, d in D if data[o, d] != "."),
    )
    @constraint(model, [o in O], sum(x[o, :]) <= data[o, "SUPPLY"])
    @constraint(model, [d in D], sum(x[:, d]) == data["DEMAND", d])
    optimize!(model)
    ## Pretty print the solution in the format of the input
    print("    ", join(lpad.(D, 7, ' ')))
    for o in O
        print("\n", o)
        for d in D
            if isapprox(value(x[o, d]), 0.0; atol = 1e-6)
                print("      .")
            else
                print(" ", lpad(value(x[o, d]), 6, ' '))
            end
        end
    end
    return
end

# ## Solution

# Let's solve and view the solution:

solve_transportation_problem(data)

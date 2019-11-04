#  Copyright 2017, Iain Dunning, Joey Huchette, Miles Lubin, and contributors
#  This Source Code Form is subject to the terms of the Mozilla Public
#  License, v. 2.0. If a copy of the MPL was not distributed with this
#  file, You can obtain one at http://mozilla.org/MPL/2.0/.
#############################################################################
# JuMP
# An algebraic modeling language for Julia
# See http://github.com/JuliaOpt/JuMP.jl
#############################################################################

using JuMP, GLPK, Test
const MOI = JuMP.MathOptInterface

"""
    example_cannery(; verbose = true)

JuMP implementation of the cannery problem from Dantzig, Linear Programming and
Extensions, Princeton University Press, Princeton, NJ, 1963.

Author: Louis Luangkesorn
Date: January 30, 2015
"""
function example_cannery(; verbose = true)
    plants = ["Seattle", "San-Diego"]
    markets = ["New-York", "Chicago", "Topeka"]

    # Capacity and demand in cases.
    capacity = [350, 600]
    demand = [300, 300, 300]

    # Distance in thousand miles.
    distance = [2.5 1.7 1.8; 2.5 1.8 1.4]

    # Cost per case per thousand miles.
    freight = 90

    num_plants = length(plants)
    num_markets = length(markets)

    cannery = Model(GLPK.Optimizer)

    @variable(cannery, ship[1:num_plants, 1:num_markets] >= 0)

    # Ship no more than plant capacity
    @constraint(cannery, capacity_con[i in 1:num_plants],
        sum(ship[i,j] for j in 1:num_markets) <= capacity[i]
    )

    # Ship at least market demand
    @constraint(cannery, demand_con[j in 1:num_markets],
        sum(ship[i,j] for i in 1:num_plants) >= demand[j]
    )

    # Minimize transporatation cost
    @objective(cannery, Min, sum(distance[i, j] * freight * ship[i, j]
        for i in 1:num_plants, j in 1:num_markets)
    )

    JuMP.optimize!(cannery)

    if verbose
        println("RESULTS:")
        for i in 1:num_plants
            for j in 1:num_markets
                println("  $(plants[i]) $(markets[j]) = $(JuMP.value(ship[i, j]))")
            end
        end
    end

    @test JuMP.termination_status(cannery) == MOI.OPTIMAL
    @test JuMP.primal_status(cannery) == MOI.FEASIBLE_POINT
    @test JuMP.objective_value(cannery) == 151200.0
end

example_cannery(verbose = false)

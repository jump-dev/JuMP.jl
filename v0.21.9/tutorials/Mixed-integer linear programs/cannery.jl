# Copyright 2017, Iain Dunning, Joey Huchette, Miles Lubin, and contributors    #src
# This Source Code Form is subject to the terms of the Mozilla Public License   #src
# v.2.0. If a copy of the MPL was not distributed with this file, You can       #src
# obtain one at https://mozilla.org/MPL/2.0/.                                   #src

# # The cannery problem

# A JuMP implementation of the cannery problem from:
#
# Dantzig, Linear Programming and Extensions, Princeton University Press,
# Princeton, NJ, 1963.
#
# It was originally contributed by Louis Luangkesorn, January 30, 2015.

using JuMP
import GLPK
import Test  #src

function example_cannery()
    ## Origin plants.
    plants = ["Seattle", "San-Diego"]
    num_plants = length(plants)
    ## Destination markets.
    markets = ["New-York", "Chicago", "Topeka"]
    num_markets = length(markets)
    ## Capacity and demand in cases.
    capacity = [350, 600]
    demand = [300, 300, 300]
    ## Distance in thousand miles.
    distance = [2.5 1.7 1.8; 2.5 1.8 1.4]
    ## Cost per case per thousand miles.
    freight = 90
    cannery = Model()
    set_optimizer(cannery, GLPK.Optimizer)
    ## Create decision variables.
    @variable(cannery, ship[1:num_plants, 1:num_markets] >= 0)
    ## Ship no more than plant capacity
    @constraint(
        cannery, capacity_con[i = 1:num_plants], sum(ship[i, :]) <= capacity[i]
    )
    ## Ship at least market demand
    @constraint(
        cannery, demand_con[j = 1:num_markets], sum(ship[:, j]) >= demand[j]
    )
    ## Minimize transportation cost
    @objective(
        cannery,
        Min,
        sum(
            distance[i, j] * freight * ship[i, j]
            for i = 1:num_plants, j = 1:num_markets
        )
    )
    optimize!(cannery)
    println("RESULTS:")
    for i = 1:num_plants
        for j = 1:num_markets
            println("  $(plants[i]) $(markets[j]) = $(value(ship[i, j]))")
        end
    end
    Test.@test termination_status(cannery) == MOI.OPTIMAL    #src
    Test.@test primal_status(cannery) == MOI.FEASIBLE_POINT  #src
    Test.@test objective_value(cannery) == 151_200.0         #src
    return
end

example_cannery()

#############################################################################
#  Copyright 2017, Iain Dunning, Joey Huchette, Miles Lubin, and contributors
#  This Source Code Form is subject to the terms of the Mozilla Public
#  License, v. 2.0. If a copy of the MPL was not distributed with this
#  file, You can obtain one at http://mozilla.org/MPL/2.0/.
#
# JuMP implementation of the cannery problem
# Dantzig, Linear Programming and Extensions,
# Princeton University Press, Princeton, NJ, 1963.
#
# Author: Louis Luangkesorn
# Date January 30, 2015
#############################################################################

using JuMP, Clp
const MOI = JuMP.MathOptInterface

solver = Clp.Optimizer

function PrintSolution(isoptimal, plants, markets, ship)
    println("RESULTS:")
    if isoptimal
      for i = 1:length(plants)
        for j in 1:length(markets)
          println("  $(plants[i]) $(markets[j]) = $(JuMP.result_value(ship[i,j]))")
        end
      end
    else
        println("  No solution")
    end
    println("")
end

function solveCannery(plants, markets, capacity, demand, distance, freight)
  numplants = length(plants)
  nummarkets = length(markets)
  cannery = Model(with_optimizer(solver))

  @variable(cannery, ship[1:numplants, 1:nummarkets] >= 0)

  # Ship no more than plant capacity
  @constraint(cannery, capacity_con[i=1:numplants],
               sum(ship[i,j] for j=1:nummarkets)<=capacity[i])

  # Ship at least market demand
  @constraint(cannery, demand_con[j=1:nummarkets],
               sum(ship[i,j] for i=1:numplants)>=demand[j])

  # Minimize transporatation cost
  @objective(cannery, Min,
              sum(distance[i,j]* freight * ship[i,j] for i=1:numplants, j=1:nummarkets))

  JuMP.optimize!(cannery)

  status = JuMP.termination_status(cannery)
  primal_status = JuMP.primal_status(cannery)
  isoptimal = status == MOI.Success && primal_status == MOI.FeasiblePoint

  PrintSolution(isoptimal, plants, markets, ship)
  return isoptimal
end


# PARAMETERS

plants = ["Seattle", "San-Diego"]
markets = ["New-York", "Chicago", "Topeka"]

# capacity and demand in cases
capacitycases = [350, 600]
demandcases = [300, 300, 300]

# distance in thousand miles
distanceKmiles = [2.5 1.7 1.8;
            2.5 1.8 1.4]

# cost per case per thousand miles
freightcost = 90

solveCannery(plants, markets, capacitycases, demandcases, distanceKmiles, freightcost)

#############################################################################
#  Copyright 2016, Iain Dunning, Joey Huchette, Miles Lubin, and contributors
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

using JuMP

function PrintSolution(status, plants, markets, ship)
    println("RESULTS:")
    if status == :Optimal
      for i = 1:length(plants)
        for j in 1:length(markets)
          println("  $(plants[i]) $(markets[j]) = $(getvalue(ship[i,j]))")
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
  cannery = Model()

  @variable(cannery, ship[1:numplants, 1:nummarkets] >= 0)

  # Ship no more than plant capacity
  @constraint(cannery, xyconstr[i=1:numplants],
                   sum{ship[i,j], j=1:nummarkets}<=capacity[i])

  # Ship at least market demand
  @constraint(cannery, xyconstr[j=1:nummarkets],
               sum{ship[i,j], i=1:numplants}>=demand[j])

  # Minimize transporatation cost
  @objective(cannery, Min,
              sum{distance[i,j]* freight * ship[i,j], i=1:numplants, j=1:nummarkets})

  solution = solve(cannery)
  if solution == :Optimal
    result = ship
    PrintSolution(solution, plants, markets, ship)
  else
    result = solution
    print("No solution")
  end
  result
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

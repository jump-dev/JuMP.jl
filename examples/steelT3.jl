################################################################################
#  Copyright 2017, Iain Dunning, Joey Huchette, Miles Lubin, and contributors
#  This Source Code Form is subject to the terms of the Mozilla Public
#  License, v. 2.0. If a copy of the MPL was not distributed with this
#  file, You can obtain one at http://mozilla.org/MPL/2.0/.
#
#   steelT3 model from
#   AMPL: A Modeling Language for Mathematical Programming, 2nd ed
#   by Robert Fourer, David Gay, and Brian W. Kernighan
#
#   Author : Louis Luangkesorn
#   Date: April 3, 2015
#
################################################################################

using JuMP, Clp
const MOI = JuMP.MathOptInterface
solver = Clp.Optimizer

# Data

T = 4
prod = ["bands", "coils"]
area = Dict([("bands", ("east", "north")),
             ("coils", ("east", "west", "export"))])

avail = [40, 40, 32,40]

rate = Dict("bands" => 200,
             "coils" =>140)

inv0 = Dict("bands" => 10,
             "coils" => 0)

prodcost = Dict("bands" => 10,
                 "coils" => 11)
invcost = Dict("bands" => 2.5,
                "coils" => 3)

revenue = Dict("bands" => Dict("east" => [25.0 26.0 27.0 27.0],
                                 "north" => [26.5 27.5 28.0 28.5]),
                "coils" => Dict("east" =>[30    35    37    39],
                                "west" => [29    32    33    35],
                                "export" => [25    25    25    28])
              )

market = Dict("bands" => Dict("east" => [2000  2000  1500  2000],
                                "north" => [4000  4000  2500  4500]),
               "coils" => Dict("east" => [1000   800  1000  1100],
                                "west" => [2000  1200  2000  2300],
                                "export" => [1000   500   500   800])
              )

# Model

Prod = Model(with_optimizer(solver))
print(area["bands"])
print(market["coils"]["west"][1])

# Decision Variables

@variable(Prod, Make[p in prod, t in 1:T] >= 0); # tons produced
@variable(Prod, Inv[p in prod, t in 0:T] >= 0); # tons inventoried
@variable(Prod, market[p][a][t] >= Sell[p in prod, a in area[p], t in 1:T] >= 0); # tons sold

@constraint(Prod, [p in prod, a in area[p], t in 1:T],
               Sell[p, a, t] - market[p][a][t] <= 0)


@constraint(Prod, [t in 1:T],
               sum((1/rate[p]) * Make[p,t] for p in prod) <= avail[t])

# Total of hours used by all products
# may not exceed hours available, in each week

@constraint(Prod, [p in prod],
               Inv[p,0] == inv0[p])
# Initial inventory must equal given value

@constraint(Prod, [p in prod, t in 1:T],
               Make[p,t] + Inv[p, t-1] == sum(Sell[p,a,t] for a in area[p]) + Inv[p,t])

# Tons produced and taken from inventory
# must equal tons sold and put into inventory


@objective(Prod, Max,
              sum( sum(
                   revenue[p][a][t] * Sell[p, a, t] -
                   prodcost[p] * Make[p,t] -
                   invcost[p]*Inv[p,t] for a in area[p]) for p in prod, t in 1:T))
#maximize Total_Profit:
# Total revenue less costs for all products in all weeks

function PrintSolution(is_optimal, area, Make, Inventory, Sell, product, Time)
    println("RESULTS:")
    if is_optimal
      for p in product
        println("Make $(p)")
        for t in 1:T
          print("$(JuMP.result_value(Make[p,t]))\t")
        end
        println()
        println("Inventory $(p)")
        for t in 1:T
          print("$(JuMP.result_value(Inventory[p,t]))\t")
        end
        println()
        for a in area[p]
          println("Sell $(p) $(a)")
          for t in 1:T
            print("$(JuMP.result_value(Sell[p,a,t])) \t")
          end
        println()
        end
      end
    else
      println("  No solution")
    end
end

JuMP.optimize!(Prod)

term_status = JuMP.termination_status(Prod)
primal_status = JuMP.primal_status(Prod)
is_optimal = status == MOI.Success && primal_status == MOI.FeasiblePoint

PrintSolution(is_optimal, area, Make, Inv, Sell, prod, T)

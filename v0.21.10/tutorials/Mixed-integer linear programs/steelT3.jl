# Copyright 2017, Iain Dunning, Joey Huchette, Miles Lubin, and contributors    #src
# This Source Code Form is subject to the terms of the Mozilla Public License   #src
# v.2.0. If a copy of the MPL was not distributed with this file, You can       #src
# obtain one at https://mozilla.org/MPL/2.0/.                                   #src

# # The SteelT3 problem

# The steelT3 model from AMPL: A Modeling Language for Mathematical Programming,
# 2nd ed by Robert Fourer, David Gay, and Brian W. Kernighan.
#
# Originally contributed by Louis Luangkesorn, April 3, 2015.

using JuMP
import GLPK
import Test

function example_steelT3(; verbose = true)
    T = 4
    prod = ["bands", "coils"]
    area = Dict(
        "bands" => ("east", "north"),
        "coils" => ("east", "west", "export"),
    )
    avail = [40, 40, 32, 40]
    rate = Dict("bands" => 200, "coils" => 140)
    inv0 = Dict("bands" => 10, "coils" => 0)
    prodcost = Dict("bands" => 10, "coils" => 11)
    invcost = Dict("bands" => 2.5, "coils" => 3)
    revenue = Dict(
        "bands" => Dict(
            "east" => [25.0, 26.0, 27.0, 27.0],
            "north" => [26.5, 27.5, 28.0, 28.5],
        ),
        "coils" => Dict(
            "east" => [30, 35, 37, 39],
            "west" => [29, 32, 33, 35],
            "export" => [25, 25, 25, 28],
        ),
    )
    market = Dict(
        "bands" => Dict(
            "east" => [2000, 2000, 1500, 2000],
            "north" => [4000, 4000, 2500, 4500],
        ),
        "coils" => Dict(
            "east" => [1000, 800, 1000, 1100],
            "west" => [2000, 1200, 2000, 2300],
            "export" => [1000, 500, 500, 800],
        ),
    )
    ## Model
    model = Model(GLPK.Optimizer)
    ## Decision Variables
    @variables(
        model,
        begin
            make[p in prod, t in 1:T] >= 0
            inventory[p in prod, t in 0:T] >= 0
            0 <= sell[p in prod, a in area[p], t in 1:T] <= market[p][a][t]
        end
    )
    @constraints(
        model,
        begin
            [p = prod, a = area[p], t = 1:T], sell[p, a, t] <= market[p][a][t]
            ## Total of hours used by all products may not exceed hours available,
            ## in each week
            [t in 1:T], sum(1 / rate[p] * make[p, t] for p in prod) <= avail[t]
            ## Initial inventory must equal given value
            [p in prod], inventory[p, 0] == inv0[p]
            ## Tons produced and taken from inventory must equal tons sold and put
            ## into inventory.
            [p in prod, t in 1:T],
            make[p, t] + inventory[p, t-1] ==
            sum(sell[p, a, t] for a in area[p]) + inventory[p, t]
        end
    )
    ## Maximize total profit: total revenue less costs for all products in all
    ## weeks.
    @objective(
        model,
        Max,
        sum(
            revenue[p][a][t] * sell[p, a, t] - prodcost[p] * make[p, t] -
            invcost[p] * inventory[p, t] for p in prod, a in area[p], t in 1:T
        )
    )
    optimize!(model)
    Test.@test termination_status(model) == MOI.OPTIMAL
    Test.@test primal_status(model) == MOI.FEASIBLE_POINT
    Test.@test objective_value(model) == 172850.0
    if verbose
        println("RESULTS:")
        for p in prod
            println("make $(p)")
            for t in 1:T
                print(value(make[p, t]), "\t")
            end
            println()
            println("Inventory $(p)")
            for t in 1:T
                print(value(inventory[p, t]), "\t")
            end
            println()
            for a in area[p]
                println("sell $(p) $(a)")
                for t in 1:T
                    print(value(sell[p, a, t]), "\t")
                end
                println()
            end
        end
    end
    return
end

example_steelT3()

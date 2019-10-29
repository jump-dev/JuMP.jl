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
    example_steelT3(; verbose = true)

steelT3 model from AMPL: A Modeling Language for Mathematical Programming, 2nd
ed by Robert Fourer, David Gay, and Brian W. Kernighan.

Author: Louis Luangkesorn
Date: April 3, 2015
"""
function example_steelT3(; verbose = true)
    T = 4
    prod = ["bands", "coils"]
    area = Dict(
        "bands" => ("east", "north"),
        "coils" => ("east", "west", "export")
    )
    avail = [40, 40, 32, 40]
    rate = Dict("bands" => 200, "coils" => 140)
    inv0 = Dict("bands" => 10, "coils" => 0)
    prodcost = Dict("bands" => 10, "coils" => 11)
    invcost = Dict("bands" => 2.5, "coils" => 3)
    revenue = Dict(
        "bands" => Dict("east" => [25.0, 26.0, 27.0, 27.0],
                        "north" => [26.5, 27.5, 28.0, 28.5]),
        "coils" => Dict("east" =>[30, 35, 37, 39],
                        "west" => [29, 32, 33, 35],
                        "export" => [25, 25, 25, 28])
    )
    market = Dict(
        "bands" => Dict(
            "east" => [2000, 2000, 1500, 2000],
            "north" => [4000, 4000, 2500, 4500]
        ),
        "coils" => Dict(
            "east" => [1000, 800, 1000, 1100],
            "west" => [2000, 1200, 2000, 2300],
            "export" => [1000, 500, 500, 800]
        )
    )

    # Model
    model = Model(GLPK.Optimizer)

    # Decision Variables
    @variables(model, begin
        make[p in prod, t in 1:T] >= 0
        inventory[p in prod, t in 0:T] >= 0
        0 <= sell[p in prod, a in area[p], t in 1:T] <= market[p][a][t]
    end)

    @constraint(model, [p in prod, a in area[p], t in 1:T],
        sell[p, a, t] - market[p][a][t] <= 0)

    # Total of hours used by all products may not exceed hours available, in
    # each week
    @constraint(model, [t in 1:T],
        sum(1 / rate[p] * make[p, t] for p in prod) <= avail[t])

    # Initial inventory must equal given value
    @constraint(model, [p in prod], inventory[p, 0] == inv0[p])

    # Tons produced and taken from inventory must equal tons sold and put into
    # inventory.
    @constraint(model, [p in prod, t in 1:T],
        make[p, t] + inventory[p, t - 1] == sum(sell[p, a, t] for a in area[p]) + inventory[p, t])

    # Maximize total profit: total revenue less costs for all products in all
    # weeks.
    @objective(model, Max, sum(
        sum(revenue[p][a][t] * sell[p, a, t] - prodcost[p] * make[p, t] -
            invcost[p] * inventory[p, t] for a in area[p])
        for p in prod, t in 1:T)
    )

    JuMP.optimize!(model)

    @test JuMP.termination_status(model) == MOI.OPTIMAL
    @test JuMP.primal_status(model) == MOI.FEASIBLE_POINT
    @test JuMP.objective_value(model) == 172850.0

    if verbose
        println("RESULTS:")
        for p in prod
            println("make $(p)")
            for t in 1:T
                print(JuMP.value(make[p, t]), "\t")
            end
            println()
            println("Inventory $(p)")
            for t in 1:T
                print(JuMP.value(inventory[p, t]), "\t")
            end
            println()
            for a in area[p]
                println("sell $(p) $(a)")
            for t in 1:T
                print(JuMP.value(sell[p, a, t]), "\t")
            end
            println()
            end
        end
    end
end

example_steelT3(verbose = false)

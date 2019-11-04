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
Author: Louis Luangkesorn
Date: February 26, 2015

JuMP implementation of the multicommodity transportation model AMPL: A Modeling
Language for Mathematical Programming, 2nd ed by Robert Fourer, David Gay, and
Brian W. Kernighan 4-1.
"""
function example_multi(; verbose = true)
    orig = ["GARY", "CLEV", "PITT"]
    dest = ["FRA", "DET", "LAN", "WIN", "STL", "FRE", "LAF"]
    prod = ["bands", "coils", "plate"]
    numorig = length(orig)
    numdest = length(dest)
    numprod = length(prod)

    # supply(prod, orig) amounts available at origins
    supply = [
        400    700    800;
        800   1600   1800;
        200    300    300
    ]

    # demand(prod, dest) amounts required at destinations
    demand = [
        300   300   100    75   650   225   250;
        500   750   400   250   950   850   500;
        100   100     0    50   200   100   250
    ]

    # limit(orig, dest) of total units from any origin to destination
    limit = [625.0 for j in 1:numorig, i in 1:numdest]

    # cost(dest, orig, prod) Shipment cost per unit
    cost = reshape([
        [
            [  30,   10,    8,   10,   11,   71,    6];
            [  22,    7,   10,    7,   21,   82,   13];
            [  19,   11,   12,   10,   25,   83,   15]
        ];
        [
            [  39,   14,   11,   14,   16,   82,    8];
            [  27,    9,   12,    9,   26,   95,   17];
            [  24,   14,   17,   13,   28,   99,   20]
        ];
        [
            [  41,   15,   12,   16,   17,   86,    8];
            [  29,    9,   13,    9,   28,   99,   18];
            [  26,   14,   17,   13,   31,  104,   20]
        ]
    ], 7, 3, 3)


    #  DECLARE MODEL
    multi = Model(GLPK.Optimizer)

    #  VARIABLES
    @variable(multi, trans[1:numorig, 1:numdest, 1:numprod] >= 0)

    #  OBJECTIVE
    @objective(multi, Max, sum(cost[j, i, p] * trans[i, j, p] for
        i in 1:numorig, j in 1:numdest, p in 1:numprod))

    #  CONSTRAINTS

    # Supply constraint
    @constraint(multi, supply_con[i in 1:numorig, p in 1:numprod],
        sum(trans[i, j, p] for j in 1:numdest) == supply[p, i])

    # Demand constraint
    @constraint(multi, demand_con[j in 1:numdest, p in 1:numprod],
        sum(trans[i, j, p] for i in 1:numorig) == demand[p, j])

    # Total shipment constraint
    @constraint(multi, total_con[i in 1:numorig, j in 1:numdest],
        sum(trans[i, j, p] for p in 1:numprod) - limit[i, j] <= 0)

    JuMP.optimize!(multi)

    @test JuMP.termination_status(multi) == MOI.OPTIMAL
    @test JuMP.primal_status(multi) == MOI.FEASIBLE_POINT
    @test JuMP.objective_value(multi) == 225700.0

    if verbose
        println("RESULTS:")
        for i in 1:length(ORIG)
            for j in 1:length(DEST)
                for p in 1:length(PROD)
                    print(" $(PROD[p]) $(ORIG[i]) $(DEST[j]) = $(JuMP.value(trans[i, j, p]))\t")
                end
                println()
            end
        end
    end
end


example_multi(verbose = false)

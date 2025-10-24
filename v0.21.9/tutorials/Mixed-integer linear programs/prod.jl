# Copyright 2017, Iain Dunning, Joey Huchette, Miles Lubin, and contributors    #src
# This Source Code Form is subject to the terms of the Mozilla Public License   #src
# v.2.0. If a copy of the MPL was not distributed with this file, You can       #src
# obtain one at https://mozilla.org/MPL/2.0/.                                   #src

# # The workforce scheduling problem
#
# This model determines a set of workforce levels that will most economically
# meet demands and inventory requirements over time. The formulation is
# motivated by the experiences of a large producer in the United States. The
# data are for three products and 13 periods.
#
# Problem taken from the Appendix C of the expanded version of Fourer, Gay, and
# Kernighan, A Modeling Language for Mathematical Programming
#
# Originally contributed by Louis Luangkesorn, February 26, 2015.

using JuMP
import GLPK
import Test

function example_prod(; verbose = true)
    ## PRODUCTION SETS AND PARAMETERS
    prd = ["18REG"  "24REG" "24PRO"]
    ## Members of the product group
    numprd = length(prd)
    pt =	[1.194,	1.509,	1.509]
    ## Crew-hours to produce 1000 units
    pc =	[2304,	2920,	2910]
    ## Nominal production cost per 1000, used
    ## to compute inventory and shortage costs
    ##
    ## TIME PERIOD SETS AND PARAMETERS
    firstperiod = 1
    ## Index of first production period to be modeled
    lastperiod  = 13
    ## Index of last production period to be modeled
    numperiods = firstperiod:lastperiod
    ## 'planning horizon' := first..last;
    ## EMPLOYMENT PARAMETERS
    ## Workers per crew
    cs = 18
    ## Regular-time hours per shift
    sl =  8
    ## Wage per hour for regular-time labor
    rtr = 16.00
    ## Wage per hour for overtime labor
    otr = 43.85
    ## Crews employed at start of first period
    iw =  8
    ## Regular working days in a production period
    dpp =	 [19.5,	19,	20,	19,	19.5,	19,	19,	20,	19,	20,	20,	18,	18]
    ## Maximum crew-hours of overtime in a period
    ol =	 [96,	96,	96,	96,	96,	96,	96,	96,	96,	96,	96,	96,	96]
    ## Lower limit on average employment in a period
    cmin =	[0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0]
    ## Upper limit on average employment in a period
    cmax =	[8,	8,	8,	8,	8,	8,	8,	8,	8,	8,	8,	8,	8]
    ## Penalty cost of hiring a crew
    hc =	 [7500,	7500,	7500,	7500,	15000,	15000,	15000,	15000,	15000,	15000,	7500,	7500,	7500]
    ## Penalty cost of laying off a crew
    lc =	 [7500,	7500,	7500,	7500,	15000,	15000,	15000,	15000,	15000,	15000,	7500,	7500,	7500]
    ## DEMAND PARAMETERS
    d18REG = [63.8,	76,	88.4,	913.8,	115,	133.8,	79.6,	111,	121.6,	470,	78.4,	99.4,	140.4,	63.8]
    d24REG = [1212,	306.2,	319,	208.4,	298,	328.2,	959.6,	257.6,	335.6,	118,	284.8,	970,	343.8,	1212]
    d24PRO = [0,	0,	0,	0,	0,	0,	0,	0,	0,	1102,	0,	0,	0,	0]
    ## Requirements (in 1000s) to be met from current production and inventory
    dem = Array[d18REG, d24REG, d24PRO]
    ## true if product will be the subject of a special promotion in the period
    pro = Array[
        [0,	0,	0,	1,	0,	0,	0,	0,	0,	1,	0,	0,	0,	0],
        [1,	0,	0,	0,	0,	0,	1,	0,	0,	0,	0,	0,	1,	1],
        [0,	0,	0,	0,	0,	0,	0,	0,	0,	1,	0,	0,	0,	0],
    ]
    ## INVENTORY AND SHORTAGE PARAMETERS
    ## Proportion of non-promoted demand that must be in inventory the previous
    ## period
    rir = 0.75
    ## Proportion of promoted demand that must be in inventory the previous
    ## period
    pir = 0.80
    ## Upper limit on number of periods that any product may sit in inventory
    life = 2
    ## Inventory cost per 1000 units is cri times nominal production cost
    cri	= [0.015,	0.015,	0.015]
    ## Shortage cost per 1000 units is crs times nominal production cost
    crs	= [1.1,	1.1,	1.1]
    ## Inventory at start of first period; age unknown
    iinv = [82,	792.2,	0]
    ## Initial inventory still available for allocation at end of period t
    iil = [
        [
            max(0, iinv[p] - sum(dem[p][v] for v in firstperiod:t))
            for t in numperiods
        ]
        for p in 1:numprd
    ]
    ## Lower limit on inventory at end of period t
    function checkpro(
        product, timeperiod, production, promotionalrate, regularrate
    )
        if production[product][timeperiod + 1] == 1
            return promotionalrate
        else
            return regularrate
        end
    end
    minv = [
        [dem[p][t + 1] * checkpro(p, t, pro, pir, rir) for t in numperiods]
        for p in 1:numprd
    ]
    ## DEFINE MODEL
    prod = Model(GLPK.Optimizer)
    ## VARIABLES
    ## Average number of crews employed in each period
    @variable(prod, Crews[0:lastperiod] >= 0)
    ## Crews hired from previous to current period
    @variable(prod, Hire[numperiods] >= 0)
    ## Crews laid off from previous to current period
    @variable(prod, Layoff[numperiods]>= 0)
    ## Production using regular-time labor, in 1000s
    @variable(prod, Rprd[1:numprd, numperiods] >= 0)
    ## Production using overtime labor, in 1000s
    @variable(prod, Oprd[1:numprd, numperiods]>= 0)
    ## a numperiods old -- produced in period (t+1)-a --
    ## and still in storage at the end of period t
    @variable(prod, Inv[1:numprd, numperiods, 1:life] >= 0)
    ## Accumulated unsatisfied demand at the end of period t
    @variable(prod, Short[1:numprd, numperiods] >= 0)
    ## CONSTRAINTS
    ## Hours needed to accomplish all regular-time production in a period must
    ## not exceed hours available on all shifts
    @constraint(
        prod,
        [t = numperiods],
        sum(pt[p] * Rprd[p, t] for p in 1:numprd) <= sl * dpp[t] * Crews[t]
    )
    ## Hours needed to accomplish all overtime production in a period must not
    ## exceed the specified overtime limit
    @constraint(
        prod,
        [t = numperiods],
        sum(pt[p] * Oprd[p, t] for p in 1:numprd)  <= ol[t]
    )
    ## Use given initial workforce
    @constraint(prod, Crews[firstperiod - 1] == iw)
    ## Workforce changes by hiring or layoffs
    @constraint(
        prod, [t in numperiods], Crews[t] == Crews[t - 1] + Hire[t] - Layoff[t]
    )
    ## Workforce must remain within specified bounds
    @constraint(prod, [t in numperiods], cmin[t] <= Crews[t])
    @constraint(prod, [t in numperiods], Crews[t] <= cmax[t])
    ## 'first demand requirement
    @constraint(
        prod,
        [p in 1:numprd],
        Rprd[p, firstperiod] + Oprd[p, firstperiod] + Short[p, firstperiod] -
            Inv[p, firstperiod, 1] == max(0, dem[p][firstperiod] - iinv[p])
    )
    ## Production plus increase in shortage plus decrease in inventory must
    ## equal demand
    for t in (firstperiod + 1):lastperiod
        @constraint(
            prod,
            [p in 1:numprd],
            Rprd[p, t] + Oprd[p, t] + Short[p,t] - Short[p,t-1] +
                sum(Inv[p, t - 1, a] - Inv[p, t, a] for a in 1:life) ==
                max(0, dem[p][t] - iil[p][t - 1])
        )
    end
    ## Inventory in storage at end of period t must meet specified minimum
    @constraint(
        prod,
        [p in 1:numprd, t in numperiods],
        sum(Inv[p, t, a] + iil[p][t] for a in 1:life) >= minv[p][t]
    )
    ## In the vth period (starting from first) no inventory may be more than v
    ## numperiods old (initial inventories are handled separately)
    @constraint(
        prod,
        [p in 1:numprd, v in 1:(life - 1), a in (v + 1):life],
        Inv[p, firstperiod + v - 1, a] == 0
    )
    ## New inventory cannot exceed production in the most recent period
    @constraint(
        prod,
        [p in 1:numprd, t in numperiods],
        Inv[p, t, 1] <= Rprd[p, t] + Oprd[p, t]
    )
    ## Inventory left from period (t+1)-p can only decrease as time goes on
    secondperiod = firstperiod + 1
    @constraint(
        prod,
        [p in 1:numprd, t in 2:lastperiod, a in 2:life],
        Inv[p, t, a] <= Inv[p, t - 1, a - 1]
    )
    ## OBJECTIVE
    ## Full regular wages for all crews employed, plus penalties for hiring and
    ## layoffs, plus wages for any overtime worked, plus inventory and shortage
    ## costs. (All other production costs are assumed to depend on initial
    ## inventory and on demands, and so are not included explicitly.)
    @objective(
        prod,
        Min,
        sum(
            rtr * sl * dpp[t] * cs * Crews[t] +
            hc[t] * Hire[t] +
            lc[t] * Layoff[t] +
            sum(
                otr * cs * pt[p] * Oprd[p, t] +
                sum(cri[p] * pc[p] * Inv[p, t, a] for a in 1:life) +
                crs[p] * pc[p] * Short[p, t]
                for p in 1:numprd
            )
           for t in numperiods
        )
    )
    ## Obtain solution
    optimize!(prod)
    Test.@test termination_status(prod) == MOI.OPTIMAL
    Test.@test primal_status(prod) == MOI.FEASIBLE_POINT
    Test.@test objective_value(prod) ≈ 4_426_822.89 atol = 1e-2
    if verbose
        println("RESULTS:")
        println("Crews")
        for t = 0:length(Crews.data) - 1
            print(" $(value(Crews[t])) ")
        end
        println()
        println("Hire")
        for t = 1:length(Hire.data)
            print(" $(value(Hire[t])) ")
        end
        println()
        println("Layoff")
        for t = 1:length(Layoff.data)
            print(" $(value(Layoff[t])) ")
        end
        println()
    end
    return
end

example_prod()

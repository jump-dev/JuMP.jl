# Copyright 2017, Iain Dunning, Joey Huchette, Miles Lubin, and contributors    #src
# This Source Code Form is subject to the terms of the Mozilla Public License   #src
# v.2.0. If a copy of the MPL was not distributed with this file, You can       #src
# obtain one at https://mozilla.org/MPL/2.0/.                                   #src

# # Mixed-integer linear programs

# These examples use the following packages:

using JuMP
import GLPK
import Test

# ## The cannery problem

# A JuMP implementation of the cannery problem from:
#
# Dantzig, Linear Programming and Extensions, Princeton University Press,
# Princeton, NJ, 1963.
#
# It was originally contributed by Louis Luangkesorn, January 30, 2015.

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
    ## Minimize transporatation cost
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

# # The diet problem

# Solve the classic "diet problem".
# Based on an [example from Gurobi](https://www.gurobi.com/documentation/9.0/examples/diet_cpp_cpp.html).

function example_diet(; verbose = true)
    function print_solution(is_optimal, foods, buy)
        println("RESULTS:")
        if is_optimal
            for food in foods
                println("  $(food) = $(value(buy[food]))")
            end
        else
            println("The solver did not find an optimal solution.")
        end
    end

    ## Nutrition guidelines
    categories = ["calories", "protein", "fat", "sodium"]
    category_data = Containers.DenseAxisArray([
        1800 2200;
        91   Inf;
        0    65;
        0    1779
        ], categories, ["min", "max"]
    )
    Test.@test category_data["protein", "min"] == 91.0
    Test.@test category_data["sodium", "max"] == 1779.0
    ## Foods
    foods = [
        "hamburger", "chicken", "hot dog", "fries", "macaroni", "pizza",
        "salad", "milk", "ice cream",
    ]
    cost = Containers.DenseAxisArray(
        [2.49, 2.89, 1.50, 1.89, 2.09, 1.99, 2.49, 0.89, 1.59],
        foods
    )
    food_data = Containers.DenseAxisArray(
        [
            410 24 26 730;
            420 32 10 1190;
            560 20 32 1800;
            380  4 19 270;
            320 12 10 930;
            320 15 12 820;
            320 31 12 1230;
            100  8 2.5 125;
            330  8 10 180
        ], foods, categories
    )
    Test.@test food_data["hamburger", "calories"] == 410.0
    Test.@test food_data["milk", "fat"] == 2.5
    ## Build model
    model = Model(GLPK.Optimizer)
    @variables(model, begin
        ## Variables for nutrition info
        category_data[c, "min"] <= nutrition[c = categories] <= category_data[c, "max"]
        ## Variables for which foods to buy
        buy[foods] >= 0
    end)
    ## Objective - minimize cost
    @objective(model, Min, sum(cost[f] * buy[f] for f in foods))
    ## Nutrition constraints
    @constraint(model, [c in categories],
        sum(food_data[f, c] * buy[f] for f in foods) == nutrition[c]
    )
    ## Solve
    if verbose
        println("Solving original problem...")
    end
    optimize!(model)
    term_status = termination_status(model)
    is_optimal = term_status == MOI.OPTIMAL
    Test.@test primal_status(model) == MOI.FEASIBLE_POINT
    Test.@test objective_value(model) ≈ 11.8288 atol = 1e-4
    if verbose
        print_solution(is_optimal, foods, buy)
    end
    ## Limit dairy (note that the problem will become infeasible).
    @constraint(model, buy["milk"] + buy["ice cream"] <= 6)
    if verbose
        println("Solving dairy-limited problem...")
    end
    optimize!(model)
    Test.@test termination_status(model) == MOI.INFEASIBLE
    Test.@test primal_status(model) == MOI.NO_SOLUTION
    if verbose
        print_solution(false, foods, buy)
    end
    return
end

example_diet()

# ## The knapsack problem

# Formulate and solve a simple knapsack problem:
#
#     max sum(p_j x_j)
#      st sum(w_j x_j) <= C
#         x binary

function example_knapsack(; verbose = true)
    profit = [5, 3, 2, 7, 4]
    weight = [2, 8, 4, 2, 5]
    capacity = 10
    model = Model(GLPK.Optimizer)
    @variable(model, x[1:5], Bin)
    ## Objective: maximize profit
    @objective(model, Max, profit' * x)
    ## Constraint: can carry all
    @constraint(model, weight' * x <= capacity)
    ## Solve problem using MIP solver
    optimize!(model)
    if verbose
        println("Objective is: ", objective_value(model))
        println("Solution is:")
        for i in 1:5
            print("x[$i] = ", value(x[i]))
            println(", p[$i]/w[$i] = ", profit[i] / weight[i])
        end
    end
    Test.@test termination_status(model) == MOI.OPTIMAL
    Test.@test primal_status(model) == MOI.FEASIBLE_POINT
    Test.@test objective_value(model) == 16.0
    return
end

example_knapsack()

# ## The SteelT3 problem

# The steelT3 model from AMPL: A Modeling Language for Mathematical Programming,
# 2nd ed by Robert Fourer, David Gay, and Brian W. Kernighan.
#
# Originally contributed by Louis Luangkesorn, April 3, 2015.

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
        "bands" => Dict(
            "east" => [25.0, 26.0, 27.0, 27.0],
            "north" => [26.5, 27.5, 28.0, 28.5],
        ),
        "coils" => Dict(
            "east" =>[30, 35, 37, 39],
            "west" => [29, 32, 33, 35],
            "export" => [25, 25, 25, 28],
        )
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
        )
    )
    ## Model
    model = Model(GLPK.Optimizer)
    ## Decision Variables
    @variables(model, begin
        make[p in prod, t in 1:T] >= 0
        inventory[p in prod, t in 0:T] >= 0
        0 <= sell[p in prod, a in area[p], t in 1:T] <= market[p][a][t]
    end)
    @constraints(model, begin
        [p = prod, a = area[p], t = 1:T], sell[p, a, t] <= market[p][a][t]
        ## Total of hours used by all products may not exceed hours available,
        ## in each week
        [t in 1:T], sum(1 / rate[p] * make[p, t] for p in prod) <= avail[t]
        ## Initial inventory must equal given value
        [p in prod], inventory[p, 0] == inv0[p]
        ## Tons produced and taken from inventory must equal tons sold and put
        ## into inventory.
        [p in prod, t in 1:T], make[p, t] + inventory[p, t - 1] == sum(sell[p, a, t] for a in area[p]) + inventory[p, t]
    end)
    ## Maximize total profit: total revenue less costs for all products in all
    ## weeks.
    @objective(
        model,
        Max,
        sum(
            revenue[p][a][t] * sell[p, a, t] -
            prodcost[p] * make[p, t] -
            invcost[p] * inventory[p, t]
            for p in prod, a in area[p], t in 1:T
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

# ## The transportation problem

# Allocation of passenger cars to trains to minimize cars required or car-miles
# run. Based on:
#
# Fourer, D.M. Gay and Brian W. Kernighan, A Modeling Language for Mathematical
# Programming, https://ampl.com/REFS/amplmod.ps.gz Appendix D.
#
# Originally contributed by Louis Luangkesorn, January 30, 2015.

function example_transp()
    ORIG = ["GARY", "CLEV", "PITT"]
    DEST = ["FRA", "DET", "LAN", "WIN", "STL", "FRE", "LAF"]
    supply = [1_400, 2_600, 2_900]
    demand = [900, 1_200, 600, 400, 1_700, 1_100, 1_000]
    Test.@test sum(supply) == sum(demand)
    cost = [
        39   14   11   14   16   82    8;
        27    9   12    9   26   95   17;
        24   14   17   13   28   99   20
    ]
    model = Model(GLPK.Optimizer)
    @variable(model, trans[1:length(ORIG), 1:length(DEST)] >= 0)
    @objective(
        model,
        Min,
        sum(
            cost[i, j] * trans[i, j]
            for i in 1:length(ORIG), j in 1:length(DEST)
        )
    )
    @constraints(model, begin
        [i in 1:length(ORIG)], sum(trans[i, :]) == supply[i]
        [j in 1:length(DEST)], sum(trans[:, j]) == demand[j]
    end)
    optimize!(model)
    Test.@test termination_status(model) == MOI.OPTIMAL
    Test.@test primal_status(model) == MOI.FEASIBLE_POINT
    Test.@test objective_value(model) == 196200.0
    println("The optimal solution is:")
    println(value.(trans))
    return
end

example_transp()

# ## The urban planning problem

# An "urban planning" problem based on an [example from puzzlor](http://www.puzzlor.com/2013-08_UrbanPlanning.html).

function example_urban_plan()
    model = Model(GLPK.Optimizer)
    ## x is indexed by row and column
    @variable(model, 0 <= x[1:5, 1:5] <= 1, Int)
    ## y is indexed by R or C, the points, and an index in 1:5. Note how JuMP
    ## allows indexing on arbitrary sets.
    rowcol = ["R", "C"]
    points = [5, 4, 3, -3, -4, -5]
    @variable(model, 0 <= y[rowcol, points, 1:5] <= 1, Int)
    ## Objective - combine the positive and negative parts
    @objective(model, Max, sum(
          3 * (y["R", 3, i] + y["C", 3, i])
        + 1 * (y["R", 4, i] + y["C", 4, i])
        + 1 * (y["R", 5, i] + y["C", 5, i])
        - 3 * (y["R", -3, i] + y["C", -3, i])
        - 1 * (y["R", -4, i] + y["C", -4, i])
        - 1 * (y["R", -5, i] + y["C", -5, i])
        for i in 1:5)
    )
    ## Constrain the number of residential lots
    @constraint(model, sum(x) == 12)
    ## Add the constraints that link the auxiliary y variables to the x variables
    for i = 1:5
        @constraints(model, begin
            ## Rows
            y["R", 5, i] <= 1 / 5 * sum(x[i, :]) # sum = 5
            y["R", 4, i] <= 1 / 4 * sum(x[i, :]) # sum = 4
            y["R", 3, i] <= 1 / 3 * sum(x[i, :]) # sum = 3
            y["R", -3, i] >= 1 - 1 / 3 * sum(x[i, :]) # sum = 2
            y["R", -4, i] >= 1 - 1 / 2 * sum(x[i, :]) # sum = 1
            y["R", -5, i] >= 1 - 1 / 1 * sum(x[i, :]) # sum = 0
            ## Columns
            y["C", 5, i] <= 1 / 5 * sum(x[:, i]) # sum = 5
            y["C", 4, i] <= 1 / 4 * sum(x[:, i]) # sum = 4
            y["C", 3, i] <= 1 / 3 * sum(x[:, i]) # sum = 3
            y["C", -3, i] >= 1 - 1 / 3 * sum(x[:, i]) # sum = 2
            y["C", -4, i] >= 1 - 1 / 2 * sum(x[:, i]) # sum = 1
            y["C", -5, i] >= 1 - 1 / 1 * sum(x[:, i]) # sum = 0
        end)
    end
    ## Solve it
    optimize!(model)
    Test.@test termination_status(model) == MOI.OPTIMAL
    Test.@test primal_status(model) == MOI.FEASIBLE_POINT
    Test.@test objective_value(model) ≈ 14.0
    return
end

example_urban_plan()

# ## The factory schedule example

# This is a Julia translation of part 5 from "Introduction to to Linear
# Programming with Python"
# available at https://github.com/benalexkeen/Introduction-to-linear-programming
#
# For 2 factories (A, B), minimize the cost of production over the course of 12
# months while meeting monthly demand. Factory B has a planned outage during
# month 5.
#
# It was originally contributed by @Crghilardi.

function example_factory_schedule()
   ## Sets in the problem:
   months, factories = 1:12, [:A, :B]
   ## This function takes a matrix and converts it to a JuMP container so we can
   ## refer to elements such as `d_max_cap[1, :A]`.
   containerize(A::Matrix) = Containers.DenseAxisArray(A, months, factories)
   ## Maximum production capacity in (month, factory) [units/month]:
   d_max_cap = containerize([
         100000	50000;
         110000	55000;
         120000	60000;
         145000	100000;
         160000	0;
         140000	70000;
         155000	60000;
         200000	100000;
         210000	100000;
         197000	100000;
         80000	120000;
         150000	150000;
   ])
   ## Minimum production capacity in (month, factory) [units/month]:
   d_min_cap = containerize([
         20000	20000;
         20000	20000;
         20000	20000;
         20000	20000;
         20000	0;
         20000	20000;
         20000	20000;
         20000	20000;
         20000	20000;
         20000	20000;
         20000	20000;
         20000	20000;
   ])
   ## Variable cost of production in (month, factory) [$/unit]:
   d_var_cost = containerize([
         10	5;
         11	4;
         12	3;
         9	5;
         8	0;
         8	6;
         5	4;
         7	6;
         9	8;
         10	11;
         8	10;
         8	12
   ])
   ## Fixed cost of production in (month, factory) # [$/month]:
   d_fixed_cost = containerize([
         500	600;
         500	600;
         500	600;
         500	600;
         500	0;
         500	600;
         500	600;
         500	600;
         500	600;
         500	600;
         500	600;
         500	600
   ])
   ## Demand in each month [units/month]:
   d_demand = [
      120_000,
      100_000,
      130_000,
      130_000,
      140_000,
      130_000,
      150_000,
      170_000,
      200_000,
      190_000,
      140_000,
      100_000,
   ]
   ## The model!
   model = Model(GLPK.Optimizer)
   ## Decision variables
   @variables(model, begin
      status[m in months, f in factories], Bin
      production[m in months, f in factories], Int
   end)
   ## The production cannot be less than minimum capacity.
   @constraint(
      model,
      [m in months, f in factories],
      production[m, f] >= d_min_cap[m, f] * status[m, f],
   )
   ## The production cannot be more tha maximum capacity.
   @constraint(
      model,
      [m in months, f in factories],
      production[m, f] <= d_max_cap[m, f] * status[m, f],
   )
   ## The production must equal demand in a given month.
   @constraint(model, [m in months], sum(production[m, :]) == d_demand[m])
   ## Factory B is shut down during month 5, so production and status are both
   ## zero.
   fix(status[5, :B], 0.0)
   fix(production[5, :B], 0.0)
   ## The objective is to minimize the cost of production across all time
   ##periods.
   @objective(
      model,
      Min,
      sum(
         d_fixed_cost[m, f] * status[m, f] + d_var_cost[m, f] * production[m, f]
         for m in months, f in factories
      )
   )
   ## Optimize the problem
   optimize!(model)
   ## Check the solution!
   Test.@testset "Check the solution against known optimal" begin
      Test.@test termination_status(model) == MOI.OPTIMAL
      Test.@test objective_value(model) == 12_906_400.0
      Test.@test value.(production)[1, :A] == 70_000
      Test.@test value.(status)[1, :A] == 1
      Test.@test value.(status)[5, :B] == 0
      Test.@test value.(production)[5, :B] == 0
   end
   println("The production schedule is:")
   println(value.(production))
   return
end

example_factory_schedule()

# ## The multi-commodity flow problem

# JuMP implementation of the multicommodity transportation model AMPL: A Modeling
# Language for Mathematical Programming, 2nd ed by Robert Fourer, David Gay, and
# Brian W. Kernighan 4-1.
#
# Originally contributed by Louis Luangkesorn, February 26, 2015.

using JuMP
import GLPK
import Test

function example_multi(; verbose = true)
    orig = ["GARY", "CLEV", "PITT"]
    dest = ["FRA", "DET", "LAN", "WIN", "STL", "FRE", "LAF"]
    prod = ["bands", "coils", "plate"]
    numorig = length(orig)
    numdest = length(dest)
    numprod = length(prod)
    ## supply(prod, orig) amounts available at origins
    supply = [
        400    700    800;
        800   1600   1800;
        200    300    300
    ]
    ## demand(prod, dest) amounts required at destinations
    demand = [
        300   300   100    75   650   225   250;
        500   750   400   250   950   850   500;
        100   100     0    50   200   100   250
    ]
    ## limit(orig, dest) of total units from any origin to destination
    limit = [625.0 for j in 1:numorig, i in 1:numdest]
    ## cost(dest, orig, prod) Shipment cost per unit
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
    ## DECLARE MODEL
    multi = Model(GLPK.Optimizer)
    ## VARIABLES
    @variable(multi, trans[1:numorig, 1:numdest, 1:numprod] >= 0)
    ## OBJECTIVE
    @objective(
        multi,
        Max,
        sum(
            cost[j, i, p] * trans[i, j, p]
            for i in 1:numorig, j in 1:numdest, p in 1:numprod
        )
    )
    ## CONSTRAINTS
    ## Supply constraint
    @constraint(
        multi,
        supply_con[i in 1:numorig, p in 1:numprod],
        sum(trans[i, j, p] for j in 1:numdest) == supply[p, i]
    )
    ## Demand constraint
    @constraint(
        multi,
        demand_con[j in 1:numdest, p in 1:numprod],
        sum(trans[i, j, p] for i in 1:numorig) == demand[p, j]
    )
    ## Total shipment constraint
    @constraint(
        multi,
        total_con[i in 1:numorig, j in 1:numdest],
        sum(trans[i, j, p] for p in 1:numprod) - limit[i, j] <= 0
    )
    optimize!(multi)
    Test.@test termination_status(multi) == MOI.OPTIMAL
    Test.@test primal_status(multi) == MOI.FEASIBLE_POINT
    Test.@test objective_value(multi) == 225_700.0
    if verbose
        println("RESULTS:")
        for i in 1:length(orig)
            for j in 1:length(dest)
                for p in 1:length(prod)
                    print(" $(prod[p]) $(orig[i]) $(dest[j]) = $(value(trans[i, j, p]))\t")
                end
                println()
            end
        end
    end
    return
end

example_multi()

# ## The workforce scheduling problem
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

# ## Solving Sudokus with MIP

# A sudoku solver that uses a MIP to find solutions.

# We have binary variables `x[i, j, k]` which, if = 1, say that cell (i, j)
# contains the number k. The constraints are:
#  1 - Each cell has one value only
#  2 - Each row contains each number exactly once
#  3 - Each column contains each number exactly once
#  4 - Each 3x3 subgrid contains each number exactly once
# We will take the initial grid as a CSV file at `filepath`, where 0s are blanks.

using JuMP
import GLPK
import Test

function example_sudoku(filepath::String)
    initial_grid = zeros(Int, 9, 9)
    open(filepath, "r") do fp
        for row in 1:9
            line = readline(fp)
            initial_grid[row, :] .= parse.(Int, split(line, ","))
        end
    end
    model = Model(GLPK.Optimizer)
    @variable(model, x[1:9, 1:9, 1:9], Bin)
    @constraints(model, begin
        ## Constraint 1 - Only one value appears in each cell
        cell[i in 1:9, j in 1:9], sum(x[i, j, :]) == 1
        ## Constraint 2 - Each value appears in each row once only
        row[i in 1:9, k in 1:9], sum(x[i, :, k]) == 1
        ## Constraint 3 - Each value appears in each column once only
        col[j in 1:9, k in 1:9], sum(x[:, j, k]) == 1
        ## Constraint 4 - Each value appears in each 3x3 subgrid once only
        subgrid[i=1:3:7, j=1:3:7, val=1:9], sum(x[i:i + 2, j:j + 2, val]) == 1
    end)
    ## Initial solution
    for row in 1:9, col in 1:9
        if initial_grid[row, col] != 0
            fix(x[row, col, initial_grid[row, col]], 1)
        end
    end
    ## Solve it
    optimize!(model)
    ## Check solution
    term_status = termination_status(model)
    is_optimal = term_status == MOI.OPTIMAL
    if is_optimal
        mip_solution = value.(x)
        sol = zeros(Int, 9, 9)
        for row in 1:9, col in 1:9, val in 1:9
            if mip_solution[row, col, val] >= 0.9
                sol[row, col] = val
            end
        end
        return sol
    else
        error("The solver did not find an optimal solution.")
    end
end

# Create an initial file. We could have set the `example_sudoku` function to
# take a matrix as input, but this example shows Julia's ability to parse files.

open("sudoku.csv", "w") do io
    write(io, """
    3, 1, 0, 0, 5, 8, 0, 0, 4
    0, 0, 9, 3, 2, 0, 0, 0, 0
    0, 2, 5, 1, 0, 4, 0, 9, 0
    0, 0, 0, 0, 0, 0, 3, 8, 9
    0, 0, 8, 0, 0, 0, 5, 0, 0
    5, 4, 6, 0, 0, 0, 0, 0, 0
    0, 8, 0, 2, 0, 3, 6, 5, 0
    0, 0, 0, 0, 7, 1, 4, 0, 0
    7, 0, 0, 4, 8, 0, 0, 2, 1
    """)
end

# Now try solving the example:

solution = example_sudoku("sudoku.csv")

Test.@test solution == [
    3 1 7 9 5 8 2 6 4;
    4 6 9 3 2 7 8 1 5;
    8 2 5 1 6 4 7 9 3;
    2 7 1 6 4 5 3 8 9;
    9 3 8 7 1 2 5 4 6;
    5 4 6 8 3 9 1 7 2;
    1 8 4 2 9 3 6 5 7;
    6 9 2 5 7 1 4 3 8;
    7 5 3 4 8 6 9 2 1
]

# The matrix can be hard to read. Print the solution properly:

function print_sudoku_solution(solution)
    println("Solution:")
    println("[-----------------------]")
    for row in 1:9
        print("[ ")
        for col in 1:9
            print(solution[row, col], " ")
            if col % 3 == 0 && col < 9
                print("| ")
            end
        end
        println("]")
        if row % 3 == 0
            println("[-----------------------]")
        end
    end
end

print_sudoku_solution(solution)

# Clean up the file we made:

rm("sudoku.csv")

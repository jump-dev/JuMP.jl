# Copyright 2017, Iain Dunning, Joey Huchette, Miles Lubin, and contributors    #src
# This Source Code Form is subject to the terms of the Mozilla Public License   #src
# v.2.0. If a copy of the MPL was not distributed with this file, You can       #src
# obtain one at https://mozilla.org/MPL/2.0/.                                   #src

# # The factory schedule example

# This is a Julia translation of part 5 from "Introduction to to Linear
# Programming with Python"
# available at https://github.com/benalexkeen/Introduction-to-linear-programming
#
# For 2 factories (A, B), minimize the cost of production over the course of 12
# months while meeting monthly demand. Factory B has a planned outage during
# month 5.
#
# It was originally contributed by @Crghilardi.

using JuMP
import GLPK
import Test

function example_factory_schedule()
    ## Sets in the problem:
    months, factories = 1:12, [:A, :B]
    ## This function takes a matrix and converts it to a JuMP container so we can
    ## refer to elements such as `d_max_cap[1, :A]`.
    containerize(A::Matrix) = Containers.DenseAxisArray(A, months, factories)
    ## Maximum production capacity in (month, factory) [units/month]:
    d_max_cap = containerize(
        [
            100000 50000
            110000 55000
            120000 60000
            145000 100000
            160000 0
            140000 70000
            155000 60000
            200000 100000
            210000 100000
            197000 100000
            80000 120000
            150000 150000
        ],
    )
    ## Minimum production capacity in (month, factory) [units/month]:
    d_min_cap = containerize(
        [
            20000 20000
            20000 20000
            20000 20000
            20000 20000
            20000 0
            20000 20000
            20000 20000
            20000 20000
            20000 20000
            20000 20000
            20000 20000
            20000 20000
        ],
    )
    ## Variable cost of production in (month, factory) [$/unit]:
    d_var_cost = containerize([
        10 5
        11 4
        12 3
        9 5
        8 0
        8 6
        5 4
        7 6
        9 8
        10 11
        8 10
        8 12
    ])
    ## Fixed cost of production in (month, factory) # [$/month]:
    d_fixed_cost = containerize(
        [
            500 600
            500 600
            500 600
            500 600
            500 0
            500 600
            500 600
            500 600
            500 600
            500 600
            500 600
            500 600
        ],
    )
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
    ## The production cannot be more that maximum capacity.
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
            d_fixed_cost[m, f] * status[m, f] +
            d_var_cost[m, f] * production[m, f] for m in months, f in factories
        )
    )
    ## Optimize the problem
    optimize!(model)
    ## Check the solution!
    Test.@testset "Check the solution against known optimal" begin
        Test.@test termination_status(model) == OPTIMAL
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

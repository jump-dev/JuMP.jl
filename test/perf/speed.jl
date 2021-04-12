#  Copyright 2017, Iain Dunning, Joey Huchette, Miles Lubin, and contributors
#  This Source Code Form is subject to the terms of the Mozilla Public
#  License, v. 2.0. If a copy of the MPL was not distributed with this
#  file, You can obtain one at https://mozilla.org/MPL/2.0/.
#############################################################################
# JuMP
# An algebraic modelling language for Julia
# See https://github.com/jump-dev/JuMP.jl
#############################################################################
# speed.jl
#
# Benchmarks model building time to test for performance regressions.
# Based on the models benchmarked in the paper:
#   Lubin, M., & Dunning, I. (2015).
#   Computing in operations research using Julia.
#   INFORMS Journal on Computing, 27(2), 238-248.
#
#############################################################################

import BenchmarkTools: @benchmark, allocs
using Random
using JuMP

"""
    p_median(num_facility, num_customer, num_location)

Implements the "p-median" facility location problem. We try to locate N
facilities such that we minimize the distance any given customer has to travel
to reach their closest facility. In this simple instance we will operate
in a 1D world with L possible locations for facilities, and customers being
located at random locations along the number line from 1 to D.

We use anonymous variables to remove the cost of name generation from the
benchmark.
"""
function p_median(num_facilities, num_customers, num_locations)
    Random.seed!(10)
    customer_locations = [rand(1:num_locations) for _ in 1:num_customers]

    model = Model()
    has_facility = @variable(model, [1:num_locations], Bin)
    is_closest = @variable(model, [1:num_locations, 1:num_customers], Bin)

    @objective(
        model,
        Min,
        sum(
            abs(customer_locations[customer] - location) *
            is_closest[location, customer] for customer in 1:num_customers,
            location in 1:num_locations
        )
    )

    for customer in 1:num_customers
        # `location` can't be closest for `customer` if there is no facility.
        @constraint(
            model,
            [location in 1:num_locations],
            is_closest[location, customer] <= has_facility[location]
        )
        # One facility must be the closest for `customer`.
        @constraint(
            model,
            sum(
                is_closest[location, customer] for location in 1:num_locations
            ) == 1
        )
    end

    # Must place all facilities.
    @constraint(model, sum(has_facility) == num_facilities)
end

println("P-Median(100 facilities, 100 customers, 5000 locations) benchmark:")
result = @benchmark p_median(100, 100, 5000)
display(result)
println()

"""
    cont5(n)

Based on a linear-Quadratic control problem (cont5_2_1) from one of Hans
Mittleman's instance collections.

We use anonymous variables to remove the cost of name generation from the
benchmark.
"""
function cont5(n)
    m = n
    n1 = n - 1
    m1 = m - 1
    dx = 1 / n
    T = 1.58
    dt = T / m
    h2 = dx^2
    a = 0.001
    yt = [0.5 * (1 - (j * dx)^2) for j in 0:n]

    model = Model()
    y = @variable(model, [0:m, 0:n], lower_bound = 0, upper_bound = 1)
    u = @variable(model, [1:m], lower_bound = -1, upper_bound = 1)

    @objective(model, Min, y[0, 0])

    # PDE
    for i in 0:m1
        for j in 1:n1
            @constraint(
                model,
                h2 * (y[i+1, j] - y[i, j]) ==
                0.5 *
                dt *
                (
                    y[i, j-1] - 2 * y[i, j] + y[i, j+1] + y[i+1, j-1] -
                    2 * y[i+1, j] + y[i+1, j+1]
                )
            )
        end
    end

    # Initial conditions.
    for j in 0:n
        @constraint(model, y[0, j] == 0)
    end

    # Boundary conditions.
    for i in 1:m
        @constraint(model, y[i, 2] - 4 * y[i, 1] + 3 * y[i, 0] == 0)
        @constraint(
            model,
            y[i, n-2] - 4 * y[i, n1] + 3 * y[i, n] ==
            (2 * dx) * (u[i] - y[i, n])
        )
    end
end

println("Cont5(n=500) benchmark:")
result = @benchmark cont5(500)
display(result)
println()

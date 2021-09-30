#############################################################################
# JuMP
# An algebraic modelling language for Julia
# See https://github.com/jump-dev/JuMP.jl
#############################################################################
# speed2.jl
#
# Runs some JuMP benchmarks to test for speed-related regressions.
# Test time to build model in memory.
#############################################################################

include("../solvers.jl")

using JuMP
using Random

function pMedian(solver, numFacility::Int, numCustomer::Int, numLocation::Int)
    Random.seed!(10)
    customerLocations = [rand(1:numLocation) for a in 1:numCustomer]

    buildTime = @elasped begin
        m = Model(solver = solver)

        # Facility locations
        @variable(m, 0 <= s[1:numLocation] <= 1)

        # Aux. variable: x_a,i = 1 iff the closest facility to a is at i
        @variable(m, 0 <= x[1:numLocation, 1:numCustomer] <= 1)

        # Objective: min distance
        @objective(
            m,
            Max,
            sum(
                abs(customerLocations[a] - i) * x[i, a] for
                a in 1:numCustomer, i in 1:numLocation
            )
        )

        # Constraints
        for a in 1:numCustomer
            # Subject to linking x with s
            for i in 1:numLocation
                @constraint(m, x[i, a] - s[i] <= 0)
            end
            # Subject to one of x must be 1
            @constraint(m, sum(x[i, a] for i in 1:numLocation) == 1)
        end

        # Subject to must allocate all facilities
        @constraint(m, sum(s[i] for i in 1:numLocation) == numFacility)
    end

    writeTime = @elapsed JuMP.build(m)

    return buildTime, writeTime
end

function cont5(solver, n)
    m = n
    n1 = n - 1
    m1 = m - 1
    dx = 1 / n
    T = 1.58
    dt = T / m
    h2 = dx^2
    a = 0.001
    yt = [0.5 * (1 - (j * dx)^2) for j in 0:n]

    buildTime = @elapsed begin
        mod = Model(solver = solver)
        @variable(mod, 0 <= y[0:m, 0:n] <= 1)
        @variable(mod, -1 <= u[1:m] <= 1)
        @objective(
            mod,
            Min,
            0.25 *
            dx *
            (
                (y[m, 0] - yt[1])^2 +
                2 * sum((y[m, j] - yt[j+1])^2 for j in 1:n1) +
                (y[m, n] - yt[n+1])^2
            ) + 0.25 * a * dt * (2 * sum(u[i]^2 for i in 1:m1) + u[m]^2)
        )

        # PDE
        for i in 0:m1
            for j in 1:n1
                @constraint(
                    mod,
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

        # IC
        for j in 0:n
            @constraint(mod, y[0, j] == 0)
        end

        # BC
        for i in 1:m
            @constraint(mod, y[i, 2] - 4 * y[i, 1] + 3 * y[i, 0] == 0)
            @constraint(
                mod,
                y[i, n-2] - 4 * y[i, n1] + 3 * y[i, n] ==
                (2 * dx) * (u[i] - y[i, n])
            )
        end
    end

    writeTime = @elapsed JuMP.build(mod)

    return buildTime, writeTime
end

function RunTests()
    # Pmedian
    if !isempty(lp_solvers)
        pmedian_build = Float64[]
        pmedian_write = Float64[]
        for runs in 1:9
            bt, wt = pMedian(first(lp_solvers), 100, 100, 5000)
            push!(pmedian_build, bt)
            push!(pmedian_write, wt)
        end
        sort!(pmedian_build)
        sort!(pmedian_write)
        print(
            "PMEDIAN BUILD MIN=",
            minimum(pmedian_build),
            "  MED=",
            pmedian_build[5],
            "\n",
        )
        print(
            "PMEDIAN INTRN MIN=",
            minimum(pmedian_write),
            "  MED=",
            pmedian_write[5],
            "\n",
        )
    else
        @warn("PMEDIAN NOT RUN!")
    end

    # Cont5
    if !isempty(quad_solvers)
        cont5_build = Float64[]
        cont5_write = Float64[]
        for runs in 1:9
            bt, wt = cont5(first(quad_solvers), 500)
            push!(cont5_build, bt)
            push!(cont5_write, wt)
        end
        sort!(cont5_build)
        sort!(cont5_write)
        print(
            "CONT5 BUILD   MIN=",
            minimum(cont5_build),
            "  MED=",
            cont5_build[5],
            "\n",
        )
        print(
            "CONT5 INTRN   MIN=",
            minimum(cont5_write),
            "  MED=",
            cont5_write[5],
            "\n",
        )
    else
        @warn("CONT5 NOT RUN!")
    end
end

RunTests()

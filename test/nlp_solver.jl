#  Copyright 2017, Iain Dunning, Joey Huchette, Miles Lubin, and contributors
#  This Source Code Form is subject to the terms of the Mozilla Public
#  License, v. 2.0. If a copy of the MPL was not distributed with this
#  file, You can obtain one at http://mozilla.org/MPL/2.0/.
#############################################################################
# JuMP
# An algebraic modelling language for Julia
# See http://github.com/JuliaOpt/JuMP.jl
#############################################################################
# Hock-Schittkowski Nonlinear Test Suite
# These tests include JuMP implementations of the models described in
#  W. Hock, K. Schittkowski, Test Examples for Nonlinear Programming
#  Codes, Lecture Notes in Economics and Mathematical Systems,
#  Springer, No, 187, 1981
# More information, including original model descriptions, is available at
# http://www.ai7.uni-bayreuth.de/downloads.htm
#############################################################################

# These tests are not executed as part of runtests because they depend on an
# NLP solver.
# TODO: Move these tests out of JuMP and into MOI so that JuMP can be tested
# separately from the solver, as it is for everything except NLP.

using Ipopt, JuMP
using Test
using MathOptInterface
const MOI = MathOptInterface

@testset "NLP solver tests" begin

    @testset "HS071" begin
        # hs071
        # Polynomial objective and constraints
        # min x1 * x4 * (x1 + x2 + x3) + x3
        # st  x1 * x2 * x3 * x4 >= 25
        #     x1^2 + x2^2 + x3^2 + x4^2 = 40
        #     1 <= x1, x2, x3, x4 <= 5
        # Start at (1,5,5,1)
        # End at (1.000..., 4.743..., 3.821..., 1.379...)
        m = Model(Ipopt.Optimizer)
        initval = [1,5,5,1]
        @variable(m, 1 <= x[i=1:4] <= 5, start=initval[i])
        @NLobjective(m, Min, x[1]*x[4]*(x[1]+x[2]+x[3]) + x[3])
        @NLconstraint(m, x[1]*x[2]*x[3]*x[4] >= 25)
        @NLconstraint(m, sum(x[i]^2 for i=1:4) == 40)

        set_silent(m)
        JuMP.optimize!(m)

        @test JuMP.has_values(m)
        @test JuMP.termination_status(m) == MOI.LOCALLY_SOLVED
        @test JuMP.primal_status(m) == MOI.FEASIBLE_POINT

        @test JuMP.value.(x) ≈ [1.000000, 4.742999, 3.821150, 1.379408] atol=1e-3
    end

    @testset "HS071 (no macros)" begin
        # hs071
        # Polynomial objective and constraints
        # min x1 * x4 * (x1 + x2 + x3) + x3
        # st  x1 * x2 * x3 * x4 >= 25
        #     x1^2 + x2^2 + x3^2 + x4^2 = 40
        #     1 <= x1, x2, x3, x4 <= 5
        # Start at (1,5,5,1)
        # End at (1.000..., 4.743..., 3.821..., 1.379...)
        m = Model(Ipopt.Optimizer)
        initval = [1,5,5,1]
        @variable(m, 1 <= x[i=1:4] <= 5, start=initval[i])
        JuMP.set_NL_objective(m, MOI.MIN_SENSE,
                              :($(x[1]) * $(x[4]) *
                                ($(x[1]) + $(x[2]) + $(x[3])) + $(x[3])))
        JuMP.add_NL_constraint(m,
                               :($(x[1]) * $(x[2]) * $(x[3]) * $(x[4]) >= 25))
        JuMP.add_NL_constraint(m, :($(x[1])^2 + $(x[2])^2 + $(x[3])^2 +
                                    $(x[4])^2 == 40))
        bad_expr = :(x[1]^2 + x[2]^2 + x[3]^2 + x[4]^2 == 40)
        @test_throws ErrorException JuMP.add_NL_constraint(m, bad_expr)

        set_silent(m)
        JuMP.optimize!(m)

        @test JuMP.has_values(m)
        @test JuMP.termination_status(m) == MOI.LOCALLY_SOLVED
        @test JuMP.primal_status(m) == MOI.FEASIBLE_POINT

        @test JuMP.value.(x) ≈ [1.000000, 4.742999, 3.821150, 1.379408] atol=1e-3
    end

    @testset "HS071 (epigraph)" begin
        # hs071, with epigraph formulation
        # Linear objective, nonlinear constraints
        # min t
        # st  t >= x1 * x4 * (x1 + x2 + x3) + x3
        #     ...
        model = Model(Ipopt.Optimizer)
        start = [1.0, 5.0, 5.0, 1.0]
        @variable(model, 1 <= x[i=1:4] <= 5, start = start[i])
        @variable(model, t, start = 100)
        @objective(model, Min, t)
        @NLconstraint(model, t >= x[1] * x[4] * (x[1] + x[2] + x[3]) + x[3])
        @NLconstraint(model, x[1] * x[2] * x[3] * x[4] >= 25)
        @NLconstraint(model, sum(x[i]^2 for i = 1:4) == 40)

        set_silent(model)
        JuMP.optimize!(model)
        @test JuMP.termination_status(model) == MOI.LOCALLY_SOLVED
        @test JuMP.primal_status(model) == MOI.FEASIBLE_POINT
        @test JuMP.value.(x) ≈
                [1.000000, 4.742999, 3.821150, 1.379408] atol=1e-3
    end

    @testset "HS109" begin
        a  = 50.176
        b1 = 0.25
        b = sin(b1)
        c = cos(b1)

        L = [0.0, 0.0, -0.55, -0.55, 196, 196, 196, -400, -400]
        U = [Inf, Inf,  0.55,  0.55, 252, 252, 252,  800,  800]

        m = Model(Ipopt.Optimizer)
        @variable(m, L[i] <= x[i=1:9] <= U[i], start = 0.0)

        @NLobjective(m, Min, 3 * x[1] + 1e-6 * x[1]^3 + 2 * x[2] + .522074e-6 * x[2]^3)

        @NLconstraint(m, x[4] - x[3] + 0.55 >= 0)
        @NLconstraint(m, x[3] - x[4] + 0.55 >= 0)
        @NLconstraint(m, 2250000 - x[1]^2 - x[8]^2 >= 0)
        @NLconstraint(m, 2250000 - x[2]^2 - x[9]^2 >= 0)
        @NLconstraint(m,
            x[5] * x[6] * sin(-x[3] - .25) + x[5] * x[7] * sin(-x[4] - .25) +
            2 * b * x[5]^2 - a * x[1] + 400 * a == 0)
        @NLconstraint(m,
            x[5] * x[6] * sin(x[3] - .25) + x[6] * x[7] * sin(x[3] - x[4] - .25) +
            2 * b * x[6]^2 - a * x[2] + 400 * a == 0)
        @NLconstraint(m,
            x[5] * x[7] * sin(x[4] - .25) + x[6] * x[7] * sin(x[4] - x[3] - .25) +
            2 * b * x[7]^2 + 881.779 * a == 0)
        @NLconstraint(m,
            a * x[8] + x[5] * x[6] * cos(-x[3] - .25) +
            x[5] * x[7] * cos(-x[4] - .25) - 200 * a - 2 * c * x[5]^2 +
            .7533e-3 * a * x[5]^2 == 0)
        @NLconstraint(m,
            a * x[9] + x[5] * x[6] * cos(x[3] - .25) +
            x[6] * x[7] * cos(x[3] - x[4] - .25) - 2 * c * x[6]^2 +
            .7533e-3 * a * x[6]^2 - 200 * a == 0)
        @NLconstraint(m,
            x[5] * x[7] * cos(x[4] - .25) + x[6] * x[7] * cos(x[4] - x[3] - .25) -
            2 * c * x[7]^2 + 22.938 * a + .7533e-3 * a * x[7]^2 == 0)

        set_silent(m)
        JuMP.optimize!(m)

        @test JuMP.has_values(m)
        @test JuMP.termination_status(m) == MOI.LOCALLY_SOLVED
        @test JuMP.primal_status(m) == MOI.FEASIBLE_POINT

        @test JuMP.objective_value(m) ≈ 5326.851310161077 atol=1e-5
    end

    @testset "HS110" begin
        m = Model(Ipopt.Optimizer)
        @variable(m, -2.001 <= x[1:10] <= 9.999, start = 9)

        @NLobjective(m, Min,
            sum(log(x[j] - 2)^2 + log(10 - x[j])^2 for j=1:10) -
            prod(x[i] for i=1:10)^0.2
        )

        set_silent(m)
        JuMP.optimize!(m)

        @test JuMP.has_values(m)
        # Ipopt returns AlmostLOCALLY_SOLVED and NearlyFEASIBLE_POINT on this instance.
        # @test JuMP.termination_status(m) == MOI.LOCALLY_SOLVED
        # @test JuMP.primal_status(m) == MOI.FEASIBLE_POINT

        @test JuMP.objective_value(m) ≈ -45.77846971 atol=1e-5
    end

    @testset "HS111" begin
        c = [-6.089, -17.164, -34.054, -5.914, -24.721, -14.986, -24.100, -10.708, -26.662, -22.179]

        m = Model(Ipopt.Optimizer)
        @variable(m, -100 <= x[1:10] <= 100, start = -2.3)

        @NLobjective(m, Min,
            sum(exp(x[j]) * (c[j] + x[j] - log( sum(exp(x[k]) for k=1:10) ) ) for j=1:10))

        @NLconstraint(m, exp(x[1]) + 2*exp(x[2]) + 2*exp(x[3]) +   exp(x[6]) + exp(x[10]) == 2)
        @NLconstraint(m, exp(x[4]) + 2*exp(x[5]) +   exp(x[6]) +   exp(x[7])              == 1)
        @NLconstraint(m, exp(x[3]) +   exp(x[7]) +   exp(x[8]) + 2*exp(x[9]) + exp(x[10]) == 1)

        set_silent(m)
        JuMP.optimize!(m)

        @test JuMP.has_values(m)
        @test JuMP.termination_status(m) == MOI.LOCALLY_SOLVED
        @test JuMP.primal_status(m) == MOI.FEASIBLE_POINT

        @test JuMP.objective_value(m) ≈ -47.76109026 atol=1e-5
    end

    @testset "HS112" begin
        c = [-6.089, -17.164, -34.054, -5.914, -24.721, -14.986, -24.100, -10.708, -26.662, -22.179]

        m = Model(Ipopt.Optimizer)
        @variable(m, x[1:10] >= 1e-6, start = 0.1)

        @NLobjective(m, Min, sum(x[j]*(c[j] + log(x[j]/sum(x[k] for k=1:10))) for j=1:10))

        @NLconstraint(m, x[1] + 2*x[2] + 2*x[3] + x[6] + x[10] == 2)
        @NLconstraint(m, x[4] + 2*x[5] + x[6] + x[7] == 1)
        @NLconstraint(m, x[3] + x[7] + x[8] + 2*x[9] + x[10] == 1)

        set_silent(m)
        JuMP.optimize!(m)

        @test JuMP.has_values(m)
        @test JuMP.termination_status(m) == MOI.LOCALLY_SOLVED
        @test JuMP.primal_status(m) == MOI.FEASIBLE_POINT

        @test JuMP.objective_value(m) ≈ -47.76109026 atol=1e-5
    end

    @testset "HS114" begin
        n = 10
        a = 0.99
        b = 0.9

        lower = [0.0001, 0.0001, 0.0001, 0.0001, 0.0001, 85, 90, 3, 1.2, 145]
        upper = [2000, 16000, 120, 5000, 2000, 93, 95, 12, 4, 162]
        start = [1745, 12000, 110, 3048, 1974, 89.2, 92.8, 8, 3.6, 145]

        m = Model(Ipopt.Optimizer)
        @variable(m, lower[i] <= x[i=1:n] <= upper[i], start = start[i])

        @NLobjective(m, Min, 5.04*x[1] + .035*x[2] + 10*x[3] + 3.36*x[5] - .063*x[4]*x[7])

        @NLconstraint(m, 35.82 - .222*x[10] - b*x[9] >= 0)
        @NLconstraint(m, -133 + 3*x[7] - a*x[10] >= 0)
        @NLconstraint(m, -(35.82 - .222*x[10] - b*x[9]) + x[9]*(1/b - b) >= 0)
        @NLconstraint(m, -(-133 + 3*x[7] - a*x[10]) + (1/a - a)*x[10] >= 0)
        @NLconstraint(m, 1.12*x[1] + .13167*x[1]*x[8] - .00667*x[1]*x[8]^2 - a*x[4] >= 0)
        @NLconstraint(m, 57.425 + 1.098*x[8] - .038*x[8]^2 + .325*x[6] - a*x[7] >= 0)
        @NLconstraint(m, -(1.12*x[1] + .13167*x[1]*x[8] - .00667*x[1]*x[8]^2 - a*x[4]) + (1/a - a)*x[4] >= 0)
        @NLconstraint(m, -(57.425 + 1.098*x[8] - .038*x[8]^2 + .325*x[6] - a*x[7]) + (1/a - a)*x[7] >= 0)
        @NLconstraint(m, 1.22*x[4] - x[1] - x[5] == 0)
        @NLconstraint(m, 98000*x[3]/(x[4]*x[9] + 1000*x[3]) - x[6] == 0)
        @NLconstraint(m, (x[2] + x[5])/x[1] - x[8] == 0)

        set_silent(m)
        JuMP.optimize!(m)

        @test JuMP.has_values(m)
        @test JuMP.termination_status(m) == MOI.LOCALLY_SOLVED
        @test JuMP.primal_status(m) == MOI.FEASIBLE_POINT

        @test JuMP.objective_value(m) ≈ -1768.80696 atol=1e-3
    end

    @testset "HS116" begin
        N = 13
        a = 0.002
        b = 1.262626
        c = 1.231059
        d = 0.03475
        e = 0.975
        f = 0.00975

        lower = [0.1, 0.1, 0.1, 0.0001, 0.1, 0.1, 0.1, 0.1, 500, 0.1, 1.0, 0.0001, 0.0001, 0.0, 0.0, 0.0]
        upper = [1.0, 1.0, 1.0, 0.1, 0.9, 0.9, 1000, 1000, 1000, 500, 150, 150, 150, Inf, Inf, Inf]
        start = [0.5  2 0.8  3 0.9  4 0.1  5 0.14  6 0.5  7 489  8 80  9 650 0.5  2 0.8  3 0.9  4 0.1  5 0.14  6 0.5  7 489  8 80  9 650]

        m = Model(Ipopt.Optimizer)
        @variable(m, lower[i] <= x[i=1:N] <= upper[i], start = start[i])
        @NLobjective(m, Min, x[11] + x[12] + x[13])

        @NLconstraints m begin
            x[3] - x[2] >= 0
            x[2] - x[1] >= 0
            1 - a * x[7] + a * x[8] >= 0
            x[11] + x[12] + x[13] >= 50
            x[13] - b * x[10] + c * x[3] * x[10] >= 0
            x[5] - d * x[2] - e * x[2] * x[5] + f * x[2]^2 >= 0
            x[6] - d * x[3] - e * x[3] * x[6] + f * x[3]^2 >= 0
            x[4] - d * x[1] - e * x[1] * x[4] + f * x[1]^2 >= 0
            x[12] - b * x[9] + c * x[2] * x[9] >= 0
            x[11] - b * x[8] + c * x[1] * x[8] >= 0
            x[5] * x[7] - x[1] * x[8] - x[4] * x[7] + x[4] * x[8] >= 0
            1 - a * (x[2] * x[9] + x[5] * x[8] - x[1] * x[8] - x[6] * x[9]) -
                  x[5] - x[6] >= 0
            x[2] * x[9] - x[3] * x[10] - x[6] * x[9] - 500 * x[2] +
                  500 * x[6] + x[2] * x[10] >= 0
            x[2] - 0.9 - a * (x[2] * x[10] - x[3] * x[10]) >= 0
            x[11] + x[12] + x[13] <= 250
        end

        set_silent(m)
        JuMP.optimize!(m)

        @test JuMP.has_values(m)
        @test JuMP.termination_status(m) == MOI.LOCALLY_SOLVED
        @test JuMP.primal_status(m) == MOI.FEASIBLE_POINT

        # This test occasionally fails, for unknown reasons.
        @test JuMP.objective_value(m) ≈ 97.588409 atol=1e-3
    end

    @testset "HS118" begin
        m = Model(Ipopt.Optimizer)

        L = zeros(15)
        L[1] =  8.0
        L[2] = 43.0
        L[3] =  3.0

        U = zeros(15)
        U[1] = 21.0
        U[2] = 57.0
        U[3] = 16.0
        for k in 1:4
            U[3*k+1] =  90.0
            U[3*k+2] = 120.0
            U[3*k+3] =  60.0
        end

        start = [20.0, 55.0, 15.0, 20.0, 60.0, 20.0, 20.0, 60.0, 20.0, 20.0,
                 60.0, 20.0, 20.0, 60.0, 20.0]

        @variable(m, L[i] <= x[i=1:15] <= U[i])

        # Initial solution (could also use 'start' keyword in @variable)
        JuMP.set_start_value.(x, start)

        @NLobjective(m, Min,
            sum(2.3     * x[3*k+1]   +
                0.0001  * x[3*k+1]^2 +
                1.7     * x[3*k+2]   +
                0.0001  * x[3*k+2]^2 +
                2.2     * x[3*k+3] +
                0.00015 * x[3*k+3]^2 for k=0:4))

        # constr1
        for j in 1:4
            @constraint(m, x[3*j+1] - x[3*j-2] + 7 <= 13)
            @constraint(m, x[3*j+1] - x[3*j-2] + 7 >=  0)
        end

        # constr2
        for j in 1:4
            @constraint(m, x[3*j+2] - x[3*j-1] + 7 <= 14)
            @constraint(m, x[3*j+2] - x[3*j-1] + 7 >=  0)
        end

        # constr3
        for j in 1:4
            @constraint(m, x[3*j+3] - x[3*j  ] + 7 <= 13)
            @constraint(m, x[3*j+3] - x[3*j  ] + 7 >=  0)
        end

        @constraint(m, x[1] + x[2] + x[3]    >= 60)
        @constraint(m, x[4] + x[5] + x[6]    >= 50)
        @constraint(m, x[7] + x[8] + x[9]    >= 70)
        @constraint(m, x[10] + x[11] + x[12] >= 85)
        @constraint(m, x[13] + x[14] + x[15] >= 100)

        set_silent(m)
        JuMP.optimize!(m)

        @test JuMP.has_values(m)
        @test JuMP.termination_status(m) == MOI.LOCALLY_SOLVED
        @test JuMP.primal_status(m) == MOI.FEASIBLE_POINT

        @test JuMP.value.(x[1:4]) ≈ [8.0, 49.0, 3.0, 1.0] atol=1e-4
        @test JuMP.objective_value(m) ≈ 664.82045 atol=1e-5
    end

    @testset "Two-sided constraints" begin
        m = Model(Ipopt.Optimizer)
        @variable(m, x)
        @NLobjective(m, Max, x)
        l = -1
        u = 1
        @NLconstraint(m, l <= x <= u)
        JuMP.optimize!(m)

        @test JuMP.has_values(m)
        @test JuMP.termination_status(m) == MOI.LOCALLY_SOLVED
        @test JuMP.primal_status(m) == MOI.FEASIBLE_POINT
        @test JuMP.objective_value(m) ≈ u atol=1e-6

        @NLobjective(m, Min, x)
        set_silent(m)
        JuMP.optimize!(m)

        @test JuMP.has_values(m)
        @test JuMP.termination_status(m) == MOI.LOCALLY_SOLVED
        @test JuMP.primal_status(m) == MOI.FEASIBLE_POINT
        @test JuMP.objective_value(m) ≈ l atol=1e-6
    end

    @testset "Two-sided constraints (no macros)" begin
        m = Model(Ipopt.Optimizer)
        @variable(m, x)
        JuMP.set_NL_objective(m, MOI.MAX_SENSE, x)
        l = -1
        u = 1
        JuMP.add_NL_constraint(m, :($l <= $x <= $u))
        set_silent(m)
        JuMP.optimize!(m)

        @test JuMP.has_values(m)
        @test JuMP.termination_status(m) == MOI.LOCALLY_SOLVED
        @test JuMP.primal_status(m) == MOI.FEASIBLE_POINT
        @test JuMP.objective_value(m) ≈ u atol=1e-6

        JuMP.set_NL_objective(m, MOI.MIN_SENSE, x)
        JuMP.optimize!(m)

        @test JuMP.has_values(m)
        @test JuMP.termination_status(m) == MOI.LOCALLY_SOLVED
        @test JuMP.primal_status(m) == MOI.FEASIBLE_POINT
        @test JuMP.objective_value(m) ≈ l atol=1e-6
    end

    @testset "Duals" begin
        m = Model(Ipopt.Optimizer)
        @variable(m, x >= 0)
        @variable(m, y <= 5)
        @variable(m, 2 <= z <= 4)
        @variable(m, 0 <= r[i=3:6] <= i)
        @NLobjective(m, Min, -((x + y)/2.0 + 3.0)/3.0 - z - r[3])
        @constraint(m, cons1, x+y >= 2)
        @constraint(m, cons2, sum(r[i] for i=3:5) <= (2 - x)/2.0)
        @NLconstraint(m, cons3, 7.0*y <= z + r[6]/1.9)

        function test_result()
            @test JuMP.has_values(m)
            @test JuMP.termination_status(m) == MOI.LOCALLY_SOLVED
            @test JuMP.primal_status(m) == MOI.FEASIBLE_POINT


            @test JuMP.value(x) ≈ 0.9774436 atol=1e-6
            @test JuMP.value(y) ≈ 1.0225563 atol=1e-6
            @test JuMP.value(z) ≈ 4.0 atol=1e-6
            @test JuMP.value(r[3]) ≈ 0.5112781 atol=1e-6
            @test JuMP.value(r[4]) ≈ 0.0 atol=1e-6
            @test JuMP.value(r[5]) ≈ 0.0 atol=1e-6
            @test JuMP.value(r[6]) ≈ 6.0 atol=1e-6
            @test JuMP.dual_status(m) == MOI.FEASIBLE_POINT
            # Reduced costs
            @test JuMP.dual(JuMP.LowerBoundRef(x)) ≈ 0.0 atol=1e-6
            @test JuMP.dual(JuMP.UpperBoundRef(y)) ≈ 0.0 atol=1e-6
            @test JuMP.dual(JuMP.UpperBoundRef(z)) ≈ -1.0714286 atol=1e-6
            @test JuMP.dual(JuMP.LowerBoundRef(r[3])) ≈ 0.0 atol=1e-6
            @test JuMP.dual(JuMP.UpperBoundRef(r[3])) ≈ 0.0 atol=1e-6
            @test JuMP.dual(JuMP.LowerBoundRef(r[4])) ≈ 1.0 atol=1e-6
            @test JuMP.dual(JuMP.UpperBoundRef(r[4])) ≈ 0.0 atol=1e-6
            @test JuMP.dual(JuMP.LowerBoundRef(r[5])) ≈ 1.0 atol=1e-6
            @test JuMP.dual(JuMP.UpperBoundRef(r[5])) ≈ 0.0 atol=1e-6
            @test JuMP.dual(JuMP.UpperBoundRef(r[6])) ≈ -0.03759398 atol=1e-6
            @test JuMP.dual(JuMP.LowerBoundRef(r[6])) ≈ 0.0 atol=1e-6

            # Constraint duals
            @test JuMP.dual(cons1) ≈ 0.333333 atol=1e-6
            @test JuMP.dual(cons2) ≈ -1.0 atol=1e-6
            @test JuMP.dual(cons3) ≈ -0.0714286 atol=1e-6
        end

        set_silent(m)
        JuMP.optimize!(m)
        test_result()
        @test JuMP.objective_value(m) ≈ -5.8446115 atol=1e-6

        # Same objective with sense/sign flipped.
        @NLobjective(m, Max, ((x + y)/2.0 + 3.0)/3.0 + z + r[3])

        JuMP.optimize!(m)
        test_result()
        @test JuMP.objective_value(m) ≈ 5.8446115 atol=1e-6
    end

    @testset "Quadratic inequality constraints, linear objective" begin
        m = Model(Ipopt.Optimizer)
        @variable(m, -2 <= x <= 2)
        @variable(m, -2 <= y <= 2)
        @objective(m, Min, x - y)
        @constraint(m, x + x^2 + x*y + y^2 <= 1)
        set_silent(m)
        JuMP.optimize!(m)

        @test JuMP.has_values(m)
        @test JuMP.termination_status(m) == MOI.LOCALLY_SOLVED
        @test JuMP.primal_status(m) == MOI.FEASIBLE_POINT
        @test JuMP.objective_value(m) ≈ -1-4/sqrt(3) atol=1e-6
        @test JuMP.value(x) + JuMP.value(y) ≈ -1/3 atol=1e-3
    end

    @testset "Quadratic inequality constraints, NL objective" begin
        m = Model(Ipopt.Optimizer)
        @variable(m, -2 <= x <= 2)
        @variable(m, -2 <= y <= 2)
        @NLobjective(m, Min, x - y)
        @constraint(m, x + x^2 + x*y + y^2 <= 1)
        set_silent(m)
        JuMP.optimize!(m)

        @test JuMP.has_values(m)
        @test JuMP.termination_status(m) == MOI.LOCALLY_SOLVED
        @test JuMP.primal_status(m) == MOI.FEASIBLE_POINT
        @test JuMP.objective_value(m) ≈ -1-4/sqrt(3) atol=1e-6
        @test JuMP.value(x) + JuMP.value(y) ≈ -1/3 atol=1e-3
    end

    @testset "Quadratic equality constraints" begin
        m = Model(Ipopt.Optimizer)
        @variable(m, 0 <= x[1:2] <= 1)
        @constraint(m, x[1]^2 + x[2]^2 == 1/2)
        @NLobjective(m, Max, x[1] - x[2])
        set_silent(m)
        JuMP.optimize!(m)

        @test JuMP.has_values(m)
        @test JuMP.termination_status(m) == MOI.LOCALLY_SOLVED
        @test JuMP.primal_status(m) == MOI.FEASIBLE_POINT
        @test JuMP.objective_value(m) ≈ sqrt(1/2) atol=1e-6
        @test JuMP.value.(x) ≈ [sqrt(1/2), 0] atol=1e-6
    end

    @testset "Fixed variables" begin
        m = Model(Ipopt.Optimizer)
        set_silent(m)
        @variable(m, x == 0)
        @variable(m, y ≥ 0)
        @objective(m, Min, y)
        @NLconstraint(m, y ≥ x^2)
        for α in 1:4
            JuMP.fix(x, α)
            JuMP.optimize!(m)
            @test JuMP.value(y) ≈ α^2 atol=1e-6
        end
    end

    @testset "ifelse" begin
        model = Model(Ipopt.Optimizer)
        @variable(model, x, start = 2)
        # The minimizer is at smooth point, so solvers should be okay.
        @NLobjective(model, Min, ifelse( x <= 1, x^2, x) )
        set_silent(model)
        JuMP.optimize!(model)

        @test JuMP.termination_status(model) == MOI.LOCALLY_SOLVED
        @test JuMP.primal_status(model) == MOI.FEASIBLE_POINT
        @test JuMP.value(x) ≈ 0.0 atol=1e-5
    end

    @testset "infeasible problem" begin
        # (Attempt to) solve an infeasible problem
        model = Model(Ipopt.Optimizer)
        n = 10
        @variable(model, 0 <= x[i=1:n] <= 1)
        @NLobjective(model, Max, x[n])
        for i in 1:n-1
            @NLconstraint(model, x[i+1]-x[i] == 0.15)
        end
        set_silent(model)
        JuMP.optimize!(model)

        @test JuMP.termination_status(model) == MOI.LOCALLY_INFEASIBLE
    end

    @testset "unbounded problem" begin
        model = Model(Ipopt.Optimizer)
        @variable(model, x >= 0)
        @NLobjective(model, Max, x)
        @NLconstraint(model, x >= 5)
        set_silent(model)
        JuMP.optimize!(model)
        @test JuMP.termination_status(model) == MOI.NORM_LIMIT
    end

    @testset "Derivatives of x^4, x < 0" begin
        model = Model(Ipopt.Optimizer)
        @variable(model, x >= -1, start = -0.5)
        @NLobjective(model, Min, x^4)
        set_silent(model)
        JuMP.optimize!(model)

        @test JuMP.termination_status(model) == MOI.LOCALLY_SOLVED
        @test JuMP.primal_status(model) == MOI.FEASIBLE_POINT
        @test JuMP.value(x) ≈ 0.0 atol=1e-2
    end

    # This test seems to be checking for a bug in expression handling.
    @testset "Entropy maximization" begin
        model = Model(Ipopt.Optimizer)
        @variable(model, x[1:4] >= 0, start = 1)
        @variable(model, z[1:4], start = 0)
        @NLexpression(model, entropy[i=1:4], -x[i] * log(x[i]))
        @NLobjective(model, Max, sum(z[i] for i = 1:2) +
                                 sum(z[i]/2 for i = 3:4))
        @NLconstraint(model, z_constr1[i=1], z[i] <= entropy[i])
        @NLconstraint(model, z_constr1_dup[i=2], z[i] <= entropy[i]) # duplicate expressions
        @NLconstraint(model, z_constr2[i=3:4], z[i] <= 2 * entropy[i])
        @constraint(model, sum(x) == 1)
        set_silent(model)
        JuMP.optimize!(model)

        @test JuMP.termination_status(model) == MOI.LOCALLY_SOLVED
        @test JuMP.primal_status(model) == MOI.FEASIBLE_POINT
        @test JuMP.value.(x) ≈ [1/4, 1/4, 1/4, 1/4] atol=1e-3
        @test JuMP.value(entropy[1]) ≈ -(1/4) * log(1/4) atol=1e-4
    end

    @testset "Changing objectives" begin
        model = Model(Ipopt.Optimizer)
        @variable(model, x >= 0)
        @variable(model, y >= 0)
        @objective(model, Max, x + y)
        @NLconstraint(model, x + 2y <= 1)
        JuMP.optimize!(model)
        @test JuMP.value(x) ≈ 1.0 atol=1e-4
        @test JuMP.value(y) ≈ 0.0 atol=1e-4
        @test JuMP.objective_value(model) ≈ 1.0 atol=1e-4

        @objective(model, Max, 2x + y)
        set_silent(model)
        JuMP.optimize!(model)
        @test JuMP.value(x) ≈ 1.0 atol=1e-4
        @test JuMP.value(y) ≈ 0.0 atol=1e-4
        @test JuMP.objective_value(model) ≈ 2.0 atol=1e-4
    end

    @testset "Setting NLobjective then objective" begin
        model = Model(Ipopt.Optimizer)
        @variable(model, x >= 0)
        @variable(model, y >= 0)
        @NLobjective(model, Max, x + y)
        @NLconstraint(model, x + 2y <= 1)
        JuMP.optimize!(model)
        @test JuMP.value(x) ≈ 1.0 atol=1e-4
        @test JuMP.value(y) ≈ 0.0 atol=1e-4
        @test JuMP.objective_value(model) ≈ 1.0 atol=1e-4

        @objective(model, Max, 2x + y)
        set_silent(model)
        JuMP.optimize!(model)
        @test JuMP.value(x) ≈ 1.0 atol=1e-4
        @test JuMP.value(y) ≈ 0.0 atol=1e-4
        @test JuMP.objective_value(model) ≈ 2.0 atol=1e-4
    end

    my_square(x) = x^2
    function my_f(x,y)
        return (x-1)^2+(y-2)^2
    end

    @testset "User-defined functions" begin
        model = Model(Ipopt.Optimizer)
        JuMP.register(model, :my_f, 2, my_f, autodiff=true)
        JuMP.register(model, :my_square, 1, my_square, autodiff=true)

        @variable(model, x[1:2] >= 0.5)
        @NLobjective(model, Min, my_f(x[1], my_square(x[2])))

        set_silent(model)
        JuMP.optimize!(model)
        @test JuMP.termination_status(model) == MOI.LOCALLY_SOLVED
        @test JuMP.primal_status(model) == MOI.FEASIBLE_POINT
        @test JuMP.value.(x) ≈ [1, sqrt(2.0)]
    end

    @testset "Univariate user-defined functions" begin
        model = Model(Ipopt.Optimizer)
        JuMP.register(model, :my_square, 1, my_square, autodiff=true)

        # Test just univariate functions because this is a path where hessians
        # are enabled.
        @variable(model, x[1:2] >= 0.5)
        @NLobjective(model, Min, my_square(x[1]-0.4) + my_square(x[2] - 2.0))
        set_silent(model)
        JuMP.optimize!(model)
        @test JuMP.termination_status(model) == MOI.LOCALLY_SOLVED
        @test JuMP.primal_status(model) == MOI.FEASIBLE_POINT
        @test JuMP.value.(x) ≈ [0.5, 2.0] atol=1e-4
    end

    @testset "Issue #927" begin
        model = Model(Ipopt.Optimizer)
        JuMP.register(model, :my_f, 2, my_f, autodiff=true)
        @variable(model, x)
        @NLobjective(model, Min, my_f(x, x))
        set_silent(model)
        JuMP.optimize!(model)
        @test JuMP.termination_status(model) == MOI.LOCALLY_SOLVED
        @test JuMP.primal_status(model) == MOI.FEASIBLE_POINT
        @test JuMP.value.(x) ≈ 1.5 atol=1e-4
    end
end

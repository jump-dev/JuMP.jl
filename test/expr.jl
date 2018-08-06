#  Copyright 2017, Iain Dunning, Joey Huchette, Miles Lubin, and contributors
#  This Source Code Form is subject to the terms of the Mozilla Public
#  License, v. 2.0. If a copy of the MPL was not distributed with this
#  file, You can obtain one at http://mozilla.org/MPL/2.0/.
#############################################################################
# JuMP
# An algebraic modeling langauge for Julia
# See http://github.com/JuliaOpt/JuMP.jl
#############################################################################
# test/expr.jl
# Testing for AffExpr and QuadExpr
#############################################################################
using JuMP, Compat.Test

@testset "Expressions" begin

    @testset "Test expression construction" begin
        maff = Model()
        @variable(maff, 0 <= x[1:5] <= 1)
        @variable(maff, 0 <= LongName <= 99)

        @testset "AffExpr" begin
            # Test affstr
            a1 = x[1] + LongName + 5
            @test string(a1) == "x[1] + LongName + 5"
            # Test like term collection
            a2 = 2*(x[2] + LongName + x[2]) + 0
            @test string(a2) == "4 x[2] + 2 LongName"
            # Test appending functionality
            push!(a1, 5.0, x[2])
            @test string(a1) == "x[1] + LongName + 5 x[2] + 5"
            append!(a1, a2)
            @test string(a1) == "x[1] + 3 LongName + 9 x[2] + 5"
            append!(a1, 2.0)
            @test string(a1) == "x[1] + 3 LongName + 9 x[2] + 7"
            append!(a1, LongName)
            @test string(a1) == "x[1] + 4 LongName + 9 x[2] + 7"
            append!(a1, 2)
            @test string(a1) == "x[1] + 4 LongName + 9 x[2] + 9"
        end

        @testset "QuadExpr" begin
            # Test string
            q1 = x[1]*x[2] + 27.2*LongName + 5
            @test string(q1) == "x[1]*x[2] + 27.2 LongName + 5"
            # Test like term collection
            q2 = x[1]*x[2] + x[2]*x[1]
            @test string(q2) == "2 x[1]*x[2]"
        end
    end

    @testset "Test getvalue(expr)" begin
        m = Model()
        @variable(m, 1 <= x[1:3] <= 2)
        setvalue(x[3], 2)
        setvalue(x[2], 2)
        setvalue(x[1], 1)
        @test isapprox(getvalue(x[1]-x[2]+2x[3]-1.0), 2.0)
        @test isapprox(getvalue(x[1]*x[1]-2x[2]*x[1]+3x[2]+1), 4.0)
    end

    @testset "Test expression iterators" begin
        m = Model()
        @variable(m, x[1:10])

        a1 = 1*x[1] + 2*x[2]
        k = 1
        for (coeff,var) in linearterms(a1)
            if k == 1
                @test coeff == 1
                @test var === x[1]
            elseif k == 2
                @test coeff == 2
                @test var === x[2]
            end
            k += 1
        end

        k = 0
        a2 = zero(AffExpr)
        for (coeff, var) in linearterms(a2)
            k += 1
        end
        @test k == 0
    end
end

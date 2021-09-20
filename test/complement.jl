#  Copyright 2017, Iain Dunning, Joey Huchette, Miles Lubin, and contributors
#  This Source Code Form is subject to the terms of the Mozilla Public
#  License, v. 2.0. If a copy of the MPL was not distributed with this
#  file, You can obtain one at https://mozilla.org/MPL/2.0/.

using JuMP
using Test

include(joinpath(@__DIR__, "utilities.jl"))

@testset "Scalar Complementarity" begin
    @testset "complements" begin
        model = Model()
        @variable(model, x >= 0)
        @constraint(model, c, complements(2x - 1, x))
        obj = constraint_object(c)
        @test obj.func == [2x - 1, x]
        @test obj.set == MOI.Complements(2)
    end

    @testset "⟂" begin
        model = Model()
        @variable(model, x >= 0)
        @variable(model, y >= 0)
        @constraint(model, c, x ⟂ y)
        obj = constraint_object(c)
        @test obj.func == [x, y]
        @test obj.set == MOI.Complements(2)
    end

    @testset "error: x--F" begin
        model = Model()
        @variable(model, x >= 0)
        @test_throws_strip(
            ErrorException(
                "In `@constraint(model, x ⟂ 2x - 1)`: second term must be a " *
                "variable.",
            ),
            @constraint(model, x ⟂ 2x - 1)
        )
    end

    @testset "error: F--F" begin
        model = Model()
        @variable(model, x >= 0)
        @test_throws_strip(
            ErrorException(
                "In `@constraint(model, x + 1 ⟂ 2x - 1)`: second term must " *
                "be a variable.",
            ),
            @constraint(model, x + 1 ⟂ 2x - 1)
        )
    end
end

@testset "Vector Complementarity" begin
    @testset "complements" begin
        model = Model()
        @variable(model, x[1:2] >= 0)
        @constraint(model, c, complements(2x .- 1, x))
        obj = constraint_object(c)
        @test obj.func == [2x .- 1; x]
        @test obj.set == MOI.Complements(4)
    end

    @testset "⟂" begin
        model = Model()
        @variable(model, x[1:2] >= 0)
        @variable(model, y[1:2] >= 0)
        @constraint(model, c, x ⟂ y)
        obj = constraint_object(c)
        @test obj.func == [x; y]
        @test obj.set == MOI.Complements(4)
    end

    @testset "error: length mismatch" begin
        model = Model()
        @variable(model, x[1:2] >= 0)
        @test_throws_strip(
            ErrorException(
                "In `@constraint(model, x ⟂ [x[1]])`: size of mapping does " *
                "not match size of variables: (2,) != (1,).",
            ),
            @constraint(model, x ⟂ [x[1]])
        )
    end

    @testset "error: x--F" begin
        model = Model()
        @variable(model, x[1:2] >= 0)
        @test_throws_strip(
            ErrorException(
                "In `@constraint(model, x ⟂ 2x .- 1)`: second term must be an " *
                "array of variables.",
            ),
            @constraint(model, x ⟂ 2x .- 1)
        )
    end

    @testset "error: F--F" begin
        model = Model()
        @variable(model, x[1:2] >= 0)
        @test_throws_strip(
            ErrorException(
                "In `@constraint(model, x .+ 1 ⟂ 2x .- 1)`: second term must " *
                "be an array of variables.",
            ),
            @constraint(model, x .+ 1 ⟂ 2x .- 1)
        )
    end
end

@testset "SparseAxisArray" begin
    @testset "complements" begin
        model = Model()
        @variable(model, x[i = 1:3; isodd(i)] >= 0)
        @constraint(model, c, complements(2x .- 1, x))
        obj = constraint_object(c)
        @test obj.func == [2x[3] - 1, 2x[1] - 1, x[3], x[1]] ||
              obj.func == [2x[1] - 1, 2x[3] - 1, x[1], x[3]]
        @test obj.set == MOI.Complements(4)
    end

    @testset "⟂" begin
        model = Model()
        @variable(model, x[i = 1:3; isodd(i)] >= 0)
        @constraint(model, c, 2x .- 1 ⟂ x)
        obj = constraint_object(c)
        @test obj.func == [2x[3] - 1, 2x[1] - 1, x[3], x[1]] ||
              obj.func == [2x[1] - 1, 2x[3] - 1, x[1], x[3]]
        @test obj.set == MOI.Complements(4)
    end

    @testset "key mismatch" begin
        model = Model()
        @variable(model, x[i = 1:3; isodd(i)] >= 0)
        @variable(model, y[i = 3:5; isodd(i)] >= 0)
        @test_throws(ErrorException, @constraint(model, 2x .- 1 ⟂ y))
    end
end

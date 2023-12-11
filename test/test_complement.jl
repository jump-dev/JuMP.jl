#  Copyright 2017, Iain Dunning, Joey Huchette, Miles Lubin, and contributors
#  This Source Code Form is subject to the terms of the Mozilla Public
#  License, v. 2.0. If a copy of the MPL was not distributed with this
#  file, You can obtain one at https://mozilla.org/MPL/2.0/.

module TestComplement

using JuMP
using Test

include(joinpath(@__DIR__, "utilities.jl"))

function test_scalar_complements()
    model = Model()
    @variable(model, x >= 0)
    @constraint(model, c, complements(2x - 1, x))
    obj = constraint_object(c)
    @test obj.func == [2x - 1, x]
    @test obj.set == MOI.Complements(2)
    return
end

function test_scalar_perp()
    model = Model()
    @variable(model, x >= 0)
    @variable(model, y >= 0)
    @constraint(model, c, x ⟂ y)
    obj = constraint_object(c)
    @test obj.func == [x, y]
    @test obj.set == MOI.Complements(2)
    return
end

function test_scalar_error_x_F()
    model = Model()
    @variable(model, x >= 0)
    @test_throws_runtime(
        ErrorException(
            "In `@constraint(model, x ⟂ 2x - 1)`: second term must be a " *
            "variable.",
        ),
        @constraint(model, x ⟂ 2x - 1)
    )
    return
end

function test_scalar_error_F_F()
    model = Model()
    @variable(model, x >= 0)
    @test_throws_runtime(
        ErrorException(
            "In `@constraint(model, x + 1 ⟂ 2x - 1)`: second term must " *
            "be a variable.",
        ),
        @constraint(model, x + 1 ⟂ 2x - 1)
    )
    return
end

function test_scalar_error_0_F()
    model = Model()
    @variable(model, x >= 0)
    @test_throws_runtime(
        ErrorException(
            "In `@constraint(model, 0 ⟂ 2x - 1)`: second term must be a " *
            "variable.",
        ),
        @constraint(model, 0 ⟂ 2x - 1)
    )
    return
end

function test_vector_complements()
    model = Model()
    @variable(model, x[1:2] >= 0)
    @constraint(model, c, complements(2x .- 1, x))
    obj = constraint_object(c)
    @test obj.func == [2x .- 1; x]
    @test obj.set == MOI.Complements(4)
    return
end

function test_vector_perp()
    model = Model()
    @variable(model, x[1:2] >= 0)
    @variable(model, y[1:2] >= 0)
    @constraint(model, c, x ⟂ y)
    obj = constraint_object(c)
    @test obj.func == [x; y]
    @test obj.set == MOI.Complements(4)
    return
end

function test_vector_error_length_mismatch()
    model = Model()
    @variable(model, x[1:2] >= 0)
    @test_throws_runtime(
        ErrorException(
            "In `@constraint(model, x ⟂ [x[1]])`: size of mapping does " *
            "not match size of variables: (2,) != (1,).",
        ),
        @constraint(model, x ⟂ [x[1]])
    )
    return
end

function test_vector_error_x_F()
    model = Model()
    @variable(model, x[1:2] >= 0)
    @test_throws_runtime(
        ErrorException(
            "In `@constraint(model, x ⟂ 2x .- 1)`: second term must be an " *
            "array of variables.",
        ),
        @constraint(model, x ⟂ 2x .- 1)
    )
end

function test_vector_error_F_F()
    model = Model()
    @variable(model, x[1:2] >= 0)
    @test_throws_runtime(
        ErrorException(
            "In `@constraint(model, x .+ 1 ⟂ 2x .- 1)`: second term must " *
            "be an array of variables.",
        ),
        @constraint(model, x .+ 1 ⟂ 2x .- 1)
    )
    return
end

function test_vector_error_0_F()
    model = Model()
    @variable(model, x[1:2] >= 0)
    y = [1.2, -1.3]
    @test_throws_runtime(
        ErrorException(
            "In `@constraint(model, y ⟂ 2x .- 1)`: second term must " *
            "be an array of variables.",
        ),
        @constraint(model, y ⟂ 2x .- 1)
    )
    return
end

function test_sparse_complements()
    model = Model()
    @variable(model, x[i = 1:3; isodd(i)] >= 0)
    @constraint(model, c, complements(2x .- 1, x))
    obj = constraint_object(c)
    @test obj.func == [2x[3] - 1, 2x[1] - 1, x[3], x[1]] ||
          obj.func == [2x[1] - 1, 2x[3] - 1, x[1], x[3]]
    @test obj.set == MOI.Complements(4)
    return
end

function test_sparse_perp()
    model = Model()
    @variable(model, x[i = 1:3; isodd(i)] >= 0)
    @constraint(model, c, 2x .- 1 ⟂ x)
    obj = constraint_object(c)
    @test obj.func == [2x[3] - 1, 2x[1] - 1, x[3], x[1]] ||
          obj.func == [2x[1] - 1, 2x[3] - 1, x[1], x[3]]
    @test obj.set == MOI.Complements(4)
    return
end

function test_sparse_key_mismatch()
    model = Model()
    @variable(model, x[i = 1:3; isodd(i)] >= 0)
    @variable(model, y[i = 3:5; isodd(i)] >= 0)
    @test_throws(ErrorException, @constraint(model, 2x .- 1 ⟂ y))
    return
end

function test_F_constant_scalar()
    model = Model()
    @variable(model, 0 <= x <= 1)
    @constraint(model, c, 0 ⟂ x)
    obj = constraint_object(c)
    @test obj.func == AffExpr[0, x]
    @test obj.set == MOI.Complements(2)
    return
end

function test_F_constant_vector()
    model = Model()
    @variable(model, 0 <= x[1:2] <= 1)
    F = [1.2, -1.3]
    @constraint(model, c, F ⟂ x)
    obj = constraint_object(c)
    @test obj.func == [F; x]
    @test obj.set == MOI.Complements(4)
    return
end

end  # module

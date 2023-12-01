#  Copyright 2017, Iain Dunning, Joey Huchette, Miles Lubin, and contributors
#  This Source Code Form is subject to the terms of the Mozilla Public
#  License, v. 2.0. If a copy of the MPL was not distributed with this
#  file, You can obtain one at https://mozilla.org/MPL/2.0/.

module TestLPMatrixData

using JuMP
using Test

function test_standard_matrix_form()
    model = Model()
    @variable(model, x >= 1, Bin)
    @variable(model, 2 <= y)
    @variable(model, 3 <= z <= 4, Int)
    @constraint(model, x == 5)
    @constraint(model, 2x + 3y <= 6)
    @constraint(model, -4y >= 5z + 7)
    @constraint(model, -1 <= x + y <= 2)
    @objective(model, Max, 1 + 2x)
    a = lp_matrix_data(model)
    @test Matrix(a.A) == [1 0 0; 0 -4 -5; 2 3 0; 1 1 0]
    @test a.b_lower == [5, 7, -Inf, -1]
    @test a.b_upper == [5, Inf, 6, 2]
    @test a.x_lower == [1, 2, 3]
    @test a.x_upper == [Inf, Inf, 4]
    @test a.c == [2, 0, 0]
    @test a.c_offset == 1
    @test a.sense == MOI.MAX_SENSE
    @test a.integers == [3]
    @test a.binaries == [1]
    @objective(model, Min, y)
    unset_binary(x)
    set_integer(x)
    b = lp_matrix_data(model)
    b.sense == MOI.MIN_SENSE
    b.c == [0, 1, 0]
    b.c_offset == 0
    @test b.integers == [1, 3]
    @test b.binaries == Int[]
    return
end

function test_standard_matrix_form_rational()
    model = GenericModel{Rational{Int}}()
    @variable(model, x >= 1)
    @variable(model, 2 <= y)
    @variable(model, 3 <= z <= 4)
    @constraint(model, x == 5)
    @constraint(model, (2 // 3)x + 3y <= 6 // 5)
    @constraint(model, -4y >= 5z + 7 // 8)
    @constraint(model, -1 <= x + y <= 2)
    @objective(model, Min, (1 // 2) - 2x + (1 // 3) * z)
    a = lp_matrix_data(model)
    @test Matrix(a.A) == [1 0 0; 0 -4 -5; (2//3) 3 0; 1 1 0]
    @test a.b_lower == [5, 7 // 8, -1 // 0, -1]
    @test a.b_upper == [5, 1 // 0, 6 // 5, 2]
    @test a.x_lower == [1, 2, 3]
    @test a.x_upper == [1 // 0, 1 // 0, 4]
    @test a.c == [-2, 0, 1 // 3]
    @test a.c_offset == 1 // 2
    @test a.sense == MOI.MIN_SENSE
    return
end

function test_standard_matrix_form_empty()
    model = Model()
    a = lp_matrix_data(model)
    @test size(a.A) == (0, 0)
    @test isempty(a.b_lower)
    @test isempty(a.b_upper)
    @test isempty(a.x_lower)
    @test isempty(a.x_upper)
    @test isempty(a.c)
    @test a.c_offset == 0
    @test a.sense == MOI.FEASIBILITY_SENSE
    return
end

function test_standard_matrix_form_bad_constraint()
    model = Model()
    @variable(model, x[1:3])
    @constraint(model, x in SecondOrderCone())
    @test_throws(
        ErrorException(
            "Unsupported constraint type in `lp_matrix_data`: $(Vector{VariableRef}) -in- $(MOI.SecondOrderCone)",
        ),
        lp_matrix_data(model),
    )
    return
end

function test_standard_matrix_form_bad_objective()
    model = Model()
    @variable(model, x)
    @objective(model, Min, [2x + 1, 3x])
    @test_throws(
        ErrorException(
            "Unsupported objective type in `lp_matrix_data`: $(Vector{AffExpr})",
        ),
        lp_matrix_data(model),
    )
    return
end

end  # module

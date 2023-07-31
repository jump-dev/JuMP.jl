#  Copyright 2017, Iain Dunning, Joey Huchette, Miles Lubin, and contributors
#  This Source Code Form is subject to the terms of the Mozilla Public
#  License, v. 2.0. If a copy of the MPL was not distributed with this
#  file, You can obtain one at https://mozilla.org/MPL/2.0/.

module TestContainersDimensionalData

using Test

using DimensionalData
using JuMP

function test_dimension_data_variable()
    model = Model()
    @variable(model, x[i = 2:4, j = ["a", "b"]], container = DimArray)
    @test x isa DimArray
    @test x[At(2), At("a")] isa VariableRef
    @test JuMP.name(x[At(4), At("b")]) == "x[4,b]"
    @test @expression(model, sum(x[At(i), At("a")] for i in 2:4)) isa AffExpr
    @constraint(model, c, sum(x[At(i), At("a")] for i in 2:4) <= 1)
    @test c isa ConstraintRef
    return
end

function test_dimension_data_expression()
    model = Model()
    B = ["a", "b"]
    @variable(model, x[i = 2:4, j = B], container = DimArray)
    @expression(
        model,
        expr[j = B],
        sum(x[At(i), At(j)] for i in 2:4),
        container = DimArray,
    )
    @test expr isa DimArray
    @test expr[At("a")] isa AffExpr
    return
end

function test_dimensional_data_missing_names()
    model = Model()
    @test @variable(model, [1:3, 1:2], container = DimArray) isa DimArray
    @test @variable(model, [i = 1:3, 1:2], container = DimArray) isa DimArray
    @test @variable(model, [1:3, j = 1:2], container = DimArray) isa DimArray
    return
end

function test_dimensional_data_sparse()
    model = Model()
    @test_throws(
        ErrorException(
            "Unable to create a `DimArray` because the container does not form " *
            "a dense rectangular array",
        ),
        @variable(model, [i = 1:3, i:2], container = DimArray),
    )
    @test_throws(
        ErrorException(
            "Unable to create a `DimArray` because the container does not form " *
            "a dense rectangular array",
        ),
        @variable(model, [i = 1:3; isodd(i)], container = DimArray),
    )
    return
end

end

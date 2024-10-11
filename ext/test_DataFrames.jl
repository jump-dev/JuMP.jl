#  Copyright 2017, Iain Dunning, Joey Huchette, Miles Lubin, and contributors
#  This Source Code Form is subject to the terms of the Mozilla Public
#  License, v. 2.0. If a copy of the MPL was not distributed with this
#  file, You can obtain one at https://mozilla.org/MPL/2.0/.

module TestContainersDataFrames

using Test

using DataFrames
using JuMP

function test_dimension_data_vector()
    model = Model()
    @variable(model, x[i = 2:4], container = DataFrame)
    @test x isa DataFrame
    @test size(x) == (3, 2)
    @test names(x) == ["i", "value"]
    return
end

function test_dimension_data_matrix()
    model = Model()
    @variable(model, x[i = 2:4, j = ["a", "b"]], container = DataFrame)
    @test x isa DataFrame
    @test size(x) == (6, 3)
    @test names(x) == ["i", "j", "value"]
    @test sum(x[x.j .== "a", :value]) isa AffExpr
    return
end

function test_dimension_data_triangle()
    model = Model()
    @variable(model, x[i = 2:4, j in i:4], container = DataFrame)
    @test x isa DataFrame
    @test size(x) == (6, 3)
    @test names(x) == ["i", "j", "value"]
    return
end

function test_dimension_data_sparse()
    model = Model()
    @variable(model, x[i in 1:4, j in 1:4; isodd(i + j)], container = DataFrame)
    @test x isa DataFrame
    @test size(x) == (8, 3)
    @test x.i == [1, 1, 2, 2, 3, 3, 4, 4]
    @test x.j == [2, 4, 1, 3, 2, 4, 1, 3]
    @test names(x) == ["i", "j", "value"]
    return
end

function test_dataframes_expression()
    model = Model()
    B = ["a", "b"]
    @variable(model, x[i = 2:4, j = B], container = DataFrame)
    @expression(
        model,
        expr[j = B],
        sum(x[x.j .== j, :value]),
        container = DataFrame,
    )
    @test expr isa DataFrame
    @test expr.j == ["a", "b"]
    expr2 = DataFrames.combine(
        DataFrames.groupby(x, :j),
        :value => sum => :value,
    )
    @test expr == expr2
    return
end

function test_data_frames_missing_names()
    model = Model()
    x = @variable(model, [1:3, 1:2], container = DataFrame)
    @test all(startswith.(names(x), ["##", "##", "value"]))
    x = @variable(model, [i in 1:3, 1:2], container = DataFrame)
    @test all(startswith.(names(x), ["i", "##", "value"]))
    x = @variable(model, [1:3, j in 1:2], container = DataFrame)
    @test all(startswith.(names(x), ["##", "j", "value"]))
    return
end

end

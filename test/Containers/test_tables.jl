#  Copyright 2017, Iain Dunning, Joey Huchette, Miles Lubin, and contributors
#  This Source Code Form is subject to the terms of the Mozilla Public
#  License, v. 2.0. If a copy of the MPL was not distributed with this
#  file, You can obtain one at https://mozilla.org/MPL/2.0/.

module TestTableInterface

using JuMP
using Test

function test_denseaxisarray()
    model = Model()
    @variable(model, x[i = 4:10, j = 2002:2022] >= 0, start = 0.0)
    @test typeof(x) <: Containers.DenseAxisArray
    start_table = Containers.rowtable(start_value, x; header = [:i1, :i2, :i3])
    T = NamedTuple{(:i1, :i2, :i3),Tuple{Int,Int,Float64}}
    @test start_table isa Vector{T}
    @test length(start_table) == length(x)
    row = first(start_table)
    @test row == (i1 = 4, i2 = 2002, i3 = 0.0)
    x_table = Containers.rowtable(x; header = [:i1, :i2, :i3])
    @test x_table[1] == (i1 = 4, i2 = 2002, i3 = x[4, 2002])
    return
end

function test_array()
    model = Model()
    @variable(model, x[1:10, 1:5] >= 0, start = 0.0)
    @test typeof(x) <: Array{VariableRef}
    start_table = Containers.rowtable(start_value, x; header = [:i1, :i2, :i3])
    T = NamedTuple{(:i1, :i2, :i3),Tuple{Int,Int,Float64}}
    @test start_table isa Vector{T}
    @test length(start_table) == length(x)
    row = first(start_table)
    @test row == (i1 = 1, i2 = 1, i3 = 0.0)
    x_table = Containers.rowtable(x; header = [:i1, :i2, :i3])
    @test x_table[1] == (i1 = 1, i2 = 1, i3 = x[1, 1])
    return
end

function test_sparseaxisarray()
    model = Model()
    @variable(model, x[i = 1:10, j = 1:5; i + j <= 8] >= 0, start = 0)
    @test typeof(x) <: Containers.SparseAxisArray
    start_table = Containers.rowtable(start_value, x; header = [:i1, :i2, :i3])
    T = NamedTuple{(:i1, :i2, :i3),Tuple{Int,Int,Float64}}
    @test start_table isa Vector{T}
    @test length(start_table) == length(x)
    @test (i1 = 1, i2 = 1, i3 = 0.0) in start_table
    x_table = Containers.rowtable(x; header = [:i1, :i2, :i3])
    @test (i1 = 1, i2 = 1, i3 = x[1, 1]) in x_table
    return
end

function test_col_name_error()
    model = Model()
    @variable(model, x[1:2, 1:2])
    @test_throws ErrorException Containers.rowtable(x; header = [:y, :a])
    @test_throws(
        ErrorException,
        Containers.rowtable(x; header = [:y, :a, :b, :c]),
    )
    @test Containers.rowtable(x; header = [:y, :a, :b]) isa Vector{<:NamedTuple}
    return
end

# Mockup of custom variable type
struct _MockVariable <: AbstractVariable
    var::ScalarVariable
end

struct _MockVariableRef <: AbstractVariableRef
    vref::VariableRef
end

JuMP.name(v::_MockVariableRef) = name(v.vref)

JuMP.owner_model(v::_MockVariableRef) = owner_model(v.vref)

JuMP.start_value(v::_MockVariableRef) = start_value(v.vref)

struct _Mock end

function JuMP.build_variable(::Function, info::VariableInfo, _::_Mock)
    return _MockVariable(ScalarVariable(info))
end

function JuMP.add_variable(model::GenericModel, x::_MockVariable, name::String)
    variable = add_variable(model, x.var, name)
    return _MockVariableRef(variable)
end

function test_custom_variable()
    model = Model()
    @variable(
        model,
        x[i = 1:3, j = 100:102] >= 0,
        _Mock(),
        container = Containers.DenseAxisArray,
        start = 0.0,
    )
    @test typeof(x) <: Containers.DenseAxisArray
    start_table = Containers.rowtable(start_value, x)
    T = NamedTuple{(:x1, :x2, :y),Tuple{Int,Int,Float64}}
    @test start_table isa Vector{T}
    @test length(start_table) == length(x)
    @test (x1 = 1, x2 = 100, y = 0.0) in start_table
    x_table = Containers.rowtable(x)
    @test (x1 = 1, x2 = 100, y = x[1, 100]) in x_table
    return
end

end

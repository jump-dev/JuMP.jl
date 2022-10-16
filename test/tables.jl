#  Copyright 2017, Iain Dunning, Joey Huchette, Miles Lubin, and contributors
#  This Source Code Form is subject to the terms of the Mozilla Public
#  License, v. 2.0. If a copy of the MPL was not distributed with this
#  file, You can obtain one at https://mozilla.org/MPL/2.0/.

module TestTableInterface

using JuMP
using Test

function runtests()
    for name in names(@__MODULE__; all = true)
        if startswith("$(name)", "test_")
            @testset "$(name)" begin
                getfield(@__MODULE__, name)()
            end
        end
    end
    return
end

function test_denseaxisarray()
    model = Model()
    @variable(model, x[i = 4:10, j = 2002:2022] >= 0, start = 0.0)
    @test typeof(x) <: Containers.DenseAxisArray
    start_table = JuMP.table(start_value, x, :solution, :index1, :index2)
    row = first(start_table)
    @test eltype(start_table) == typeof(row)
    @test row.solution == 0
    @test row.index1 == 4
    @test propertynames(row) == (:index1, :index2, :solution)
    rows = collect(start_table)
    @test length(rows) == length(start_table)
    var_tbl = JuMP.table(x, :variable, :index1, :index2)
    @test typeof(first(var_tbl).variable) <: VariableRef
    return
end

function test_array()
    model = Model()
    @variable(model, x[1:10, 1:5] >= 0, start = 0.0)
    @test typeof(x) <: Array{VariableRef}
    start_table = JuMP.table(start_value, x, :solution, :index1, :index2)
    row = first(start_table)
    @test eltype(start_table) == typeof(row)
    @test row.solution == 0
    @test row.index1 == 1
    @test propertynames(row) == (:index1, :index2, :solution)
    rows = collect(start_table)
    @test length(rows) == length(start_table)
    return
end

function test_sparseaxisarray()
    model = Model()
    @variable(model, x[i = 1:10, j = 1:5; i + j <= 8] >= 0, start = 0)
    @test typeof(x) <: Containers.SparseAxisArray
    start_table = JuMP.table(start_value, x, :solution, :index1, :index2)
    row = first(start_table)
    @test eltype(start_table) == typeof(row)
    @test row.solution == 0.0
    @test row.index1 == 1
    @test propertynames(row) == (:index1, :index2, :solution)
    rows = collect(start_table)
    @test length(rows) == length(start_table)
    return
end

function test_col_name_error()
    model = Model()
    @variable(model, x[1:2, 1:2])
    @test_throws ErrorException table(x, :y, :a)
    @test_throws ErrorException table(x, :y, :a, :b, :c)
    @test table(x, :y, :a, :b) isa Vector{<:NamedTuple}
    return
end

# Mockup of custom variable type
struct _MockVariable <: JuMP.AbstractVariable
    var::JuMP.ScalarVariable
end

struct _MockVariableRef <: JuMP.AbstractVariableRef
    vref::VariableRef
end

JuMP.name(v::_MockVariableRef) = JuMP.name(v.vref)

JuMP.owner_model(v::_MockVariableRef) = JuMP.owner_model(v.vref)

JuMP.start_value(v::_MockVariableRef) = JuMP.start_value(v.vref)

struct _Mock end

function JuMP.build_variable(::Function, info::JuMP.VariableInfo, _::_Mock)
    return _MockVariable(JuMP.ScalarVariable(info))
end

function JuMP.add_variable(model::Model, x::_MockVariable, name::String)
    variable = JuMP.add_variable(model, x.var, name)
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
    start_table = JuMP.table(start_value, x, :solution, :index1, :index2)
    row = first(start_table)
    @test eltype(start_table) == typeof(row)
    @test row.solution == 0.0
    @test row.index1 == 1
    @test propertynames(row) == (:index1, :index2, :solution)
    return
end

end

TestTableInterface.runtests()

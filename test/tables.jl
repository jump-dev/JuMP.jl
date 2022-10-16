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
    tbl = JuMP.table(start_value, x, :solution, :index1, :index2)
    tblrow = first(tbl)
    @test eltype(tbl) == typeof(tblrow)
    @test tblrow.solution == 0
    @test tblrow.index1 == 4
    @test propertynames(tblrow) == (:index1, :index2, :solution)
    rows = collect(tbl)
    @test length(rows) == length(tbl)
    var_tbl = JuMP.table(x, :variable, :index1, :index2)
    @test typeof(first(var_tbl).variable) <: VariableRef
    return
end

function test_array()
    model = Model()
    @variable(model, x[1:10, 1:5] >= 0, start = 0.0)
    @test typeof(x) <: Array{VariableRef}
    tbl = JuMP.table(start_value, x, :solution, :index1, :index2)
    tblrow = first(tbl)
    @test eltype(tbl) == typeof(tblrow)
    @test tblrow.solution == 0
    @test tblrow.index1 == 1
    @test propertynames(tblrow) == (:index1, :index2, :solution)
    rows = collect(tbl)
    @test length(rows) == length(tbl)
    return
end

function test_sparseaxisarray()
    model = Model()
    @variable(model, x[i = 1:10, j = 1:5; i + j <= 8] >= 0, start = 0)
    @test typeof(x) <: Containers.SparseAxisArray
    tbl = JuMP.table(start_value, x, :solution, :index1, :index2)
    tblrow = first(tbl)
    @test eltype(tbl) == typeof(tblrow)
    @test tblrow.solution == 0.0
    @test tblrow.index1 == 1
    @test propertynames(tblrow) == (:index1, :index2, :solution)
    rows = collect(tbl)
    @test length(rows) == length(tbl)
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

JuMP.value(v::_MockVariableRef) = JuMP.value(v.vref)

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
    tbl = JuMP.table(start_value, x, :solution, :index1, :index2)
    tblrow = first(tbl)
    @test eltype(tbl) == typeof(tblrow)
    @test tblrow.solution == 0.0
    @test tblrow.index1 == 1
    @test propertynames(tblrow) == (:index1, :index2, :solution)
    return
end

end

TestTableInterface.runtests()

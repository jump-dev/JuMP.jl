module TestTableInterface

using JuMP
using Tables
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
    @variable(model, x[i = 4:10, j = 2002:2022] >= 0)
    @objective(model, Min, sum(x))
    set_optimizer(
        model,
        () -> MOI.Utilities.MockOptimizer(
            MOI.Utilities.Model{Float64}();
            eval_objective_value = false,
        ),
    )
    optimize!(model)
    mockoptimizer = JuMP.unsafe_backend(model)
    MOI.set(mockoptimizer, MOI.TerminationStatus(), MOI.OPTIMAL)
    MOI.set(mockoptimizer, MOI.ResultCount(), 1)
    MOI.set(mockoptimizer, MOI.PrimalStatus(), MOI.FEASIBLE_POINT)
    MOI.set(mockoptimizer, MOI.DualStatus(), MOI.FEASIBLE_POINT)

    for ind in eachindex(x)
        MOI.set(
            mockoptimizer,
            MOI.VariablePrimal(),
            JuMP.optimizer_index(x[ind]),
            0.0,
        )
    end

    tbl = JuMP.solution_table(x, :solution, :index1, :index2)
    @test Tables.istable(typeof(tbl))
    @test Tables.rowaccess(typeof(tbl))

    tblrow = first(tbl)
    @test eltype(tbl) == typeof(tblrow)
    @test Tables.getcolumn(tblrow, :index1) == 4
    @test Tables.getcolumn(tblrow, 1) == 4
    @test tblrow.index1 == 4
    @test propertynames(tblrow) == [:index1, :index2, :solution]

    rows = collect(tbl)
    @test length(rows) == length(tbl)
end

function test_array()
    model = Model()
    @variable(model, x[1:10, 1:5] >= 0)
    @test typeof(x) <: Array{VariableRef}

    @objective(model, Min, sum(x))
    set_optimizer(
        model,
        () -> MOI.Utilities.MockOptimizer(
            MOI.Utilities.Model{Float64}();
            eval_objective_value = false,
        ),
    )
    optimize!(model)
    mockoptimizer = JuMP.unsafe_backend(model)
    MOI.set(mockoptimizer, MOI.TerminationStatus(), MOI.OPTIMAL)
    MOI.set(mockoptimizer, MOI.ResultCount(), 1)
    MOI.set(mockoptimizer, MOI.PrimalStatus(), MOI.FEASIBLE_POINT)
    MOI.set(mockoptimizer, MOI.DualStatus(), MOI.FEASIBLE_POINT)

    for ind in eachindex(x)
        MOI.set(
            mockoptimizer,
            MOI.VariablePrimal(),
            JuMP.optimizer_index(x[ind]),
            0.0,
        )
    end

    tbl = JuMP.solution_table(x, :solution, :index1, :index2)

    @test Tables.istable(typeof(tbl))
    @test Tables.rowaccess(typeof(tbl))

    tblrow = first(tbl)
    @test eltype(tbl) == typeof(tblrow)
    @test Tables.getcolumn(tblrow, :index1) == 1
    @test Tables.getcolumn(tblrow, 1) == 1
    @test tblrow.index1 == 1
    @test propertynames(tblrow) == [:index1, :index2, :solution]

    rows = collect(tbl)
    @test length(rows) == length(tbl)
end

function test_sparseaxisarray()
    model = Model()
    @variable(model, x[i = 1:10, j = 1:5; i + j <= 8] >= 0)
    @test typeof(x) <: Containers.SparseAxisArray
    @objective(model, Min, sum(x))
    set_optimizer(
        model,
        () -> MOI.Utilities.MockOptimizer(
            MOI.Utilities.Model{Float64}();
            eval_objective_value = false,
        ),
    )
    optimize!(model)
    mockoptimizer = JuMP.unsafe_backend(model)
    MOI.set(mockoptimizer, MOI.TerminationStatus(), MOI.OPTIMAL)
    MOI.set(mockoptimizer, MOI.ResultCount(), 1)
    MOI.set(mockoptimizer, MOI.PrimalStatus(), MOI.FEASIBLE_POINT)
    MOI.set(mockoptimizer, MOI.DualStatus(), MOI.FEASIBLE_POINT)

    for ind in eachindex(x)
        MOI.set(
            mockoptimizer,
            MOI.VariablePrimal(),
            JuMP.optimizer_index(x[ind]),
            0.0,
        )
    end

    tbl = JuMP.solution_table(x, :solution, :index1, :index2)
    @test Tables.istable(typeof(tbl))
    @test Tables.rowaccess(typeof(tbl))

    tblrow = first(tbl)
    @test eltype(tbl) == typeof(tblrow)
    @test Tables.getcolumn(tblrow, :index1) == 1
    @test Tables.getcolumn(tblrow, 1) == 1
    @test tblrow.index1 == 1
    @test propertynames(tblrow) == [:index1, :index2, :solution]

    rows = collect(tbl)
    @test length(rows) == length(tbl)
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
        container = Containers.DenseAxisArray
    )

    @objective(model, Min, 0)
    set_optimizer(
        model,
        () -> MOI.Utilities.MockOptimizer(
            MOI.Utilities.Model{Float64}();
            eval_objective_value = false,
        ),
    )
    optimize!(model)
    mockoptimizer = JuMP.unsafe_backend(model)
    MOI.set(mockoptimizer, MOI.TerminationStatus(), MOI.OPTIMAL)
    MOI.set(mockoptimizer, MOI.ResultCount(), 1)
    MOI.set(mockoptimizer, MOI.PrimalStatus(), MOI.FEASIBLE_POINT)
    MOI.set(mockoptimizer, MOI.DualStatus(), MOI.FEASIBLE_POINT)

    for ind in eachindex(x)
        MOI.set(
            mockoptimizer,
            MOI.VariablePrimal(),
            JuMP.optimizer_index(x[ind].vref),
            0.0,
        )
    end

    tbl = JuMP.solution_table(x, :solution, :index1, :index2)
    @test Tables.istable(typeof(tbl))
    @test Tables.rowaccess(typeof(tbl))

    tblrow = first(tbl)
    @test eltype(tbl) == typeof(tblrow)
    @test Tables.getcolumn(tblrow, :index1) == 1
    @test Tables.getcolumn(tblrow, 1) == 1
    @test tblrow.index1 == 1
    @test propertynames(tblrow) == [:index1, :index2, :solution]
end

end

TestTableInterface.runtests()

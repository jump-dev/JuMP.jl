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
    @variable(model, x[i = 4:10, j=2002:2022] >= 0)
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
        MOI.set(mockoptimizer, MOI.VariablePrimal(), JuMP.optimizer_index(x[ind]), 0.0)
    end
    
    t = JuMP.table(x, :solution, :index1, :index2)
    
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
        MOI.set(mockoptimizer, MOI.VariablePrimal(), JuMP.optimizer_index(x[ind]), 0.0)
    end
    
    t = JuMP.table(x, :solution, :index1, :index2)
end

function test_sparseaxisarray()
    model = Model()
    @variable(model, x[i=1:10, j=1:5; i + j <= 8] >= 0)
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
        MOI.set(mockoptimizer, MOI.VariablePrimal(), JuMP.optimizer_index(x[ind]), 0.0)
    end
    
    t = JuMP.table(x, :solution, :index1, :index2)
end

end
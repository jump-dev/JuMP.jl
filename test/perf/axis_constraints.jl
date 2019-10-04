using JuMP
using BenchmarkTools

function sum_iterate(con_refs)
    x = 0.0
    for con_ref in con_refs
        x += dual(con_ref)
    end
    return x
end

function sum_index(con_refs)
    x = 0.0
    for i in eachindex(con_refs)
        x += dual(con_refs[i])
    end
    return x
end

function dense_axis_constraints(n)
    model = Model()
    mock = MOIU.MockOptimizer(MOIU.Model{Float64}(),
                              eval_variable_constraint_dual=false)
    MOIU.set_mock_optimize!(mock,
        mock -> MOIU.mock_optimize!(mock, zeros(n),
            (MOI.SingleVariable, MOI.EqualTo{Float64}) => ones(n - 1)))
    MOIU.reset_optimizer(model, mock)

    @variable(model, x[1:n])
    set = MOI.EqualTo(0.0)
    con_refs = @time @constraint(model, [i = 2:n], x[i] in set)
    optimize!(model)
    @assert sum_iterate(con_refs) == n - 1
    @btime sum_iterate($con_refs)
    @assert sum_index(con_refs) == n - 1
    @btime sum_index($con_refs)
    return
end
function sparse_axis_constraints(n)
    model = Model()
    mock = MOIU.MockOptimizer(MOIU.Model{Float64}(),
                              eval_variable_constraint_dual=false)
    MOIU.set_mock_optimize!(mock,
        mock -> MOIU.mock_optimize!(mock, zeros(n),
            (MOI.SingleVariable, MOI.EqualTo{Float64}) => ones(div(n, 2))))
    MOIU.reset_optimizer(model, mock)

    @variable(model, x[1:n])
    set = MOI.EqualTo(0.0)
    con_refs = @time @constraint(model, [i = 1:n; iseven(i)], x[i] in set)
    optimize!(model)
    @assert sum_iterate(con_refs) == div(n, 2)
    @btime sum_iterate($con_refs)
    @assert sum_index(con_refs) == div(n, 2)
    @btime sum_index($con_refs)
    return
end

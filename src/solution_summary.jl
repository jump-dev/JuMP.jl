#  Copyright 2017, Iain Dunning, Joey Huchette, Miles Lubin, and contributors
#  This Source Code Form is subject to the terms of the Mozilla Public
#  License, v. 2.0. If a copy of the MPL was not distributed with this
#  file, You can obtain one at https://mozilla.org/MPL/2.0/.

struct _SolutionSummary{T}
    result::Int
    verbose::Bool
    solver::String
    # Status
    termination_status::MOI.TerminationStatusCode
    primal_status::MOI.ResultStatusCode
    dual_status::MOI.ResultStatusCode
    raw_status::String
    result_count::Int
    has_values::Bool
    has_duals::Bool
    # Candidate solution
    objective_value::Union{Missing,T,Vector{T}}
    objective_bound::Union{Missing,T,Vector{T}}
    relative_gap::Union{Missing,T}
    dual_objective_value::Union{Missing,T}
    primal_solution::Union{Missing,Dict{String,T}}
    dual_solution::Union{Missing,Dict{String,T}}
    # Work counters
    solve_time::Union{Missing,Float64}
    barrier_iterations::Union{Missing,Int}
    simplex_iterations::Union{Missing,Int}
    node_count::Union{Missing,Int}
    # The default construction `_SolutionSummary(...)` wouldn't work
    # as `T` would be unbound if `missing` is passed for all possible fields
    # hence `Aqua` complains.
    function _SolutionSummary{T}(args...) where {T}
        return new{T}(args...)
    end
end

"""
    solution_summary(model::GenericModel; result::Int = 1, verbose::Bool = false)

Return a struct that can be used print a summary of the solution in result
`result`.

If `verbose=true`, write out the primal solution for every variable and the
dual solution for every constraint, excluding those with empty names.

## Example

When called at the REPL, the summary is automatically printed:
```jldoctest
julia> model = Model();

julia> solution_summary(model)
* Solver : No optimizer attached.

* Status
  Result count       : 0
  Termination status : OPTIMIZE_NOT_CALLED
  Message from the solver:
  "optimize not called"

* Candidate solution (result #1)
  Primal status      : NO_SOLUTION
  Dual status        : NO_SOLUTION

* Work counters
```

Use `print` to force the printing of the summary from inside a function:
```jldoctest
julia> model = Model();

julia> function foo(model)
           print(solution_summary(model))
           return
       end
foo (generic function with 1 method)

julia> foo(model)
* Solver : No optimizer attached.

* Status
  Result count       : 0
  Termination status : OPTIMIZE_NOT_CALLED
  Message from the solver:
  "optimize not called"

* Candidate solution (result #1)
  Primal status      : NO_SOLUTION
  Dual status        : NO_SOLUTION

* Work counters
```
"""
function solution_summary(
    model::GenericModel{T};
    result::Int = 1,
    verbose::Bool = false,
) where {T}
    num_results = result_count(model)
    has_primal = has_values(model; result = result)
    has_dual = has_duals(model; result = result)
    return _SolutionSummary{T}(
        result,
        verbose,
        solver_name(model),
        termination_status(model),
        primal_status(model; result = result),
        dual_status(model; result = result),
        raw_status(model),
        num_results,
        has_primal,
        has_dual,
        _try_get(m -> objective_value(m; result = result), model),
        _try_get(objective_bound, model),
        _try_get(relative_gap, model),
        _try_get(m -> dual_objective_value(m; result = result), model),
        verbose && has_primal ? _get_solution_dict(model, result) : missing,
        verbose && has_dual ? _get_constraint_dict(model, result) : missing,
        _try_get(solve_time, model),
        _try_get(barrier_iterations, model),
        _try_get(simplex_iterations, model),
        _try_get(node_count, model),
    )
end

"""
    Base.show([io::IO], summary::SolutionSummary; verbose::Bool = false)

Write a summary of the solution results to `io` (or to `stdout` if `io` is not
given).
"""
function Base.show(io::IO, summary::_SolutionSummary)
    println(io, "* Solver : ", summary.solver)
    println(io)
    _show_status_summary(io, summary)
    _show_candidate_solution_summary(io, summary)
    if summary.result == 1
        _show_work_counters_summary(io, summary)
    end
    return
end

function _show_status_summary(io::IO, summary::_SolutionSummary)
    println(io, "* Status")
    println(io, "  Result count       : ", summary.result_count)
    println(io, "  Termination status : ", summary.termination_status)
    if summary.result == 1
        println(io, "  Message from the solver:")
        println(io, "  \"", summary.raw_status, "\"")
    end
    println(io)
    return
end

function _show_candidate_solution_summary(io::IO, summary::_SolutionSummary)
    println(io, "* Candidate solution (result #$(summary.result))")
    println(io, "  Primal status      : ", summary.primal_status)
    println(io, "  Dual status        : ", summary.dual_status)
    _print_if_not_missing(
        io,
        "  Objective value    : ",
        summary.objective_value,
    )
    if summary.result == 1
        _print_if_not_missing(
            io,
            "  Objective bound    : ",
            summary.objective_bound,
        )
        _print_if_not_missing(
            io,
            "  Relative gap       : ",
            summary.relative_gap,
        )
    end
    _print_if_not_missing(
        io,
        "  Dual objective value : ",
        summary.dual_objective_value,
    )
    if summary.verbose && summary.has_values
        println(io, "  Primal solution :")
        for variable_name in sort(collect(keys(summary.primal_solution)))
            _print_if_not_missing(
                io,
                "    $(variable_name) : ",
                summary.primal_solution[variable_name],
            )
        end
    end
    if summary.verbose && summary.has_duals
        println(io, "  Dual solution :")
        for constraint_name in sort(collect(keys(summary.dual_solution)))
            _print_if_not_missing(
                io,
                "    $(constraint_name) : ",
                summary.dual_solution[constraint_name],
            )
        end
    end
    return
end

function _show_work_counters_summary(io::IO, summary::_SolutionSummary)
    println(io)
    println(io, "* Work counters")
    _print_if_not_missing(io, "  Solve time (sec)   : ", summary.solve_time)
    _print_if_not_missing(
        io,
        "  Simplex iterations : ",
        summary.simplex_iterations,
    )
    _print_if_not_missing(
        io,
        "  Barrier iterations : ",
        summary.barrier_iterations,
    )
    _print_if_not_missing(io, "  Node count         : ", summary.node_count)
    return
end

function _get_solution_dict(model, result)
    dict = Dict{String,value_type(typeof(model))}()
    for x in all_variables(model)
        variable_name = name(x)
        if !isempty(variable_name)
            dict[variable_name] = value(x; result = result)
        end
    end
    return dict
end

function _get_constraint_dict(model, result)
    dict = Dict{String,value_type(typeof(model))}()
    for (F, S) in list_of_constraint_types(model)
        for constraint in all_constraints(model, F, S)
            constraint_name = name(constraint)
            if !isempty(constraint_name)
                dict[constraint_name] = dual(constraint; result = result)
            end
        end
    end
    return dict
end

function _try_get(f, model)
    try
        return f(model)
    catch
        return missing
    end
end

_print_if_not_missing(io, header, ::Missing) = nothing
_print_if_not_missing(io, header, value::Real) = println(io, header, value)
function _print_if_not_missing(io, header, value::AbstractFloat)
    println(io, header, Printf.@sprintf("%.5e", value))
    return
end

function _print_if_not_missing(io, header, value::Vector{<:Real})
    array = join([Printf.@sprintf("%.5e", v) for v in value], ",")
    println(io, header, "[", array, "]")
    return
end

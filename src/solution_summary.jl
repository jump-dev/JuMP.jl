#  Copyright 2017, Iain Dunning, Joey Huchette, Miles Lubin, and contributors
#  This Source Code Form is subject to the terms of the Mozilla Public
#  License, v. 2.0. If a copy of the MPL was not distributed with this
#  file, You can obtain one at https://mozilla.org/MPL/2.0/.

import Printf

struct _SolutionSummary
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
    objective_value::Union{Missing,Float64}
    objective_bound::Union{Missing,Float64}
    dual_objective_value::Union{Missing,Float64}
    primal_solution::Union{Missing,Dict{String,Float64}}
    dual_solution::Union{Missing,Dict{String,Float64}}
    # Work counters
    solve_time::Union{Missing,Float64}
    barrier_iterations::Union{Missing,Int}
    simplex_iterations::Union{Missing,Int}
    node_count::Union{Missing,Int}
end

"""
    solution_summary(model::Model; verbose::Bool = false)

Return a struct that can be used print a summary of the solution.

If `verbose=true`, write out the primal solution for every variable and the
dual solution for every constraint, excluding those with empty names.

## Examples

When called at the REPL, the summary is automatically printed:
```julia
julia> solution_summary(model)
[...]
```

Use `print` to force the printing of the summary from inside a function:
```julia
function foo(model)
    print(solution_summary(model))
    return
end
```
"""
function solution_summary(model::Model; verbose::Bool = false)
    return _SolutionSummary(
        verbose,
        solver_name(model),
        termination_status(model),
        primal_status(model),
        dual_status(model),
        raw_status(model),
        result_count(model),
        has_values(model),
        has_duals(model),
        _try_get(objective_value, model),
        _try_get(objective_bound, model),
        _try_get(dual_objective_value, model),
        verbose ? _get_solution_dict(model) : missing,
        verbose ? _get_constraint_dict(model) : missing,
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
    _show_work_counters_summary(io, summary)
    return
end

function _show_status_summary(io::IO, summary::_SolutionSummary)
    println(io, "* Status")
    println(io, "  Termination status : ", summary.termination_status)
    println(io, "  Primal status      : ", summary.primal_status)
    println(io, "  Dual status        : ", summary.dual_status)
    if summary.verbose
        println(io, "  Result count       : ", summary.result_count)
        println(io, "  Has duals          : ", summary.has_duals)
    end
    println(io, "  Message from the solver:")
    println(io, "  \"", summary.raw_status, "\"")
    println(io)
    return
end

function _show_candidate_solution_summary(io::IO, summary::_SolutionSummary)
    println(io, "* Candidate solution")
    _print_if_not_missing(
        io,
        "  Objective value      : ",
        summary.objective_value,
    )
    _print_if_not_missing(
        io,
        "  Objective bound      : ",
        summary.objective_bound,
    )
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
    println(io)
    return
end

function _show_work_counters_summary(io::IO, summary::_SolutionSummary)
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

function _get_solution_dict(model)
    dict = Dict{String,Float64}()
    if has_values(model)
        for x in all_variables(model)
            variable_name = name(x)
            if !isempty(variable_name)
                dict[variable_name] = value(x)
            end
        end
    end
    return dict
end

function _get_constraint_dict(model)
    dict = Dict{String,Float64}()
    if has_duals(model)
        for (F, S) in list_of_constraint_types(model)
            for constraint in all_constraints(model, F, S)
                constraint_name = name(constraint)
                if !isempty(constraint_name)
                    dict[constraint_name] = dual(constraint)
                end
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
_print_if_not_missing(io, header, value::Int) = println(io, header, value)
function _print_if_not_missing(io, header, value::Real)
    println(io, header, Printf.@sprintf("%.5e", value))
    return
end

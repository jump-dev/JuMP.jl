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
    primal_solution::Union{Missing,Dict{String,Union{Missing,T}}}
    dual_solution::Union{Missing,Dict{String,Any}}
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
solution_summary(; result = 1, verbose = false)
├ solver_name          : No optimizer attached.
├ Termination
│ ├ termination_status : OPTIMIZE_NOT_CALLED
│ ├ result_count       : 0
│ └ raw_status         : optimize not called
└ Solution (result = 1)
  ├ primal_status        : NO_SOLUTION
  └ dual_status          : NO_SOLUTION
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
solution_summary(; result = 1, verbose = false)
├ solver_name          : No optimizer attached.
├ Termination
│ ├ termination_status : OPTIMIZE_NOT_CALLED
│ ├ result_count       : 0
│ └ raw_status         : optimize not called
└ Solution (result = 1)
  ├ primal_status        : NO_SOLUTION
  └ dual_status          : NO_SOLUTION
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
    branches = Pair{String,Any}[
        "solver_name          : "=>summary.solver,
        "Termination"=>Pair{String,Any}[
            "termination_status : "=>summary.termination_status,
            "result_count       : "=>summary.result_count,
            "raw_status         : "=>summary.raw_status,
            "objective_bound    : "=>summary.objective_bound,
        ],
        "Solution (result = $(summary.result))"=>Pair{String,Any}[
            "primal_status        : "=>summary.primal_status,
            "dual_status          : "=>summary.dual_status,
            "objective_value      : "=>summary.objective_value,
            "dual_objective_value : "=>summary.dual_objective_value,
        ],
        "Work counters"=>Pair{String,Any}[
            "solve_time (sec)   : "=>summary.solve_time,
            "simplex_iterations : "=>summary.simplex_iterations,
            "barrier_iterations : "=>summary.barrier_iterations,
            "node_count         : "=>summary.node_count,
        ],
    ]
    if summary.result == 1
        push!(
            last(branches[3]),
            "relative_gap         : " => summary.relative_gap,
        )
    end
    if summary.verbose && summary.has_values
        primal_solution = Pair{String,Any}[
            "$name : " => coalesce(
                summary.primal_solution[name],
                "multiple variables with the same name",
            ) for name in sort(collect(keys(summary.primal_solution)))
        ]
        push!(last(branches[3]), "value" => primal_solution)
    end
    if summary.verbose && summary.has_duals
        dual_solution = Pair{String,Any}[
            "$name : " => coalesce(
                summary.dual_solution[name],
                "multiple constraints with the same name",
            ) for name in sort(collect(keys(summary.dual_solution)))
        ]
        push!(last(branches[3]), "dual" => dual_solution)
    end
    if summary.result != 1
        branches = branches[3:3]
    end
    _print_tree(
        io,
        "solution_summary(; result = $(summary.result), verbose = $(summary.verbose))" =>
            branches,
    )
    return
end

function _replace_prefix(prefix)
    if isempty(prefix)
        return prefix
    elseif endswith(prefix, "├ ")
        return string(chop(prefix; tail = 2), "│ ")
    else
        return string(chop(prefix; tail = 2), "  ")
    end
end

_should_keep(::Missing) = false

_should_keep(::Any) = true

_should_keep(x::Pair{String,<:Any}) = _should_keep(last(x))

_should_keep(x::Vector) = any(_should_keep, x)

_format(x::Any) = x

_format(x::AbstractFloat) = Printf.@sprintf("%.5e", x)

function _format(x::Vector{<:Real})
    return string("[", join((Printf.@sprintf("%.5e", v) for v in x), ","), "]")
end

function _print_tree(io::IO, args::Pair{String}, prefix = "")
    if !isempty(prefix)
        println(io)
    end
    if !(last(args) isa Vector{Pair{String,Any}})  # Leaf node
        print(io, prefix, first(args), _format(last(args)))
        return
    end
    branches = filter(_should_keep, last(args))
    print(io, prefix, first(args))
    for (i, branch) in enumerate(branches)
        suffix = i == length(branches) ? "└ " : "├ "
        _print_tree(io, branch, _replace_prefix(prefix) * suffix)
    end
    return
end

function _get_solution_dict(model, result)
    dict = Dict{String,Union{Missing,value_type(typeof(model))}}()
    for x in all_variables(model)
        variable_name = name(x)
        if isempty(variable_name)
            continue
        elseif haskey(dict, variable_name)
            dict[variable_name] = missing
        else
            dict[variable_name] = value(x; result = result)
        end
    end
    return dict
end

function _get_constraint_dict(model, result)
    dict = Dict{String,Any}()
    for (F, S) in list_of_constraint_types(model)
        for constraint in all_constraints(model, F, S)
            constraint_name = name(constraint)
            if isempty(constraint_name)
                continue
            elseif haskey(dict, constraint_name)
                dict[constraint_name] = missing
            else
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

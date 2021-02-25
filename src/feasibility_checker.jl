#  Copyright 2017, Iain Dunning, Joey Huchette, Miles Lubin, and contributors
#  This Source Code Form is subject to the terms of the Mozilla Public
#  License, v. 2.0. If a copy of the MPL was not distributed with this
#  file, You can obtain one at https://mozilla.org/MPL/2.0/.

using LinearAlgebra

"""
    primal_feasibility_report(
        model::Model,
        [point::AbstractDict{VariableRef,Float64}];
        atol::Float64 = 0.0,
        skip_missing::Bool = false,
    )::Dict{Any,Float64}

Given a dictionary `point`, which maps variables to primal values, return a
dictionary mapping the constraint reference of each constraint in `model` to the
distance between the point and the nearest feasible point, if the distance is
greater than `atol`.

## Notes

 * If `skip_missing = true`, constraints containing variables that are not in
   `point` will be ignored.
 * If no point is provided, the primal solution from the last time the model was
   solved is used.

## Examples

```jldoctest; setup=:(using JuMP)
julia> model = Model();

julia> @variable(model, 0.5 <= x <= 1);

julia> primal_feasibility_report(model, Dict(x => 0.2))
Dict{Any,Float64} with 1 entry:
  x ≥ 0.5 => 0.3
```
"""
function primal_feasibility_report(
    model::Model,
    point::AbstractDict{VariableRef,Float64};
    atol::Float64 = 0.0,
    skip_missing::Bool = false,
)
    function point_f(x)
        fx = get(point, x, missing)
        if ismissing(fx) && !skip_missing
            error(
                "point does not contain a value for variable $x. Provide a " *
                "value, or pass `skip_missing = true`.",
            )
        end
        return fx
    end
    violated_constraints = Dict{Any,Float64}()
    for (F, S) in list_of_constraint_types(model)
        _add_infeasible_constraints(
            model, F, S, violated_constraints, point_f, atol
        )
    end
    if num_nl_constraints(model) > 0
        if skip_missing
            error(
                "`skip_missing = true` is not allowed when nonlinear " *
                "constraints are present.",
            )
        end
        _add_infeasible_nonlinear_constraints(
            model, violated_constraints, point_f, atol
        )
    end
    return violated_constraints
end

function primal_feasibility_report(model::Model; kwargs...)
    if !has_values(model)
        error(
            "No primal solution is available. You must provide a point at " *
            "which to check feasibility."
        )
    end
    point = Dict(v => value(v) for v in all_variables(model))
    return primal_feasibility_report(model, point; kwargs...)
end

function _add_infeasible_constraints(
    model::Model,
    ::Type{F},
    ::Type{S},
    violated_constraints::Dict{Any,Float64},
    point_f::Function,
    atol::Float64,
) where {F,S}
    for con in all_constraints(model, F, S)
        obj = constraint_object(con)
        d = _distance_to_set(value.(obj.func, point_f), obj.set)
        if d > atol
            violated_constraints[con] = d
        end
    end
    return
end

function _add_infeasible_nonlinear_constraints(
    model::Model,
    violated_constraints::Dict{Any,Float64},
    point_f::Function,
    atol::Float64,
)
    evaluator = NLPEvaluator(model)
    MOI.initialize(evaluator, Symbol[])
    g = zeros(num_nl_constraints(model))
    MOI.eval_constraint(evaluator, g, point_f.(all_variables(model)))
    for (i, con) in enumerate(model.nlp_data.nlconstr)
        d = max(0.0, con.lb - g[i], g[i] - con.ub)
        if d > atol
            cref = ConstraintRef(
                model,
                NonlinearConstraintIndex(i),
                ScalarShape(),
            )
            violated_constraints[cref] = d
        end
    end
    return
end

function _distance_to_set(::Any, set::MOI.AbstractSet)
    error(
        "Feasibility checker for set type $(typeof(set)) has not been " *
        "implemented yet."
    )
end

_distance_to_set(::Missing, ::MOI.AbstractSet) = 0.0

###
### MOI.AbstractScalarSets
###

function _distance_to_set(x::T, set::MOI.LessThan{T}) where {T<:Real}
    return max(x - set.upper, zero(T))
end

function _distance_to_set(x::T, set::MOI.GreaterThan{T}) where {T<:Real}
    return max(set.lower - x, zero(T))
end

function _distance_to_set(x::T, set::MOI.EqualTo{T}) where {T<:Number}
    return abs(set.value - x)
end

function _distance_to_set(x::T, set::MOI.Interval{T}) where {T<:Real}
    return max(x - set.upper, set.lower - x, zero(T))
end

function _distance_to_set(x::T, ::MOI.ZeroOne) where {T<:Real}
    return min(abs(x - zero(T)), abs(x - one(T)))
end

function _distance_to_set(x::T, ::MOI.Integer) where {T<:Real}
    return abs(x - round(x))
end

function _distance_to_set(x::T, set::MOI.Semicontinuous{T}) where {T<:Real}
    return min(max(x - set.upper, set.lower - x, zero(T)), abs(x))
end

function _distance_to_set(x::T, set::MOI.Semiinteger{T}) where {T<:Real}
    d = max(
        ceil(set.lower) - x,
        x - floor(set.upper),
        abs(x - round(x)),
    )
    return min(d, abs(x))
end

###
### MOI.AbstractVectorSets
###

function _check_dimension(v::AbstractVector, s)
    if length(v) != MOI.dimension(s)
        throw(DimensionMismatch("Mismatch between value and set"))
    end
    return
end

function _distance_to_set(
    x::Vector{T},
    set::MOI.Nonnegatives,
) where {T<:Real}
    _check_dimension(x, set)
    return LinearAlgebra.norm(max(-xi, zero(T)) for xi in x)
end

function _distance_to_set(
    x::Vector{T},
    set::MOI.Nonpositives,
) where {T<:Real}
    _check_dimension(x, set)
    return LinearAlgebra.norm(max(xi, zero(T)) for xi in x)
end

function _distance_to_set(
    x::Vector{T},
    set::MOI.Zeros
) where {T<:Number}
    _check_dimension(x, set)
    return LinearAlgebra.norm(x)
end


function _distance_to_set(
    x::Vector{T},
    set::MOI.Reals
) where {T<:Real}
    _check_dimension(x, set)
    return zero(T)
end


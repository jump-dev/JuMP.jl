#  Copyright 2017, Iain Dunning, Joey Huchette, Miles Lubin, and contributors
#  This Source Code Form is subject to the terms of the Mozilla Public
#  License, v. 2.0. If a copy of the MPL was not distributed with this
#  file, You can obtain one at https://mozilla.org/MPL/2.0/.

function _last_primal_solution(model::GenericModel)
    if !has_values(model)
        error(
            "No primal solution is available. You must provide a point at " *
            "which to check feasibility.",
        )
    end
    return Dict(v => value(v) for v in all_variables(model))
end

"""
    primal_feasibility_report(
        model::GenericModel{T},
        point::AbstractDict{GenericVariableRef{T},T} = _last_primal_solution(model),
        atol::T = zero(T),
        skip_missing::Bool = false,
    )::Dict{Any,T}

Given a dictionary `point`, which maps variables to primal values, return a
dictionary whose keys are the constraints with an infeasibility greater than the
supplied tolerance `atol`. The value corresponding to each key is the respective
infeasibility. Infeasibility is defined as the distance between the primal
value of the constraint (see `MOI.ConstraintPrimal`) and the nearest point by
Euclidean distance in the corresponding set.

## Notes

 * If `skip_missing = true`, constraints containing variables that are not in
   `point` will be ignored.
 * If `skip_missing = false` and a partial primal solution is provided, an error
   will be thrown.
 * If no point is provided, the primal solution from the last time the model was
   solved is used.

## Example

```jldoctest
julia> model = Model();

julia> @variable(model, 0.5 <= x <= 1);

julia> primal_feasibility_report(model, Dict(x => 0.2))
Dict{Any,Float64} with 1 entry:
  x ≥ 0.5 => 0.3
```
"""
function primal_feasibility_report(
    model::GenericModel{T},
    point::AbstractDict{GenericVariableRef{T},T} = _last_primal_solution(model);
    atol::T = zero(T),
    skip_missing::Bool = false,
) where {T}
    return primal_feasibility_report(
        model;
        atol = atol,
        skip_missing = skip_missing,
    ) do x
        fx = get(point, x, missing)
        if ismissing(fx) && !skip_missing
            error(
                "point does not contain a value for variable $x. Provide a " *
                "value, or pass `skip_missing = true`.",
            )
        end
        return fx
    end
end

"""
    primal_feasibility_report(
        point::Function,
        model::GenericModel{t};
        atol::T = zero(T),
        skip_missing::Bool = false,
    ) where {T}

A form of `primal_feasibility_report` where a function is passed as the first
argument instead of a dictionary as the second argument.

## Example

```jldoctest
julia> model = Model();

julia> @variable(model, 0.5 <= x <= 1);

julia> primal_feasibility_report(model) do v
           return value(v)
       end
Dict{Any,Float64} with 1 entry:
    x ≥ 0.5 => 0.3
```
"""
function primal_feasibility_report(
    point::Function,
    model::GenericModel{T};
    atol::T = zero(T),
    skip_missing::Bool = false,
) where {T}
    violated_constraints = Dict{Any,T}()
    for (F, S) in list_of_constraint_types(model)
        _add_infeasible_constraints(
            model,
            F,
            S,
            violated_constraints,
            point,
            atol,
        )
    end
    if num_nonlinear_constraints(model) > 0
        if skip_missing
            error(
                "`skip_missing = true` is not allowed when nonlinear " *
                "constraints are present.",
            )
        end
        _add_infeasible_nonlinear_constraints(
            model,
            violated_constraints,
            point,
            atol,
        )
    end
    return violated_constraints
end

function _add_infeasible_constraints(
    model::GenericModel{T},
    ::Type{F},
    ::Type{S},
    violated_constraints::Dict{Any,T},
    point_f::Function,
    atol::T,
) where {T,F,S}
    for con in all_constraints(model, F, S)
        obj = constraint_object(con)
        d = _distance_to_set(value.(point_f, obj.func), obj.set, T)
        if d > atol
            violated_constraints[con] = d
        end
    end
    return
end

function _add_infeasible_nonlinear_constraints(
    model::GenericModel{T},
    violated_constraints::Dict{Any,T},
    point_f::Function,
    atol::T,
) where {T}
    evaluator = NLPEvaluator(model)
    MOI.initialize(evaluator, Symbol[])
    g = zeros(num_nonlinear_constraints(model))
    MOI.eval_constraint(evaluator, g, point_f.(all_variables(model)))
    for (i, (index, constraint)) in enumerate(evaluator.model.constraints)
        d = _distance_to_set(g[i], constraint.set, T)
        if d > atol
            cref = ConstraintRef(model, index, ScalarShape())
            violated_constraints[cref] = d
        end
    end
    return
end

_distance_to_set(::Missing, set, ::Type{T}) where {T} = zero(T)
_distance_to_set(point, set, ::Type) = MOI.Utilities.distance_to_set(point, set)

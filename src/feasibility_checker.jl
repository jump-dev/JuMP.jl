#  Copyright 2017, Iain Dunning, Joey Huchette, Miles Lubin, and contributors
#  This Source Code Form is subject to the terms of the Mozilla Public
#  License, v. 2.0. If a copy of the MPL was not distributed with this
#  file, You can obtain one at https://mozilla.org/MPL/2.0/.

"""
    primal_feasibility_report(
        model::Model,
        [point::Dict{VariableRef,Float64}];
        atol::Float64 = 1e-8,
    )::Union{Nothing,Dict{Any,Float64}}

Check the primal feasibility of `model` at the point given by the dictionary
`point` mapping variables to primal values. If a variable is not given in
`point`, the value is assumed to be `0.0`.

If the point is feasible, this function returns `nothing`. If infeasible, this
function returns a dictionary mapping the constriant reference of each violated
constraint to the distance it is from being feasible.

To obtain `point` from a solution of a solved model, use:
```julia
point = Dict(v => value(v) for v in all_variables(model))
```
"""
function primal_feasibility_report(
    model::Model,
    point::Dict{VariableRef,Float64};
    atol::Float64 = 1e-8,
)
    point_f = x -> get(point, x, 0.0)
    violated_constraints = Dict{Any,Float64}()
    for (F, S) in list_of_constraint_types(model)
        # This loop is not type-stable because `F` and `S` change, but it should
        # be fine because no one should be using this in performance-critical
        # code.
        for con in all_constraints(model, F, S)
            obj = constraint_object(con)
            d = _distance_to_set(value(obj.func, point_f), obj.set)
            if d > atol
                violated_constraints[con] = d
            end
        end
    end
    return length(violated_constraints) == 0 ? nothing : violated_constraints
end

function _distance_to_set(::Any, set::MOI.AbstractSet)
    error(
        "Feasibility checker for set type $(typeof(set)) has not been " *
        "implemented yet."
    )
end

function _distance_to_set(x::Float64, set::MOI.LessThan{Float64})
    return max(x - set.upper, 0.0)
end

function _distance_to_set(x::Float64, set::MOI.GreaterThan{Float64})
    return max(set.lower - x, 0.0)
end

function _distance_to_set(x::Float64, set::MOI.EqualTo{Float64})
    return abs(set.value - x)
end

function _distance_to_set(x::Float64, set::MOI.Interval{Float64})
    return max(x - set.upper, set.lower - x, 0.0)
end

function _distance_to_set(x::Float64, ::MOI.ZeroOne)
    return min(abs(x - 0.0), abs(x - 1.0))
end

function _distance_to_set(x::Float64, ::MOI.Integer)
    return abs(x - round(Int, x))
end

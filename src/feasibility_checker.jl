#  Copyright 2017, Iain Dunning, Joey Huchette, Miles Lubin, and contributors
#  This Source Code Form is subject to the terms of the Mozilla Public
#  License, v. 2.0. If a copy of the MPL was not distributed with this
#  file, You can obtain one at https://mozilla.org/MPL/2.0/.

using LinearAlgebra

"""
    primal_feasibility_report(
        model::Model,
        point::AbstractDict{VariableRef,Float64};
        atol::Float64,
    )::Union{Nothing,Dict{Any,Float64}}

Check the primal feasibility of `model` at the point given by the dictionary
`point` mapping variables to primal values. If a variable is not given in
`point`, the value is assumed to be `0.0`.

A point is classed as feasible if the distance between the point and the set is
less than or equal to `atol`.

If the point is feasible, this function returns `nothing`. If infeasible, this
function returns a dictionary mapping the constraint reference of each violated
constraint to the distance from feasibility.

To obtain `point` from a solution of a solved model, use:
```julia
point = Dict(v => value(v) for v in all_variables(model))
```
"""
function primal_feasibility_report(
    model::Model,
    point::AbstractDict{VariableRef,Float64};
    atol::Float64,
)
    point_f = x -> get(point, x, 0.0)
    violated_constraints = Dict{Any,Float64}()
    for (F, S) in list_of_constraint_types(model)
        # This loop is not type-stable because `F` and `S` change, but it should
        # be fine because no one should be using this in performance-critical
        # code.
        for con in all_constraints(model, F, S)
            obj = constraint_object(con)
            d = _distance_to_set(value.(obj.func, point_f), obj.set)
            if d > atol
                violated_constraints[con] = d
            end
        end
    end
    if num_nl_constraints(model) > 0
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
    end
    return length(violated_constraints) == 0 ? nothing : violated_constraints
end

function _distance_to_set(::Any, set::MOI.AbstractSet)
    error(
        "Feasibility checker for set type $(typeof(set)) has not been " *
        "implemented yet."
    )
end

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


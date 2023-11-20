#  Copyright 2017, Iain Dunning, Joey Huchette, Miles Lubin, and contributors
#  This Source Code Form is subject to the terms of the Mozilla Public
#  License, v. 2.0. If a copy of the MPL was not distributed with this
#  file, You can obtain one at http://mozilla.org/MPL/2.0/.

"""
    LPMatrixData{T}

The struct returned by [`lp_matrix_data`](@ref). See [`lp_matrix_data`](@ref)
for a description of the public fields.
"""
struct LPMatrixData{T}
    A::SparseArrays.SparseMatrixCSC{T,Int}
    b_lower::Vector{T}
    b_upper::Vector{T}
    x_lower::Vector{T}
    x_upper::Vector{T}
    c::Vector{T}
    c_offset::T
    sense::MOI.OptimizationSense
    variables::Vector{GenericVariableRef{T}}
    affine_constraints::Vector{ConstraintRef}
    variable_constraints::Vector{ConstraintRef}
end

"""
    lp_matrix_data(model::GenericModel{T})

Given a JuMP model of a linear program, return an [`LPMatrixData{T}`](@ref)
struct storing data for an equivalent linear program in the form:
```math
\\begin{aligned}
\\min & c^\\top x + c_0\\
      & b_l \\le A x \\le b_u \\
      & x_l \\le x \\le x_u
\\end{aligned}
```

## Fields

The struct returned by [`lp_matrix_data`](@ref) has the fields:

 * `A::SparseArrays.SparseMatrixCSC{T,Int}`: the constraint matrix in sparse
   matrix form.
 * `b_lower::Vector{T}`: the dense vector of row lower bounds. If missing, the
   value of `typemin(T)` is used.
 * `b_upper::Vector{T}`: the dense vector of row upper bounds. If missing, the
   value of `typemax(T)` is used.
 * `x_lower::Vector{T}`: the dense vector of variable lower bounds. If missing,
   the value of `typemin(T)` is used.
 * `x_upper::Vector{T}`: the dense vector of variable upper bounds. If missing,
   the value of `typemax(T)` is used.
 * `c::Vector{T}`: the dense vector of linear objective coefficiennts
 * `c_offset::T`: the constant term in the objective function.
 * `sense::MOI.OptimizationSense`: the objective sense of the model.
 * `variables::Vector{GenericVariableRef{T}}`: a vector of [`GenericVariableRef`](@ref),
   corresponding to order of the columns in the matrix form.
 * `affine_constraints::Vector{ConstraintRef}`: a vector of [`ConstraintRef`](@ref),
   corresponding to the order of rows in the matrix form.

## Limitations

The models supported by [`lp_matrix_data`](@ref) are intentionally limited
to linear programs.

If your model has integrality, use [`relax_integrality`](@ref) to remove integer
restrictions before calling [`lp_matrix_data`](@ref).
"""
function lp_matrix_data(model::GenericModel{T}) where {T}
    variables = all_variables(model)
    columns = Dict(var => i for (i, var) in enumerate(variables))
    n = length(columns)
    cache = (;
        x_l = fill(typemin(T), n),
        x_u = fill(typemax(T), n),
        c = zeros(T, n),
        c_offset = Ref{T}(zero(T)),
        b_l = T[],
        b_u = T[],
        I = Int[],
        J = Int[],
        V = T[],
        variable_to_column = columns,
        bound_constraints = ConstraintRef[],
        affine_constraints = ConstraintRef[],
    )
    for (F, S) in list_of_constraint_types(model)
        _fill_standard_form(model, F, S, cache)
    end
    _fill_standard_form(model, objective_function_type(model), cache)
    return LPMatrixData(
        SparseArrays.sparse(cache.I, cache.J, cache.V, length(cache.b_l), n),
        cache.b_l,
        cache.b_u,
        cache.x_l,
        cache.x_u,
        cache.c,
        cache.c_offset[],
        MOI.get(model, MOI.ObjectiveSense()),
        variables,
        cache.affine_constraints,
        cache.bound_constraints,
    )
end

_bounds(s::MOI.LessThan{T}) where {T} = (typemin(T), s.upper)
_bounds(s::MOI.GreaterThan{T}) where {T} = (s.lower, typemax(T))
_bounds(s::MOI.EqualTo{T}) where {T} = (s.value, s.value)
_bounds(s::MOI.Interval{T}) where {T} = (s.lower, s.upper)

function _fill_standard_form(
    model::GenericModel{T},
    ::Type{GenericVariableRef{T}},
    ::Type{S},
    cache::Any,
) where {
    T,
    S<:Union{MOI.LessThan{T},MOI.GreaterThan{T},MOI.EqualTo{T},MOI.Interval{T}},
}
    for c in all_constraints(model, GenericVariableRef{T}, S)
        push!(cache.bound_constraints, c)
        c_obj = constraint_object(c)
        i = cache.variable_to_column[c_obj.func]
        l, u = _bounds(c_obj.set)
        cache.x_l[i] = max(cache.x_l[i], l)
        cache.x_u[i] = min(cache.x_u[i], u)
    end
    return
end

function _fill_standard_form(
    model::GenericModel{T},
    ::Type{F},
    ::Type{S},
    cache::Any,
) where {
    T,
    F<:GenericAffExpr{T,GenericVariableRef{T}},
    S<:Union{MOI.LessThan{T},MOI.GreaterThan{T},MOI.EqualTo{T},MOI.Interval{T}},
}
    for c in all_constraints(model, F, S)
        push!(cache.affine_constraints, c)
        c_obj = constraint_object(c)
        @assert iszero(c_obj.func.constant)
        row = length(cache.b_l) + 1
        l, u = _bounds(c_obj.set)
        push!(cache.b_l, l)
        push!(cache.b_u, u)
        for (x, coef) in c_obj.func.terms
            push!(cache.I, row)
            push!(cache.J, cache.variable_to_column[x])
            push!(cache.V, coef)
        end
    end
    return
end

function _fill_standard_form(
    ::GenericModel{T},
    ::Type{F},
    ::Type{S},
    ::Any,
) where {T,F,S}
    return error("Unsupported constraint type in `lp_matrix_data`: $F -in- $S")
end

function _fill_standard_form(
    model::GenericModel{T},
    ::Type{GenericVariableRef{T}},
    cache::Any,
) where {T}
    cache.c[cache.variable_to_column[objective_function(model)]] = one(T)
    cache.c_offset[] = zero(T)
    return
end

function _fill_standard_form(
    model::GenericModel{T},
    ::Type{GenericAffExpr{T,GenericVariableRef{T}}},
    cache::Any,
) where {T}
    f = objective_function(model)
    for (k, v) in f.terms
        cache.c[cache.variable_to_column[k]] += v
    end
    cache.c_offset[] = f.constant
    return
end

function _fill_standard_form(::GenericModel{T}, ::Type{F}, ::Any) where {T,F}
    return error("Unsupported objective type in `lp_matrix_data`: $F")
end

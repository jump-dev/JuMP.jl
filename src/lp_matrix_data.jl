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
    integers::Vector{Int}
    binaries::Vector{Int}
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
\\min & c^\\top x + c_0 \\\\
      & b_l \\le A x \\le b_u \\\\
      & x_l \\le x \\le x_u
\\end{aligned}
```
where elements in `x` may be continuous, integer, or binary variables.

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
 * `c::Vector{T}`: the dense vector of linear objective coefficients
 * `c_offset::T`: the constant term in the objective function.
 * `sense::MOI.OptimizationSense`: the objective sense of the model.
 * `integers::Vector{Int}`: the sorted list of column indices that are integer
   variables.
 * `binaries::Vector{Int}`: the sorted list of column indices that are binary
   variables.
 * `variables::Vector{GenericVariableRef{T}}`: a vector of [`GenericVariableRef`](@ref),
   corresponding to order of the columns in the matrix form.
 * `affine_constraints::Vector{ConstraintRef}`: a vector of [`ConstraintRef`](@ref),
   corresponding to the order of rows in the matrix form.

## Limitations

The models supported by [`lp_matrix_data`](@ref) are intentionally limited
to linear programs.

## Example

```jldoctest
julia> model = Model();

julia> @variable(model, x[1:2] >= 0);

julia> @constraint(model, x[1] + 2 * x[2] <= 1);

julia> @objective(model, Max, x[2]);

julia> data = lp_matrix_data(model);

julia> data.A
1×2 SparseArrays.SparseMatrixCSC{Float64, Int64} with 2 stored entries:
 1.0  2.0

julia> data.b_lower
1-element Vector{Float64}:
 -Inf

julia> data.b_upper
1-element Vector{Float64}:
 1.0

julia> data.x_lower
2-element Vector{Float64}:
 0.0
 0.0

julia> data.x_upper
2-element Vector{Float64}:
 Inf
 Inf

julia> data.c
2-element Vector{Float64}:
 0.0
 1.0

julia> data.c_offset
0.0

julia> data.sense
MAX_SENSE::OptimizationSense = 1
```
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
        integers = Int[],
        binaries = Int[],
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
        sort!(cache.integers),
        sort!(cache.binaries),
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
    ::Type{GenericVariableRef{T}},
    ::Type{MOI.Integer},
    cache::Any,
) where {T}
    for c in all_constraints(model, GenericVariableRef{T}, MOI.Integer)
        c_obj = constraint_object(c)
        push!(cache.integers, cache.variable_to_column[c_obj.func])
    end
    return
end

function _fill_standard_form(
    model::GenericModel{T},
    ::Type{GenericVariableRef{T}},
    ::Type{MOI.ZeroOne},
    cache::Any,
) where {T}
    for c in all_constraints(model, GenericVariableRef{T}, MOI.ZeroOne)
        c_obj = constraint_object(c)
        push!(cache.binaries, cache.variable_to_column[c_obj.func])
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

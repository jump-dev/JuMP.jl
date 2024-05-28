#  Copyright 2017, Iain Dunning, Joey Huchette, Miles Lubin, and contributors
#  This Source Code Form is subject to the terms of the Mozilla Public
#  License, v. 2.0. If a copy of the MPL was not distributed with this
#  file, You can obtain one at https://mozilla.org/MPL/2.0/.
#############################################################################
# JuMP
# An algebraic modeling language for Julia
# See https://github.com/jump-dev/JuMP.jl
#############################################################################

"""
    AbstractVectorSet

An abstract type for defining new sets in JuMP.

Implement `moi_set(::AbstractVectorSet, dim::Int)` to convert the type into an
MOI set.

See also: [`moi_set`](@ref).
"""
abstract type AbstractVectorSet end

# Used in `@variable(model, [1:n] in s)`
function build_variable(
    error_fn::Function,
    variables::Vector{<:AbstractVariable},
    set::AbstractVectorSet,
)
    return VariablesConstrainedOnCreation(
        variables,
        moi_set(set, length(variables)),
    )
end

# Used in `@constraint(model, func in set)`
function build_constraint(
    error_fn::Function,
    func::AbstractVector,
    set::AbstractVectorSet,
)
    return build_constraint(error_fn, func, moi_set(set, length(func)))
end

function build_constraint(
    error_fn::Function,
    f::AbstractVector{<:AbstractJuMPScalar},
    ::Nonnegatives,
    extra::Union{MOI.AbstractVectorSet,AbstractVectorSet},
)
    return build_constraint(error_fn, f, extra)
end

function build_constraint(
    error_fn::Function,
    f::AbstractVector{<:AbstractJuMPScalar},
    ::Nonpositives,
    extra::Union{MOI.AbstractVectorSet,AbstractVectorSet},
)
    new_f = _MA.operate!!(*, -1, f)
    return build_constraint(error_fn, new_f, extra)
end

# Handle the case `@constraint(model, X >= 0, Set())`.
function _MA.operate!!(
    ::typeof(_MA.sub_mul),
    x::AbstractArray{<:AbstractJuMPScalar},
    y::Int,
)
    if !iszero(y)
        error(
            "Operation `sub_mul` between `$(typeof(x))` and `$(typeof(y))` " *
            "is not allowed. This most often happens when you write a " *
            "constraint like `x >= y` where `x` is an array and `y` is a " *
            "constant. Use the broadcast syntax `x .- y >= 0` instead.",
        )
    end
    return x
end

# Handle the case `@constraint(model, 0 <= X, Set())`.
function _MA.operate!!(
    ::typeof(_MA.sub_mul),
    y::Int,
    x::AbstractArray{<:AbstractJuMPScalar},
)
    if !iszero(y)
        error(
            "Operation `sub_mul` between `$(typeof(y))` and `$(typeof(x))` " *
            "is not allowed. This most often happens when you write a " *
            "constraint like `x >= y` where `x` is a constant and `y` is an " *
            "array. Use the broadcast syntax `x .- y >= 0` instead.",
        )
    end
    return _MA.operate!!(*, -1, x)
end

"""
    SecondOrderCone

Second order cone object that can be used to constrain the euclidean norm of a
vector `x` to be less than or equal to a nonnegative scalar `t`. This is a
shortcut for the `MOI.SecondOrderCone`.

## Example

The following constrains ``\\|(x-1, x-2)\\|_2 \\le t`` and ``t \\ge 0``:
```jldoctest
julia> model = Model();

julia> @variable(model, x)
x

julia> @variable(model, t)
t

julia> @constraint(model, [t, x-1, x-2] in SecondOrderCone())
[t, x - 1, x - 2] ∈ MathOptInterface.SecondOrderCone(3)
```
"""
struct SecondOrderCone <: AbstractVectorSet end
moi_set(::SecondOrderCone, dim::Int) = MOI.SecondOrderCone(dim)

"""
    RotatedSecondOrderCone

Rotated second order cone object that can be used to constrain the square of the
euclidean norm of a vector `x` to be less than or equal to ``2tu`` where `t` and
`u` are nonnegative scalars. This is a shortcut for the
`MOI.RotatedSecondOrderCone`.

## Example

The following constrains ``\\|(x-1, x-2)\\|^2_2 \\le 2tx`` and ``t, x \\ge 0``:
```jldoctest
julia> model = Model();

julia> @variable(model, x)
x

julia> @variable(model, t)
t

julia> @constraint(model, [t, x, x-1, x-2] in RotatedSecondOrderCone())
[t, x, x - 1, x - 2] ∈ MathOptInterface.RotatedSecondOrderCone(4)
```
"""
struct RotatedSecondOrderCone <: AbstractVectorSet end
moi_set(::RotatedSecondOrderCone, dim::Int) = MOI.RotatedSecondOrderCone(dim)

"""
    SOS1(weights = Real[])

The SOS1 (Special Ordered Set of Type 1) set constrains a vector `x` to the set
where at most one variable can take a non-zero value, and all other elements are
zero.

The `weights` vector, if specified, induces an ordering of the variables; as
such, it should contain unique values. The `weights` vector must have the same
number of elements as the vector `x`, and the element `weights[i]` corresponds
to element `x[i]`. If not provided, the `weights` vector defaults to
`weights[i] = i`.

This is a shortcut for the [`MOI.SOS1`](@ref) set.

## Example

```jldoctest
julia> model = Model();

julia> @variable(model, x[1:3] in SOS1([4.1, 3.2, 5.0]))
3-element Vector{VariableRef}:
 x[1]
 x[2]
 x[3]

julia> print(model)
Feasibility
Subject to
 [x[1], x[2], x[3]] ∈ MathOptInterface.SOS1{Float64}([4.1, 3.2, 5.0])
```
"""
struct SOS1{T} <: AbstractVectorSet
    weights::Vector{T}
    function SOS1{T}(weights::AbstractVector = T[]) where {T}
        return new{T}(convert(Vector{T}, weights))
    end
end

SOS1(weights::AbstractVector = Int[]) = SOS1{eltype(weights)}(weights)

function moi_set(set::SOS1{T}, dim::Int) where {T}
    if length(set.weights) == 0
        return MOI.SOS1{T}(collect(1:dim))
    elseif length(set.weights) == dim
        return MOI.SOS1{T}(set.weights)
    else
        error("Weight vector in SOS1 is not of length $(dim).")
    end
end

"""
    SOS2(weights = Real[])

The SOS2 (Special Ordered Set of Type 2) set constrains a vector `x` to the set
where at most two variables can take a non-zero value, and all other elements
are zero. In addition, the two non-zero values must be consecutive given the
ordering of the `x` vector induced by `weights`.

The `weights` vector, if specified, induces an ordering of the variables; as
such, it must contain unique values. The `weights` vector must have the same
number of elements as the vector `x`, and the element `weights[i]` corresponds
to element `x[i]`. If not provided, the `weights` vector defaults to
`weights[i] = i`.

This is a shortcut for the [`MOI.SOS2`](@ref) set.

## Example

```jldoctest
julia> model = Model();

julia> @variable(model, x[1:3] in SOS2([4.1, 3.2, 5.0]))
3-element Vector{VariableRef}:
 x[1]
 x[2]
 x[3]

julia> print(model)
Feasibility
Subject to
 [x[1], x[2], x[3]] ∈ MathOptInterface.SOS2{Float64}([4.1, 3.2, 5.0])
```
"""
struct SOS2{T} <: AbstractVectorSet
    weights::Vector{T}
    function SOS2{T}(weights::AbstractVector = T[]) where {T}
        return new{T}(convert(Vector{T}, weights))
    end
end

# `Int` is chosen as a placeholder, and it is replaced by the `value_type`
# converted by `model_convert` when adding to the model.
SOS2(weights::AbstractVector = Int[]) = SOS2{eltype(weights)}(weights)

function moi_set(set::SOS2{T}, dim::Int) where {T}
    if length(set.weights) == 0
        return MOI.SOS2{T}(collect(1:dim))
    elseif length(set.weights) == dim
        return MOI.SOS2{T}(set.weights)
    else
        error("Weight vector in SOS2 is not of length $(dim).")
    end
end

"""
    AbstractScalarSet

An abstract type for defining new scalar sets in JuMP.

Implement `moi_set(::AbstractScalarSet)` to convert the type into an MOI set.

See also: [`moi_set`](@ref).
"""
abstract type AbstractScalarSet end

function build_variable(
    error_fn::Function,
    variable::AbstractVariable,
    set::AbstractScalarSet,
)
    return VariableConstrainedOnCreation(variable, moi_set(set))
end

function build_variable(
    error_fn::Function,
    variables::AbstractArray{<:AbstractVariable},
    sets::AbstractArray{<:AbstractScalarSet},
)
    return build_variable.(error_fn, variables, sets)
end

function build_variable(
    error_fn::Function,
    variables::AbstractArray{<:AbstractVariable},
    set::AbstractScalarSet,
)
    return build_variable.(error_fn, variables, Ref(set))
end

function build_constraint(
    error_fn::Function,
    func::AbstractJuMPScalar,
    set::AbstractScalarSet,
)
    return build_constraint(error_fn, func, moi_set(set))
end

"""
    Semicontinuous(lower, upper)

A short-cut for the [`MOI.Semicontinuous`](@ref) set.

This short-cut is useful because it automatically promotes `lower` and `upper`
to the same type, and converts them into the element type supported by the JuMP
model.

## Example

```jldoctest
julia> model = Model();

julia> @variable(model, x in Semicontinuous(1, 2))
x

julia> print(model)
Feasibility
Subject to
 x ∈ MathOptInterface.Semicontinuous{Int64}(1, 2)
```
"""
struct Semicontinuous{T} <: AbstractScalarSet
    lower::T
    upper::T
    function Semicontinuous(lower, upper)
        new_lower, new_upper = promote(lower, upper)
        return new{typeof(new_lower)}(new_lower, new_upper)
    end
end

function moi_set(set::Semicontinuous{T}) where {T}
    return MOI.Semicontinuous{T}(set.lower, set.upper)
end

"""
    Semiinteger(lower, upper)

A short-cut for the [`MOI.Semiinteger`](@ref) set.

This short-cut is useful because it automatically promotes `lower` and `upper`
to the same type, and converts them into the element type supported by the JuMP
model.

## Example

```jldoctest
julia> model = Model();

julia> @variable(model, x in Semiinteger(3, 5))
x

julia> print(model)
Feasibility
Subject to
 x ∈ MathOptInterface.Semiinteger{Int64}(3, 5)
```
"""
struct Semiinteger{T} <: AbstractScalarSet
    lower::T
    upper::T
    function Semiinteger(lower, upper)
        new_lower, new_upper = promote(lower, upper)
        return new{typeof(new_lower)}(new_lower, new_upper)
    end
end

function moi_set(set::Semiinteger{T}) where {T}
    return MOI.Semiinteger{T}(set.lower, set.upper)
end

"""
    Parameter(value)

A short-cut for the [`MOI.Parameter`](@ref) set.

## Example

```jldoctest
julia> model = Model();

julia> @variable(model, x in Parameter(2))
x

julia> print(model)
Feasibility
Subject to
 x ∈ MathOptInterface.Parameter{Float64}(2.0)
```
"""
struct Parameter{T} <: AbstractScalarSet
    value::T
end

function moi_set(set::Parameter{T}) where {T}
    return MOI.Parameter{T}(set.value)
end

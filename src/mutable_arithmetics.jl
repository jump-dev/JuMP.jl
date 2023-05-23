#  Copyright 2017, Iain Dunning, Joey Huchette, Miles Lubin, and contributors
#  This Source Code Form is subject to the terms of the Mozilla Public
#  License, v. 2.0. If a copy of the MPL was not distributed with this
#  file, You can obtain one at https://mozilla.org/MPL/2.0/.
#############################################################################
# JuMP
# An algebraic modeling language for Julia
# See https://github.com/jump-dev/JuMP.jl
#############################################################################
# src/mutable_arithmetics.jl
# Implements the mutable arithmetics api defined in MutableArithmetics.jl for
# `GenericAffExpr` and `GenericQuadExpr`.
#############################################################################

const _GenericAffOrQuadExpr{C,V} =
    Union{GenericAffExpr{C,V},GenericQuadExpr{C,V}}

function _change_coef(::Type{T}, ::Type{GenericAffExpr{C,V}}) where {T,C,V}
    return GenericAffExpr{T,V}
end

function _change_coef(::Type{T}, ::Type{GenericQuadExpr{C,V}}) where {T,C,V}
    return GenericQuadExpr{T,V}
end

# The default fallbacks to calling `op(zero(x), zero(y))` which produces
# allocations. The compiler could avoid these allocations at runtime with
# constant propagation as the types `x` and `y` are known at compile time, but
# apparently it does not.
function _MA.promote_operation(
    ::Union{typeof(+),typeof(-),typeof(*)},
    C::Type{<:_Constant},
    V::Type{<:AbstractVariableRef},
)
    return GenericAffExpr{_complex_convert_type(value_type(V), C),V}
end

function _MA.promote_operation(
    ::Union{typeof(+),typeof(-),typeof(*)},
    V::Type{<:AbstractVariableRef},
    C::Type{<:_Constant},
)
    return GenericAffExpr{_complex_convert_type(value_type(V), C),V}
end

function _MA.promote_operation(
    op::Union{typeof(+),typeof(-),typeof(*)},
    ::Type{S},
    A::Type{<:_GenericAffOrQuadExpr{T,V}},
) where {S<:_Constant,T,V}
    U = _MA.promote_operation(op, S, T)
    return _change_coef(U, A)
end

function _MA.promote_operation(
    op::Union{typeof(+),typeof(-),typeof(*)},
    A::Type{<:_GenericAffOrQuadExpr{T,V}},
    ::Type{S},
) where {S<:_Constant,T,V}
    U = _MA.promote_operation(op, T, S)
    return _change_coef(U, A)
end

function _MA.promote_operation(
    ::Union{typeof(+),typeof(-)},
    ::Type{V},
    ::Type{V},
) where {V<:AbstractVariableRef}
    return GenericAffExpr{value_type(V),V}
end

function _MA.promote_operation(
    ::Union{typeof(+),typeof(-)},
    ::Type{V},
    S::Type{<:_GenericAffOrQuadExpr{C,V}},
) where {C,V<:AbstractVariableRef}
    return S
end

function _MA.promote_operation(
    ::Union{typeof(+),typeof(-)},
    S::Type{<:_GenericAffOrQuadExpr{C,V}},
    ::Type{V},
) where {C,V<:AbstractVariableRef}
    return S
end

function _MA.promote_operation(
    ::Union{typeof(+),typeof(-)},
    ::Type{A},
    ::Type{A},
) where {A<:_GenericAffOrQuadExpr}
    return A
end

function _MA.promote_operation(
    ::Union{typeof(+),typeof(-)},
    ::Type{<:GenericAffExpr{T,V}},
    ::Type{<:GenericQuadExpr{T,V}},
) where {T,V}
    return GenericQuadExpr{T,V}
end

function _MA.promote_operation(
    ::Union{typeof(+),typeof(-)},
    ::Type{<:GenericQuadExpr{T,V}},
    ::Type{<:GenericAffExpr{T,V}},
) where {T,V}
    return GenericQuadExpr{T,V}
end

function _MA.promote_operation(
    ::typeof(*),
    ::Type{V},
    ::Type{V},
) where {V<:AbstractVariableRef}
    return GenericQuadExpr{value_type(V),V}
end

function _MA.promote_operation(
    ::typeof(*),
    ::Type{V},
    ::Type{GenericAffExpr{T,V}},
) where {T,V<:AbstractVariableRef}
    return GenericQuadExpr{T,V}
end

function _MA.promote_operation(
    ::typeof(*),
    ::Type{GenericAffExpr{T,V}},
    ::Type{V},
) where {T,V<:AbstractVariableRef}
    return GenericQuadExpr{T,V}
end

function _MA.isequal_canonical(x::T, y::T) where {T<:AbstractJuMPScalar}
    return isequal_canonical(x, y)
end

# `SparseArrays/src/linalg.jl` convert numbers to JuMP expressions. MA calls
# `scaling` to convert them back to numbers.
function _MA.scaling(aff::GenericAffExpr{C}) where {C}
    if !isempty(aff.terms)
        throw(InexactError("Cannot convert `$aff` to `$C`."))
    end
    return _MA.scaling(aff.constant)
end

function _MA.scaling(quad::GenericQuadExpr{C}) where {C}
    if !isempty(quad.terms)
        throw(InexactError("Cannot convert `$quad` to `$C`."))
    end
    return _MA.scaling(quad.aff)
end

_MA.mutability(::Type{<:_GenericAffOrQuadExpr}) = _MA.IsMutable()

function _MA.mutable_copy(expr::_GenericAffOrQuadExpr)
    return map_coefficients(_MA.copy_if_mutable, expr)
end

function _MA.operate!(::typeof(zero), aff::GenericAffExpr)
    empty!(aff.terms)
    aff.constant = _MA.zero!!(aff.constant)
    return aff
end

function _MA.operate!(::typeof(one), aff::GenericAffExpr)
    empty!(aff.terms)
    aff.constant = _MA.one!!(aff.constant)
    return aff
end

function _MA.operate!(
    op::Union{typeof(zero),typeof(one)},
    quad::GenericQuadExpr,
)
    _MA.operate!(op, quad.aff)
    empty!(quad.terms)
    return quad
end

function _MA.operate!(::typeof(*), expr::_GenericAffOrQuadExpr, α::_Constant)
    if iszero(α)
        return _MA.operate!(zero, expr)
    else
        return map_coefficients_inplace!(x -> _MA.mul!!(x, α), expr)
    end
end

function _MA.operate!(::typeof(+), expr::_GenericAffOrQuadExpr, x)
    return add_to_expression!(expr, x)
end

function _MA.operate!(::typeof(-), expr::_GenericAffOrQuadExpr, x)
    return add_to_expression!(expr, -1, x)
end

const _Scalar = Union{AbstractJuMPScalar,_Constant}

# `add_to_expression!` is responsible to implement all methods of up to 3 arguments.
# in addition to `add_to_expression(::GenericQuadExpr, ::Real, ::AbstractVariableRef, ::AbstractVariableRef)`.
function _MA.operate!(
    ::typeof(_MA.add_mul),
    expr::_GenericAffOrQuadExpr,
    x::_Scalar,
)
    return add_to_expression!(expr, x)
end

function _MA.operate!(
    ::typeof(_MA.add_mul),
    expr::_GenericAffOrQuadExpr,
    x::_Scalar,
    y::_Scalar,
)
    return add_to_expression!(expr, x, y)
end

function _MA.operate!(
    ::typeof(_MA.sub_mul),
    expr::_GenericAffOrQuadExpr{T},
    x::_Scalar,
) where {T}
    return add_to_expression!(expr, -one(T), x)
end

function _MA.operate!(
    ::typeof(_MA.sub_mul),
    expr::_GenericAffOrQuadExpr,
    x::_Scalar,
    y::_Scalar,
)
    return add_to_expression!(expr, -x, y)
end

# It is less costly to negate a constant than a JuMP scalar
function _MA.operate!(
    ::typeof(_MA.sub_mul),
    expr::_GenericAffOrQuadExpr,
    x::AbstractJuMPScalar,
    y::_Constant,
)
    return add_to_expression!(expr, x, -y)
end

# If `x` could be a transposed vector and `y` a vector, they are not subtypes
# of `_Scalar` but their product is.
function _MA.operate!(op::_MA.AddSubMul, expr::_GenericAffOrQuadExpr, x, y)
    return _MA.operate!(op, expr, x * y)
end

# If there are more arguments, we multiply the constants together.
@generated function _add_sub_mul_reorder!(
    op::_MA.AddSubMul,
    expr::_GenericAffOrQuadExpr,
    args::Vararg{Any,N},
) where {N}
    n = length(args)
    @assert n ≥ 3
    varidx = findall(t -> (t <: AbstractJuMPScalar), collect(args))
    allscalar = all(t -> (t <: _Constant), args[setdiff(1:n, varidx)])
    # We need to get down to only two factors
    # If there are only constants and one JuMP expressions, then we multiply
    # the constants together. Otherwise we multiply all factors except the
    # last one, there may be a better thing to do here.
    idx = (allscalar && length(varidx) == 1) ? varidx[1] : n
    coef = Expr(:call, :*, [:(args[$i]) for i in setdiff(1:n, idx)]...)
    return :(_MA.operate!(op, expr, $coef, args[$idx]))
end

function _MA.operate!(
    op::_MA.AddSubMul,
    expr::_GenericAffOrQuadExpr,
    x,
    y,
    z,
    other_args::Vararg{Any,N},
) where {N}
    return _add_sub_mul_reorder!(op, expr, x, y, z, other_args...)
end

# `add_mul(expr, args...)` defaults to `muladd(args..., expr)` which gives
# `*(args...) + expr`. If `expr isa AbstractJuMPScalar`, this reorders the terms.
# The following implementation avoids this issue and is also more efficient.
function _MA.add_mul(lhs::AbstractJuMPScalar, x::_Scalar, y::_Scalar)
    T = _MA.promote_operation(_MA.add_mul, typeof(lhs), typeof(x), typeof(y))
    expr = _MA.operate(convert, T, lhs)
    return _MA.operate!(_MA.add_mul, expr, x, y)
end

function _MA.add_mul(
    lhs::AbstractJuMPScalar,
    x::_Scalar,
    y::_Scalar,
    args::Vararg{_Scalar,N},
) where {N}
    T = _MA.promote_operation(
        _MA.add_mul,
        typeof(lhs),
        typeof(x),
        typeof(y),
        typeof.(args)...,
    )
    expr = _MA.operate(convert, T, lhs)
    return _MA.operate!(_MA.add_mul, expr, x, y, args...)
end

function _MA.sub_mul(lhs::AbstractJuMPScalar, x::_Scalar, y::_Scalar)
    T = _MA.promote_operation(_MA.sub_mul, typeof(lhs), typeof(x), typeof(y))
    expr = _MA.operate(convert, T, lhs)
    return _MA.operate!(_MA.sub_mul, expr, x, y)
end

function _MA.sub_mul(
    lhs::AbstractJuMPScalar,
    x::_Scalar,
    y::_Scalar,
    args::Vararg{_Scalar,N},
) where {N}
    T = _MA.promote_operation(
        _MA.sub_mul,
        typeof(lhs),
        typeof(x),
        typeof(y),
        typeof.(args)...,
    )
    expr = _MA.operate(convert, T, lhs)
    return _MA.operate!(_MA.sub_mul, expr, x, y, args...)
end

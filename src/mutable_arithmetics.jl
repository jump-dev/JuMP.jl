#  Copyright 2017, Iain Dunning, Joey Huchette, Miles Lubin, and contributors
#  This Source Code Form is subject to the terms of the Mozilla Public
#  License, v. 2.0. If a copy of the MPL was not distributed with this
#  file, You can obtain one at http://mozilla.org/MPL/2.0/.
#############################################################################
# JuMP
# An algebraic modeling language for Julia
# See http://github.com/JuliaOpt/JuMP.jl
#############################################################################
# src/mutable_arithmetics.jl
# Implements the mutable arithmetics api defined in MutableArithmetics.jl for
# `GenericAffExpr` and `GenericQuadExpr`.
#############################################################################

const _GenericAffOrQuadExpr{C, V} = Union{GenericAffExpr{C, V}, GenericQuadExpr{C, V}}

_MA.isequal_canonical(x::T, y::T) where {T<:AbstractJuMPScalar} = isequal_canonical(x, y)

# `SparseArrays/src/linalg.jl` convert numbers to JuMP expressions. MA calls
# `scaling` to convert them back to numbers.
function _MA.scaling(aff::GenericAffExpr{C}) where C
    if !isempty(aff.terms)
        throw(InexactError("Cannot convert `$aff` to `$C`."))
    end
    return _MA.scaling(aff.constant)
end
function _MA.scaling(quad::GenericQuadExpr{C}) where C
    if !isempty(quad.terms)
        throw(InexactError("Cannot convert `$quad` to `$C`."))
    end
    return _MA.scaling(quad.aff)
end

_MA.mutability(::Type{<:_GenericAffOrQuadExpr}) = _MA.IsMutable()
function _MA.mutable_copy(expr::_GenericAffOrQuadExpr)
    return map_coefficients(_MA.copy_if_mutable, expr)
end

function _MA.mutable_operate!(::typeof(zero), aff::GenericAffExpr)
    empty!(aff.terms)
    aff.constant = _MA.zero!(aff.constant)
    return aff
end
function _MA.mutable_operate!(::typeof(one), aff::GenericAffExpr)
    empty!(aff.terms)
    aff.constant = _MA.one!(aff.constant)
    return aff
end
function _MA.mutable_operate!(op::Union{typeof(zero), typeof(one)}, quad::GenericQuadExpr)
    _MA.mutable_operate!(op, quad.aff)
    empty!(quad.terms)
    return quad
end

function _MA.mutable_operate!(::typeof(*), expr::_GenericAffOrQuadExpr, α::_Constant)
    if iszero(α)
        return _MA.mutable_operate!(zero, expr)
    else
        return map_coefficients_inplace!(x -> _MA.mul!(x, α), expr)
    end
end

function _MA.mutable_operate!(::typeof(+), expr::_GenericAffOrQuadExpr, x)
    return add_to_expression!(expr, x)
end
function _MA.mutable_operate!(::typeof(-), expr::_GenericAffOrQuadExpr, x)
    return add_to_expression!(expr, -1, x)
end

# `add_to_expression!` is responsible to implement all methods of up to 3 arguments.
# in addition to `add_to_expression(::GenericQuadExpr, ::Real, ::AbstractVariableRef, ::AbstractVariableRef)`.
function _MA.mutable_operate!(::typeof(_MA.add_mul), expr::_GenericAffOrQuadExpr, x)
    return add_to_expression!(expr, x)
end
function _MA.mutable_operate!(::typeof(_MA.add_mul), expr::_GenericAffOrQuadExpr, x, y)
    return add_to_expression!(expr, x, y)
end
# If there are more arguments, we multiply the rest to make the number of argument drop.
function _MA.mutable_operate!(op::typeof(_MA.add_mul), expr::_GenericAffOrQuadExpr, x, y, z, args::Vararg{Any, N}) where N
    return _MA.mutable_operate!(op, expr, x, *(y, z, args...))
end
# This method is missing for both affine and quadratic expressions.
function add_to_expression!(expr::_GenericAffOrQuadExpr, α::_Constant, β::_Constant)
    return add_to_expression!(expr, *(α, β))
end

const _AffineLike = Union{AbstractVariableRef, GenericAffExpr, _Constant}

const _Scalar = Union{AbstractJuMPScalar, _Constant}

# `add_mul(expr, args...)` defaults to `muladd(args..., expr)` which gives
# `*(args...) + expr`. If `expr isa AbstractJuMPScalar`, this reorders the terms.
# The following implementation avoids this issue and is also more efficient.
function _MA.add_mul(lhs::AbstractJuMPScalar, x::_Scalar, y::_Scalar)
    T = _MA.promote_operation(_MA.add_mul, typeof(lhs), typeof(x), typeof(y))
    return _add_mul_to_type(T, lhs, x, y)
end
function _MA.add_mul(lhs::AbstractJuMPScalar, x::_Scalar, y::_Scalar, args::Vararg{_Scalar, N}) where N
    T = _MA.promote_operation(_MA.add_mul, typeof(lhs), typeof(x), typeof(y), typeof.(args)...)
    return _add_mul_to_type(T, lhs, x, y, args...)
end
function _add_mul_to_type(::Type{T}, expr::T, factors::Vararg{Any, N}) where {T, N}
    return _MA.mutable_operate!(_MA.add_mul, expr, factors...)
end
function _add_mul_to_type(::Type{T}, lhs, factors::Vararg{Any, N}) where {T, N}
    expr = zero(T)
    _MA.mutable_operate!(+, expr, lhs)
    return _MA.mutable_operate!(_MA.add_mul, expr, factors...)
end

#############################################################################
# Helpers to initialize memory for GenericAffExpr/GenericQuadExpr
#############################################################################

# TODO Add it to MA API to allow linear algebra generic implementation to exploit this

#_sizehint_expr!(a::GenericAffExpr, n::Int) = sizehint!(a, n)
#
## TODO: Why do we allocate the same size for the quadratic and affine parts?
#function _sizehint_expr!(q::GenericQuadExpr, n::Int)
#        sizehint!(q.terms,  length(q.terms) + n)
#        _sizehint_expr!(q.aff, n)
#        nothing
#end
#_sizehint_expr!(q, n) = nothing

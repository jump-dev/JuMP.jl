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

const GenericAffOrQuadExpr{C, V} = Union{GenericAffExpr{C, V}, GenericQuadExpr{C, V}}

MA.isequal_canonical(x::T, y::T) where {T<:AbstractJuMPScalar} = isequal_canonical(x, y)

# `SparseArrays/src/linalg.jl` convert numbers to JuMP expressions. MA calls
# `scaling` to convert them back to numbers.
function MA.scaling(aff::GenericAffExpr{C}) where C
    if !isempty(aff.terms)
        throw(InexactError("Cannot convert `$aff` to `$C`."))
    end
    return MA.scaling(aff.constant)
end
function MA.scaling(quad::GenericQuadExpr{C}) where C
    if !isempty(quad.terms)
        throw(InexactError("Cannot convert `$quad` to `$C`."))
    end
    return MA.scaling(quad.aff)
end

MA.mutability(::Type{<:GenericAffOrQuadExpr}) = MA.IsMutable()
function MA.mutable_copy(expr::GenericAffOrQuadExpr)
    return map_coefficients(MA.copy_if_mutable, expr)
end

function MA.mutable_operate!(::typeof(zero), aff::GenericAffExpr)
    empty!(aff.terms)
    aff.constant = MA.zero!(aff.constant)
    return aff
end
function MA.mutable_operate!(::typeof(one), aff::GenericAffExpr)
    empty!(aff.terms)
    aff.constant = MA.one!(aff.constant)
    return aff
end
function MA.mutable_operate!(op::Union{typeof(zero), typeof(one)}, quad::GenericQuadExpr)
    MA.mutable_operate!(op, quad.aff)
    empty!(quad.terms)
    return quad
end

function MA.mutable_operate!(::typeof(*), expr::GenericAffOrQuadExpr, α::Constant)
    if iszero(α)
        return MA.mutable_operate!(zero, expr)
    else
        return map_coefficients_inplace!(x -> MA.mul!(x, α), expr)
    end
end

function MA.mutable_operate!(::typeof(+), expr::GenericAffOrQuadExpr, x)
    return add_to_expression!(expr, x)
end
function MA.mutable_operate!(::typeof(-), expr::GenericAffOrQuadExpr, x)
    return add_to_expression!(expr, -1, x)
end

# `add_to_expression!` is responsible to implement all methods of up to 3 arguments.
# in addition to `add_to_expression(::GenericQuadExpr, ::Real, ::AbstractVariableRef, ::AbstractVariableRef)`.
function MA.mutable_operate!(::typeof(MA.add_mul), expr::GenericAffOrQuadExpr, x)
    return add_to_expression!(expr, x)
end
function MA.mutable_operate!(::typeof(MA.add_mul), expr::GenericAffOrQuadExpr, x, y)
    return add_to_expression!(expr, x, y)
end
# If there are more arguments, we multiply the rest to make the number of argument drop.
function MA.mutable_operate!(op::typeof(MA.add_mul), expr::GenericAffOrQuadExpr, x, y, z, args::Vararg{Any, N}) where N
    return MA.mutable_operate!(op, expr, x, *(y, z, args...))
end
# This method is missing for both affine and quadratic expressions.
function add_to_expression!(expr::GenericAffOrQuadExpr, α::Constant, β::Constant)
    return add_to_expression!(expr, *(α, β))
end

const _AffineLike = Union{AbstractVariableRef, GenericAffExpr, Constant}

const _Scalar = Union{AbstractJuMPScalar, Constant}

# `add_mul(expr, args...)` defaults to `muladd(args..., expr)` which gives
# `*(args...) + expr`. If `expr isa AbstractJuMPScalar`, this reorders the terms.
# The following implementation avoids this issue and is also more efficient.
function MA.add_mul(lhs::AbstractJuMPScalar, x::_Scalar, y::_Scalar)
    T = MA.promote_operation(MA.add_mul, typeof(lhs), typeof(x), typeof(y))
    return _add_mul_to_type(T, lhs, x, y)
end
function MA.add_mul(lhs::AbstractJuMPScalar, x::_Scalar, y::_Scalar, args::Vararg{_Scalar, N}) where N
    T = MA.promote_operation(MA.add_mul, typeof(lhs), typeof(x), typeof(y), typeof.(args)...)
    return _add_mul_to_type(T, lhs, x, y, args...)
end
function _add_mul_to_type(::Type{T}, expr::T, factors::Vararg{Any, N}) where {T, N}
    return MA.mutable_operate!(MA.add_mul, expr, factors...)
end
function _add_mul_to_type(::Type{T}, lhs, factors::Vararg{Any, N}) where {T, N}
    expr = zero(T)
    MA.mutable_operate!(+, expr, lhs)
    return MA.mutable_operate!(MA.add_mul, expr, factors...)
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

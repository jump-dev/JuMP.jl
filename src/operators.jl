#  Copyright 2017, Iain Dunning, Joey Huchette, Miles Lubin, and contributors
#  This Source Code Form is subject to the terms of the Mozilla Public
#  License, v. 2.0. If a copy of the MPL was not distributed with this
#  file, You can obtain one at http://mozilla.org/MPL/2.0/.
#############################################################################
# JuMP
# An algebraic modeling language for Julia
# See http://github.com/JuliaOpt/JuMP.jl
#############################################################################

const _JuMPTypes = Union{AbstractJuMPScalar, NonlinearExpression}
_float(x::Number) = convert(Float64, x)
_float(J::UniformScaling) = _float(J.λ)

# Overloads
#
# Different objects that must all interact:
# 1. _Constant
# 2. AbstractVariableRef
# 4. GenericAffExpr
# 5. GenericQuadExpr

# _Constant
# _Constant--_Constant obviously already taken care of!
# _Constant--VariableRef
Base.:+(lhs::_Constant, rhs::AbstractVariableRef) = _build_aff_expr(_float(lhs), 1.0, rhs)
Base.:-(lhs::_Constant, rhs::AbstractVariableRef) = _build_aff_expr(_float(lhs), -1.0, rhs)
function Base.:*(lhs::_Constant, rhs::AbstractVariableRef)
    if iszero(lhs)
        return zero(GenericAffExpr{Float64, typeof(rhs)})
    else
        return _build_aff_expr(0.0, _float(lhs), rhs)
    end
end
# _Constant--_GenericAffOrQuadExpr
function Base.:+(lhs::_Constant, rhs::_GenericAffOrQuadExpr)
    result = _MA.mutable_copy(rhs)
    add_to_expression!(result, lhs)
    return result
end
function Base.:-(lhs::_Constant, rhs::_GenericAffOrQuadExpr)
    result = -rhs
    add_to_expression!(result, lhs)
    return result
end
function Base.:*(lhs::_Constant, rhs::_GenericAffOrQuadExpr)
    if iszero(lhs)
        return zero(rhs)
    else
        α = _constant_to_number(lhs)
        return map_coefficients(c -> α * c, rhs)
    end
end

# AbstractVariableRef (or, AbstractJuMPScalar)
# TODO: What is the role of AbstractJuMPScalar??
Base.:+(lhs::AbstractJuMPScalar) = lhs
Base.:-(lhs::AbstractVariableRef) = _build_aff_expr(0.0, -1.0, lhs)
Base.:*(lhs::AbstractJuMPScalar) = lhs # make this more generic so extensions don't have to define unary multiplication for our macros
# AbstractVariableRef--_Constant
Base.:+(lhs::AbstractVariableRef, rhs::_Constant) = (+)( rhs,lhs)
Base.:-(lhs::AbstractVariableRef, rhs::_Constant) = (+)(-rhs,lhs)
Base.:*(lhs::AbstractVariableRef, rhs::_Constant) = (*)(rhs,lhs)
Base.:/(lhs::AbstractVariableRef, rhs::_Constant) = (*)(1.0/rhs,lhs)
# AbstractVariableRef--AbstractVariableRef
Base.:+(lhs::V, rhs::V) where {V <: AbstractVariableRef} = _build_aff_expr(0.0, 1.0, lhs, 1.0, rhs)
function Base.:-(lhs::V, rhs::V) where {V <: AbstractVariableRef}
    if lhs == rhs
        return zero(GenericAffExpr{Float64, V})
    else
        return _build_aff_expr(0.0, 1.0, lhs, -1.0, rhs)
    end
end
function Base.:*(lhs::V, rhs::V) where {V <: AbstractVariableRef}
    GenericQuadExpr(GenericAffExpr{Float64, V}(), UnorderedPair(lhs, rhs) => 1.0)
end
# AbstractVariableRef--GenericAffExpr
function Base.:+(lhs::V, rhs::GenericAffExpr{C, V}) where {C, V <: AbstractVariableRef}
    # For the variables to have the proper order in the result, we need to add the lhs first.
    result = zero(rhs)
    result.constant = rhs.constant
    sizehint!(result, length(linear_terms(rhs)) + 1)
    add_to_expression!(result, one(C), lhs)
    for (coef, var) in linear_terms(rhs)
        add_to_expression!(result, coef, var)
    end
    return result
end

function Base.:-(lhs::V, rhs::GenericAffExpr{C,V}) where {C,V <: AbstractVariableRef}
    # For the variables to have the proper order in the result, we need to add the lhs first.
    result = zero(rhs)
    result.constant = -rhs.constant
    sizehint!(result, length(linear_terms(rhs)) + 1)
    add_to_expression!(result, one(C), lhs)
    for (coef, var) in linear_terms(rhs)
        add_to_expression!(result, -coef, var)
    end
    return result
end

function Base.:*(lhs::V, rhs::GenericAffExpr{C, V}) where {C, V <: AbstractVariableRef}
    if !iszero(rhs.constant)
        result = GenericQuadExpr{C, V}(lhs*rhs.constant)
    else
        result = zero(GenericQuadExpr{C, V})
    end
    for (coef, var) in linear_terms(rhs)
        add_to_expression!(result, coef, lhs, var)
    end
    return result
end
Base.:/(lhs::AbstractVariableRef, rhs::GenericAffExpr) = error("Cannot divide a variable by an affine expression")
# AbstractVariableRef--GenericQuadExpr
Base.:+(v::AbstractVariableRef, q::GenericQuadExpr) = GenericQuadExpr(v+q.aff, copy(q.terms))
function Base.:-(v::AbstractVariableRef, q::GenericQuadExpr)
    result = -q
    # This makes an unnecessary copy of aff, but it's important for v to appear
    # first.
    result.aff = v + result.aff
    return result
end

# GenericAffExpr
Base.:+(lhs::GenericAffExpr) = lhs
Base.:-(lhs::GenericAffExpr) = map_coefficients(-, lhs)
# GenericAffExpr--_Constant
Base.:+(lhs::GenericAffExpr, rhs::_Constant) = (+)(rhs,lhs)
Base.:-(lhs::GenericAffExpr, rhs::_Constant) = (+)(-rhs,lhs)
Base.:*(lhs::GenericAffExpr, rhs::_Constant) = (*)(rhs,lhs)
Base.:/(lhs::GenericAffExpr, rhs::_Constant) = map_coefficients(c -> c/rhs, lhs)
function Base.:^(lhs::Union{AbstractVariableRef, GenericAffExpr}, rhs::Integer)
    if rhs == 2
        return lhs*lhs
    elseif rhs == 1
        return convert(GenericQuadExpr{Float64, variable_ref_type(lhs)}, lhs)
    elseif rhs == 0
        return one(GenericQuadExpr{Float64, variable_ref_type(lhs)})
    else
        error("Only exponents of 0, 1, or 2 are currently supported. Are you trying to build a nonlinear problem? Make sure you use @NLconstraint/@NLobjective.")
    end
end
Base.:^(lhs::Union{AbstractVariableRef, GenericAffExpr}, rhs::_Constant) = error("Only exponents of 0, 1, or 2 are currently supported. Are you trying to build a nonlinear problem? Make sure you use @NLconstraint/@NLobjective.")
# GenericAffExpr--AbstractVariableRef
function Base.:+(lhs::GenericAffExpr{C, V}, rhs::V) where {C, V <: AbstractVariableRef}
    return add_to_expression!(copy(lhs), one(C), rhs)
end
function Base.:-(lhs::GenericAffExpr{C, V}, rhs::V) where {C, V <: AbstractVariableRef}
    return add_to_expression!(copy(lhs), -one(C), rhs)
end
# Don't fall back on AbstractVariableRef*GenericAffExpr to preserve lhs/rhs
# consistency (appears in printing).
function Base.:*(lhs::GenericAffExpr{C, V}, rhs::V) where {C, V <: AbstractVariableRef}
    if !iszero(lhs.constant)
        result = GenericQuadExpr{C, V}(lhs.constant*rhs)
    else
        result = zero(GenericQuadExpr{C, V})
    end
    for (coef, var) in linear_terms(lhs)
        add_to_expression!(result, coef, var, rhs)
    end
    return result
end
Base.:/(lhs::GenericAffExpr, rhs::AbstractVariableRef) = error("Cannot divide affine expression by a variable")
# AffExpr--AffExpr
function Base.:+(lhs::GenericAffExpr{C, V}, rhs::GenericAffExpr{C, V}) where {C, V<:_JuMPTypes}
    if length(linear_terms(lhs)) > 50 || length(linear_terms(rhs)) > 50
        if length(linear_terms(lhs)) > 1
            operator_warn(owner_model(first(linear_terms(lhs))[2]))
        end
    end
    return add_to_expression!(copy(lhs), rhs)
end

function Base.:-(lhs::GenericAffExpr{C, V}, rhs::GenericAffExpr{C, V}) where {C, V<:_JuMPTypes}
    result = copy(lhs)
    result.constant -= rhs.constant
    sizehint!(result, length(linear_terms(lhs)) + length(linear_terms(rhs)))
    for (coef, var) in linear_terms(rhs)
        add_to_expression!(result, -coef, var)
    end
    return result
end

function Base.:*(lhs::GenericAffExpr{C, V}, rhs::GenericAffExpr{C, V}) where {C, V<:_JuMPTypes}
    result = zero(GenericQuadExpr{C, V})
    add_to_expression!(result, lhs, rhs)
    return result
end
# GenericAffExpr--GenericQuadExpr
Base.:+(a::GenericAffExpr, q::GenericQuadExpr) = GenericQuadExpr(a+q.aff, copy(q.terms))
function Base.:-(a::GenericAffExpr, q::GenericQuadExpr)
    result = -q
    # This makes an unnecessary copy of aff, but it's important for a to appear
    # first.
    result.aff = a + result.aff
    return result
end

# GenericQuadExpr
Base.:+(lhs::GenericQuadExpr) = lhs
Base.:-(lhs::GenericQuadExpr) = map_coefficients(-, lhs)
# GenericQuadExpr--_Constant
# We don't do `+rhs` as `UniformScaling` does not support unary `+`
Base.:+(lhs::GenericQuadExpr, rhs::_Constant) = (+)(rhs, lhs)
Base.:-(lhs::GenericQuadExpr, rhs::_Constant) = (+)(-rhs, lhs)
Base.:*(lhs::GenericQuadExpr, rhs::_Constant) = (*)(rhs, lhs)
Base.:/(lhs::GenericQuadExpr, rhs::_Constant) = (*)(inv(rhs), lhs)
# GenericQuadExpr--AbstractVariableRef
Base.:+(q::GenericQuadExpr, v::AbstractVariableRef) = GenericQuadExpr(q.aff+v, copy(q.terms))
Base.:-(q::GenericQuadExpr, v::AbstractVariableRef) = GenericQuadExpr(q.aff-v, copy(q.terms))
Base.:*(q::GenericQuadExpr, v::AbstractVariableRef) = error("Cannot multiply a quadratic expression by a variable")
Base.:/(q::GenericQuadExpr, v::AbstractVariableRef) = error("Cannot divide a quadratic expression by a variable")
# GenericQuadExpr--GenericAffExpr
Base.:+(q::GenericQuadExpr, a::GenericAffExpr) = GenericQuadExpr(q.aff+a, copy(q.terms))
Base.:-(q::GenericQuadExpr, a::GenericAffExpr) = GenericQuadExpr(q.aff-a, copy(q.terms))
Base.:*(q::GenericQuadExpr, a::GenericAffExpr) = error("Cannot multiply a quadratic expression by an aff. expression")
Base.:/(q::GenericQuadExpr, a::GenericAffExpr) = error("Cannot divide a quadratic expression by an aff. expression")
# GenericQuadExpr--GenericQuadExpr
function Base.:+(q1::GenericQuadExpr, q2::GenericQuadExpr)
    result = copy(q1)
    for (coef, var1, var2) in quad_terms(q2)
        add_to_expression!(result, coef, var1, var2)
    end
    for (coef, var) in linear_terms(q2)
        add_to_expression!(result, coef, var)
    end
    result.aff.constant += q2.aff.constant
    return result
end
function Base.:-(q1::GenericQuadExpr, q2::GenericQuadExpr)
    result = copy(q1)
    for (coef, var1, var2) in quad_terms(q2)
        add_to_expression!(result, -coef, var1, var2)
    end
    for (coef, var) in linear_terms(q2)
        add_to_expression!(result, -coef, var)
    end
    result.aff.constant -= q2.aff.constant
    return result
end

Base.:(==)(lhs::GenericAffExpr, rhs::GenericAffExpr) = (lhs.terms == rhs.terms) && (lhs.constant == rhs.constant)
Base.:(==)(lhs::GenericQuadExpr, rhs::GenericQuadExpr) = (lhs.terms == rhs.terms) && (lhs.aff == rhs.aff)

# Base Julia's generic fallback vecdot, aka dot, requires that dot, aka LinearAlgebra.dot, be defined
# for scalars, so instead of defining them one-by-one, we will
# fallback to the multiplication operator
LinearAlgebra.dot(lhs::_JuMPTypes, rhs::_JuMPTypes) = lhs*rhs
LinearAlgebra.dot(lhs::_JuMPTypes, rhs::_Constant) = lhs*rhs
LinearAlgebra.dot(lhs::_Constant, rhs::_JuMPTypes) = lhs*rhs

Base.promote_rule(V::Type{<:AbstractVariableRef}, R::Type{<:Real}) = GenericAffExpr{Float64, V}
Base.promote_rule(V::Type{<:AbstractVariableRef}, ::Type{<:GenericAffExpr{T}}) where {T} = GenericAffExpr{T, V}
Base.promote_rule(V::Type{<:AbstractVariableRef}, ::Type{<:GenericQuadExpr{T}}) where {T} = GenericQuadExpr{T, V}
Base.promote_rule(::Type{GenericAffExpr{S, V}}, R::Type{<:Real}) where {S, V} = GenericAffExpr{promote_type(S, R), V}
Base.promote_rule(::Type{<:GenericAffExpr{S, V}}, ::Type{<:GenericQuadExpr{T, V}}) where {S, T, V} = GenericQuadExpr{promote_type(S, T), V}
Base.promote_rule(::Type{GenericQuadExpr{S, V}}, R::Type{<:Real}) where {S, V} = GenericQuadExpr{promote_type(S, R), V}

Base.transpose(x::AbstractJuMPScalar) = x

# only real scalars are supported
Base.conj(x::AbstractJuMPScalar) = x

# Can remove the following code once == overloading is removed

function LinearAlgebra.issymmetric(x::Matrix{T}) where {T <: _JuMPTypes}
    (n = size(x,1)) == size(x,2) || return false
    for i in 1:n, j in (i+1):n
        isequal(x[i,j], x[j,i]) || return false
    end
    true
end

###############################################################################
# nonlinear function fallbacks for JuMP built-in types
###############################################################################

const op_hint = "Are you trying to build a nonlinear problem? Make sure you use @NLconstraint/@NLobjective."
for (func,_) in Calculus.symbolic_derivatives_1arg(), typ in [:AbstractVariableRef,:GenericAffExpr,:GenericQuadExpr]
    errstr = "$func is not defined for type $typ. $op_hint"
    if isdefined(Base, func)
        @eval Base.$(func)(::$typ) = error($errstr)
    end
end

Base.:*(::T, ::S) where {T <: GenericQuadExpr, S <: Union{AbstractVariableRef, GenericAffExpr, GenericQuadExpr}} =
    error( "*(::$T,::$S) is not defined. $op_hint")
Base.:*(lhs::GenericQuadExpr, rhs::GenericQuadExpr) =
    error( "*(::GenericQuadExpr,::GenericQuadExpr) is not defined. $op_hint")
Base.:*(::S, ::T) where {T <: GenericQuadExpr,
                         S <: Union{AbstractVariableRef, GenericAffExpr, GenericQuadExpr}} =
    error( "*(::$S,::$T) is not defined. $op_hint")
Base.:/(::S, ::T) where {S <: Union{_Constant, AbstractVariableRef, GenericAffExpr, GenericQuadExpr},
                         T <: Union{AbstractVariableRef, GenericAffExpr, GenericQuadExpr}} =
    error( "/(::$S,::$T) is not defined. $op_hint")

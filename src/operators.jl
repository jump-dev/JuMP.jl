#  Copyright 2017, Iain Dunning, Joey Huchette, Miles Lubin, and contributors
#  This Source Code Form is subject to the terms of the Mozilla Public
#  License, v. 2.0. If a copy of the MPL was not distributed with this
#  file, You can obtain one at https://mozilla.org/MPL/2.0/.
#############################################################################
# JuMP
# An algebraic modeling language for Julia
# See https://github.com/jump-dev/JuMP.jl
#############################################################################

const _JuMPTypes = Union{AbstractJuMPScalar,NonlinearExpression}

_float_type(::Type{<:Real}) = Float64
_float_type(::Type{UniformScaling{T}}) where {T} = _float_type(T)
_float_type(::Type{<:Complex}) = Complex{Float64}

_float(x::Real) = convert(Float64, x)
_float(x::Complex) = convert(Complex{Float64}, x)
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
function Base.:+(lhs::_Constant, rhs::AbstractVariableRef)
    constant = _float(lhs)
    return _build_aff_expr(constant, one(constant), rhs)
end
function Base.:-(lhs::_Constant, rhs::AbstractVariableRef)
    constant = _float(lhs)
    return _build_aff_expr(constant, -one(constant), rhs)
end
function Base.:*(lhs::_Constant, rhs::AbstractVariableRef)
    coef = _float(lhs)
    if iszero(coef)
        return zero(GenericAffExpr{typeof(coef),typeof(rhs)})
    else
        return _build_aff_expr(zero(coef), coef, rhs)
    end
end
# _Constant--_GenericAffOrQuadExpr
function Base.:+(lhs::_Constant, rhs::_GenericAffOrQuadExpr)
    # If `lhs` is complex and `rhs` has real coefficients then the conversion is needed
    T = _MA.promote_operation(+, _float_type(typeof(lhs)), typeof(rhs))
    result = _MA.mutable_copy(convert(T, rhs))
    add_to_expression!(result, lhs)
    return result
end
function Base.:-(lhs::_Constant, rhs::_GenericAffOrQuadExpr)
    # If `lhs` is complex and `rhs` has real coefficients then the conversion is needed
    T = _MA.promote_operation(+, _float_type(typeof(lhs)), typeof(rhs))
    result = convert(T, -rhs)
    add_to_expression!(result, lhs)
    return result
end
function Base.:*(lhs::_Constant, rhs::_GenericAffOrQuadExpr)
    if iszero(lhs)
        # If `lhs` is complex and `rhs` has real coefficients, `zero(rhs)` would not work
        return zero(
            _MA.promote_operation(*, _float_type(typeof(lhs)), typeof(rhs)),
        )
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
Base.:+(lhs::AbstractVariableRef, rhs::_Constant) = (+)(rhs, lhs)
Base.:-(lhs::AbstractVariableRef, rhs::_Constant) = (+)(-rhs, lhs)
Base.:*(lhs::AbstractVariableRef, rhs::_Constant) = (*)(rhs, lhs)
Base.:/(lhs::AbstractVariableRef, rhs::_Constant) = (*)(1.0 / rhs, lhs)
# AbstractVariableRef--AbstractVariableRef
function Base.:+(lhs::V, rhs::V) where {V<:AbstractVariableRef}
    return _build_aff_expr(0.0, 1.0, lhs, 1.0, rhs)
end
function Base.:-(lhs::V, rhs::V) where {V<:AbstractVariableRef}
    if lhs == rhs
        return zero(GenericAffExpr{Float64,V})
    else
        return _build_aff_expr(0.0, 1.0, lhs, -1.0, rhs)
    end
end
function Base.:*(lhs::V, rhs::V) where {V<:AbstractVariableRef}
    return GenericQuadExpr(
        GenericAffExpr{Float64,V}(),
        UnorderedPair(lhs, rhs) => 1.0,
    )
end
# AbstractVariableRef--GenericAffExpr
function Base.:+(
    lhs::V,
    rhs::GenericAffExpr{C,V},
) where {C,V<:AbstractVariableRef}
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

function Base.:-(
    lhs::V,
    rhs::GenericAffExpr{C,V},
) where {C,V<:AbstractVariableRef}
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

function Base.:*(
    lhs::V,
    rhs::GenericAffExpr{C,V},
) where {C,V<:AbstractVariableRef}
    if !iszero(rhs.constant)
        result = GenericQuadExpr{C,V}(lhs * rhs.constant)
    else
        result = zero(GenericQuadExpr{C,V})
    end
    for (coef, var) in linear_terms(rhs)
        add_to_expression!(result, coef, lhs, var)
    end
    return result
end
function Base.:/(lhs::AbstractVariableRef, rhs::GenericAffExpr)
    return error("Cannot divide a variable by an affine expression")
end
# AbstractVariableRef--GenericQuadExpr
function Base.:+(v::AbstractVariableRef, q::GenericQuadExpr)
    return GenericQuadExpr(v + q.aff, copy(q.terms))
end
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
Base.:+(lhs::GenericAffExpr, rhs::_Constant) = (+)(rhs, lhs)
Base.:-(lhs::GenericAffExpr, rhs::_Constant) = (+)(-rhs, lhs)
Base.:*(lhs::GenericAffExpr, rhs::_Constant) = (*)(rhs, lhs)
function Base.:/(lhs::GenericAffExpr, rhs::_Constant)
    return map_coefficients(c -> c / rhs, lhs)
end

function Base.:^(lhs::AbstractVariableRef, rhs::Integer)
    if rhs == 2
        return lhs * lhs
    elseif rhs == 1
        return convert(GenericQuadExpr{Float64,variable_ref_type(lhs)}, lhs)
    elseif rhs == 0
        return one(GenericQuadExpr{Float64,variable_ref_type(lhs)})
    else
        error(
            "Only exponents of 0, 1, or 2 are currently supported. Are you " *
            "trying to build a nonlinear problem? Make sure you use " *
            "@NLconstraint/@NLobjective.",
        )
    end
end

function Base.:^(lhs::GenericAffExpr{T}, rhs::Integer) where {T}
    if rhs == 2
        return lhs * lhs
    elseif rhs == 1
        return convert(GenericQuadExpr{T,variable_ref_type(lhs)}, lhs)
    elseif rhs == 0
        return one(GenericQuadExpr{T,variable_ref_type(lhs)})
    else
        error(
            "Only exponents of 0, 1, or 2 are currently supported. Are you " *
            "trying to build a nonlinear problem? Make sure you use " *
            "@NLconstraint/@NLobjective.",
        )
    end
end

function Base.:^(lhs::Union{AbstractVariableRef,GenericAffExpr}, rhs::_Constant)
    return error(
        "Only exponents of 0, 1, or 2 are currently supported. Are you trying to build a nonlinear problem? Make sure you use @NLconstraint/@NLobjective.",
    )
end
# GenericAffExpr--AbstractVariableRef
function Base.:+(
    lhs::GenericAffExpr{C,V},
    rhs::V,
) where {C,V<:AbstractVariableRef}
    return add_to_expression!(copy(lhs), one(C), rhs)
end
function Base.:-(
    lhs::GenericAffExpr{C,V},
    rhs::V,
) where {C,V<:AbstractVariableRef}
    return add_to_expression!(copy(lhs), -one(C), rhs)
end
# Don't fall back on AbstractVariableRef*GenericAffExpr to preserve lhs/rhs
# consistency (appears in printing).
function Base.:*(
    lhs::GenericAffExpr{C,V},
    rhs::V,
) where {C,V<:AbstractVariableRef}
    if !iszero(lhs.constant)
        result = GenericQuadExpr{C,V}(lhs.constant * rhs)
    else
        result = zero(GenericQuadExpr{C,V})
    end
    for (coef, var) in linear_terms(lhs)
        add_to_expression!(result, coef, var, rhs)
    end
    return result
end
function Base.:/(lhs::GenericAffExpr, rhs::AbstractVariableRef)
    return error("Cannot divide affine expression by a variable")
end
# AffExpr--AffExpr

_copy_convert_coef(::Type{C}, aff::GenericAffExpr{C}) where {C} = copy(aff)
function _copy_convert_coef(::Type{T}, aff::GenericAffExpr{C,V}) where {T,C,V}
    return convert(GenericAffExpr{T,V}, aff)
end

_copy_convert_coef(::Type{C}, quad::GenericQuadExpr{C}) where {C} = copy(quad)
function _copy_convert_coef(::Type{T}, quad::GenericQuadExpr{C,V}) where {T,C,V}
    return convert(GenericQuadExpr{T,V}, quad)
end

function Base.:+(
    lhs::GenericAffExpr{S,V},
    rhs::GenericAffExpr{T,V},
) where {S,T,V<:_JuMPTypes}
    if length(linear_terms(lhs)) > 50 || length(linear_terms(rhs)) > 50
        if length(linear_terms(lhs)) > 1
            operator_warn(owner_model(first(linear_terms(lhs))[2]))
        end
    end
    return add_to_expression!(
        _copy_convert_coef(_MA.promote_operation(+, S, T), lhs),
        rhs,
    )
end

function Base.:-(
    lhs::GenericAffExpr{S,V},
    rhs::GenericAffExpr{T,V},
) where {S,T,V<:_JuMPTypes}
    result = _copy_convert_coef(_MA.promote_operation(-, S, T), lhs)
    result.constant -= rhs.constant
    sizehint!(result, length(linear_terms(lhs)) + length(linear_terms(rhs)))
    for (coef, var) in linear_terms(rhs)
        add_to_expression!(result, -coef, var)
    end
    return result
end

function Base.:*(
    lhs::GenericAffExpr{S,V},
    rhs::GenericAffExpr{T,V},
) where {S,T,V<:_JuMPTypes}
    result = zero(GenericQuadExpr{_MA.promote_sum_mul(S, T),V})
    add_to_expression!(result, lhs, rhs)
    return result
end
# GenericAffExpr--GenericQuadExpr
function Base.:+(a::GenericAffExpr, q::GenericQuadExpr)
    return GenericQuadExpr(a + q.aff, copy(q.terms))
end
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
function Base.:+(q::GenericQuadExpr, v::AbstractVariableRef)
    return GenericQuadExpr(q.aff + v, copy(q.terms))
end
function Base.:-(q::GenericQuadExpr, v::AbstractVariableRef)
    return GenericQuadExpr(q.aff - v, copy(q.terms))
end
function Base.:*(q::GenericQuadExpr, v::AbstractVariableRef)
    return error("Cannot multiply a quadratic expression by a variable")
end
function Base.:/(q::GenericQuadExpr, v::AbstractVariableRef)
    return error("Cannot divide a quadratic expression by a variable")
end
# GenericQuadExpr--GenericAffExpr
function Base.:+(q::GenericQuadExpr, a::GenericAffExpr)
    return GenericQuadExpr(q.aff + a, copy(q.terms))
end
function Base.:-(q::GenericQuadExpr, a::GenericAffExpr)
    return GenericQuadExpr(q.aff - a, copy(q.terms))
end
function Base.:*(q::GenericQuadExpr, a::GenericAffExpr)
    return error("Cannot multiply a quadratic expression by an aff. expression")
end
function Base.:/(q::GenericQuadExpr, a::GenericAffExpr)
    return error("Cannot divide a quadratic expression by an aff. expression")
end
# GenericQuadExpr--GenericQuadExpr
function Base.:+(q1::GenericQuadExpr{S}, q2::GenericQuadExpr{T}) where {S,T}
    result = _copy_convert_coef(_MA.promote_operation(+, S, T), q1)
    for (coef, var1, var2) in quad_terms(q2)
        add_to_expression!(result, coef, var1, var2)
    end
    for (coef, var) in linear_terms(q2)
        add_to_expression!(result, coef, var)
    end
    result.aff.constant += q2.aff.constant
    return result
end
function Base.:-(q1::GenericQuadExpr{S}, q2::GenericQuadExpr{T}) where {S,T}
    result = _copy_convert_coef(_MA.promote_operation(-, S, T), q1)
    for (coef, var1, var2) in quad_terms(q2)
        add_to_expression!(result, -coef, var1, var2)
    end
    for (coef, var) in linear_terms(q2)
        add_to_expression!(result, -coef, var)
    end
    result.aff.constant -= q2.aff.constant
    return result
end

function Base.:(==)(lhs::GenericAffExpr, rhs::GenericAffExpr)
    return (lhs.terms == rhs.terms) && (lhs.constant == rhs.constant)
end
function Base.:(==)(lhs::GenericQuadExpr, rhs::GenericQuadExpr)
    return (lhs.terms == rhs.terms) && (lhs.aff == rhs.aff)
end

# Base Julia's generic fallback vecdot, aka dot, requires that dot, aka LinearAlgebra.dot, be defined
# for scalars, so instead of defining them one-by-one, we will
# fallback to the multiplication operator
LinearAlgebra.dot(lhs::_JuMPTypes, rhs::_JuMPTypes) = lhs * rhs
LinearAlgebra.dot(lhs::_JuMPTypes, rhs::_Constant) = lhs * rhs
LinearAlgebra.dot(lhs::_Constant, rhs::_JuMPTypes) = lhs * rhs

function Base.promote_rule(V::Type{<:AbstractVariableRef}, R::Type{<:Number})
    return GenericAffExpr{_float_type(R),V}
end
function Base.promote_rule(
    V::Type{<:AbstractVariableRef},
    ::Type{<:GenericAffExpr{T}},
) where {T}
    return GenericAffExpr{T,V}
end
function Base.promote_rule(
    V::Type{<:AbstractVariableRef},
    ::Type{<:GenericQuadExpr{T}},
) where {T}
    return GenericQuadExpr{T,V}
end
function Base.promote_rule(
    ::Type{GenericAffExpr{S,V}},
    R::Type{<:Number},
) where {S,V}
    return GenericAffExpr{promote_type(S, R),V}
end
function Base.promote_rule(
    ::Type{<:GenericAffExpr{S,V}},
    ::Type{<:GenericAffExpr{T,V}},
) where {S,T,V}
    return GenericAffExpr{promote_type(S, T),V}
end
function Base.promote_rule(
    ::Type{<:GenericAffExpr{S,V}},
    ::Type{<:GenericQuadExpr{T,V}},
) where {S,T,V}
    return GenericQuadExpr{promote_type(S, T),V}
end
function Base.promote_rule(
    ::Type{GenericQuadExpr{S,V}},
    R::Type{<:Number},
) where {S,V}
    return GenericQuadExpr{promote_type(S, R),V}
end

Base.transpose(x::AbstractJuMPScalar) = x

Base.conj(x::VariableRef) = x
# Can remove the following code once == overloading is removed

function LinearAlgebra.issymmetric(x::Matrix{T}) where {T<:_JuMPTypes}
    (n = size(x, 1)) == size(x, 2) || return false
    for i in 1:n, j in (i+1):n
        isequal(x[i, j], x[j, i]) || return false
    end
    return true
end

###############################################################################
# nonlinear function fallbacks for JuMP built-in types
###############################################################################

const _OP_HINT = "Are you trying to build a nonlinear problem? Make sure you use @NLconstraint/@NLobjective."
for f in MOI.Nonlinear.DEFAULT_UNIVARIATE_OPERATORS
    if f in (:+, :-, :abs) || !isdefined(Base, f)
        continue
    end
    for T in (:AbstractVariableRef, :GenericAffExpr, :GenericQuadExpr)
        if f == :abs2 && (T == :AbstractVariableRef || T == :GenericAffExpr)
            continue
        end
        error_string = "$f is not defined for type $T. $_OP_HINT"
        @eval Base.$(f)(::$T) = error($error_string)
    end
end

function Base.:*(
    ::T,
    ::S,
) where {
    T<:GenericQuadExpr,
    S<:Union{AbstractVariableRef,GenericAffExpr,GenericQuadExpr},
}
    return error("*(::$T,::$S) is not defined. $_OP_HINT")
end
function Base.:*(lhs::GenericQuadExpr, rhs::GenericQuadExpr)
    return error(
        "*(::GenericQuadExpr,::GenericQuadExpr) is not defined. $_OP_HINT",
    )
end
function Base.:*(
    ::S,
    ::T,
) where {
    T<:GenericQuadExpr,
    S<:Union{AbstractVariableRef,GenericAffExpr,GenericQuadExpr},
}
    return error("*(::$S,::$T) is not defined. $_OP_HINT")
end
function Base.:/(
    ::S,
    ::T,
) where {
    S<:Union{_Constant,AbstractVariableRef,GenericAffExpr,GenericQuadExpr},
    T<:Union{AbstractVariableRef,GenericAffExpr,GenericQuadExpr},
}
    return error("/(::$S,::$T) is not defined. $_OP_HINT")
end

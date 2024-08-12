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

_complex_convert_type(::Type{T}, ::Type{<:Real}) where {T} = T
function _complex_convert_type(
    ::Type{T},
    ::Type{<:LinearAlgebra.UniformScaling{S}},
) where {T,S}
    return _complex_convert_type(T, S)
end
_complex_convert_type(::Type{T}, ::Type{<:Complex}) where {T} = Complex{T}

_complex_convert(::Type{T}, x::Real) where {T} = convert(T, x)
_complex_convert(::Type{T}, x::Complex) where {T} = convert(Complex{T}, x)
function _complex_convert(::Type{T}, J::LinearAlgebra.UniformScaling) where {T}
    return _complex_convert(T, J.λ)
end

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
    constant = _complex_convert(value_type(typeof(rhs)), lhs)
    return _build_aff_expr(constant, one(constant), rhs)
end
function Base.:-(lhs::_Constant, rhs::AbstractVariableRef)
    constant = _complex_convert(value_type(typeof(rhs)), lhs)
    return _build_aff_expr(constant, -one(constant), rhs)
end
function Base.:*(lhs::_Constant, rhs::AbstractVariableRef)
    coef = _complex_convert(value_type(typeof(rhs)), lhs)
    if iszero(coef)
        return zero(GenericAffExpr{typeof(coef),typeof(rhs)})
    else
        return _build_aff_expr(zero(coef), coef, rhs)
    end
end
# _Constant--_GenericAffOrQuadExpr
function Base.:+(lhs::_Constant, rhs::_GenericAffOrQuadExpr)
    # If `lhs` is complex and `rhs` has real coefficients then the conversion is needed
    T = _MA.promote_operation(
        +,
        _complex_convert_type(value_type(variable_ref_type(rhs)), typeof(lhs)),
        typeof(rhs),
    )
    result = _MA.mutable_copy(convert(T, rhs))
    add_to_expression!(result, lhs)
    return result
end
function Base.:-(lhs::_Constant, rhs::_GenericAffOrQuadExpr)
    # If `lhs` is complex and `rhs` has real coefficients then the conversion is needed
    T = _MA.promote_operation(
        +,
        _complex_convert_type(value_type(variable_ref_type(rhs)), typeof(lhs)),
        typeof(rhs),
    )
    result = convert(T, -rhs)
    add_to_expression!(result, lhs)
    return result
end
function Base.:*(lhs::_Constant, rhs::_GenericAffOrQuadExpr)
    T = value_type(variable_ref_type(rhs))
    if iszero(lhs)
        # If `lhs` is complex and `rhs` has real coefficients, `zero(rhs)` would not work
        return zero(
            _MA.promote_operation(
                *,
                _complex_convert_type(T, typeof(lhs)),
                typeof(rhs),
            ),
        )
    else
        return map_coefficients(Base.Fix1(*, _complex_convert(T, lhs)), rhs)
    end
end

# AbstractVariableRef (or, AbstractJuMPScalar)
# TODO: What is the role of AbstractJuMPScalar??
Base.:+(lhs::AbstractJuMPScalar) = lhs
function Base.:-(lhs::AbstractVariableRef)
    T = value_type(typeof(lhs))
    return _build_aff_expr(zero(T), -one(T), lhs)
end
Base.:*(lhs::AbstractJuMPScalar) = lhs # make this more generic so extensions don't have to define unary multiplication for our macros
# AbstractVariableRef--_Constant
Base.:+(lhs::AbstractVariableRef, rhs::_Constant) = (+)(rhs, lhs)
Base.:-(lhs::AbstractVariableRef, rhs::_Constant) = (+)(-rhs, lhs)
Base.:*(lhs::AbstractVariableRef, rhs::_Constant) = (*)(rhs, lhs)
function Base.:/(lhs::AbstractVariableRef, rhs::_Constant)
    T = value_type(typeof(lhs))
    return (*)(inv(convert(T, rhs)), lhs)
end
# AbstractVariableRef--AbstractVariableRef
function Base.:+(lhs::V, rhs::V) where {V<:AbstractVariableRef}
    T = value_type(V)
    return _build_aff_expr(zero(T), one(T), lhs, one(T), rhs)
end
function Base.:-(lhs::V, rhs::V) where {V<:AbstractVariableRef}
    T = value_type(V)
    if lhs == rhs
        return zero(GenericAffExpr{T,V})
    else
        return _build_aff_expr(zero(T), one(T), lhs, -one(T), rhs)
    end
end
function Base.:*(lhs::V, rhs::V) where {V<:AbstractVariableRef}
    T = value_type(V)
    return GenericQuadExpr(
        GenericAffExpr{T,V}(),
        UnorderedPair(lhs, rhs) => one(T),
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

function Base.:^(lhs::V, rhs::Integer) where {V<:AbstractVariableRef}
    if rhs == 0
        return one(value_type(V))
    elseif rhs == 1
        return lhs
    elseif rhs == 2
        return lhs * lhs
    else
        return GenericNonlinearExpr{V}(:^, Any[lhs, rhs])
    end
end

function Base.:^(lhs::GenericAffExpr{T,V}, rhs::Integer) where {T,V}
    if rhs == 0
        return one(T)
    elseif rhs == 1
        return lhs
    elseif rhs == 2
        return lhs * lhs
    else
        return GenericNonlinearExpr{V}(:^, Any[lhs, rhs])
    end
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
# AffExpr--AffExpr

_copy_convert_coef(::Type{C}, aff::GenericAffExpr{C}) where {C} = copy(aff)
function _copy_convert_coef(::Type{T}, aff::GenericAffExpr{C,V}) where {T,C,V}
    return convert(GenericAffExpr{T,V}, aff)
end

_copy_convert_coef(::Type{C}, quad::GenericQuadExpr{C}) where {C} = copy(quad)
function _copy_convert_coef(::Type{T}, quad::GenericQuadExpr{C,V}) where {T,C,V}
    return convert(GenericQuadExpr{T,V}, quad)
end

"""
    operator_warn(model::AbstractModel)
    operator_warn(model::GenericModel)

This function is called on the model whenever two affine expressions are added
together without using `destructive_add!`, and at least one of the two
expressions has more than 50 terms.

For the case of `Model`, if this function is called more than 20,000 times then
a warning is generated once.

This method should only be implemented by developers creating JuMP extensions.
It should never be called by users of JuMP.
"""
operator_warn(::AbstractModel) = nothing

function operator_warn(model::GenericModel)
    model.operator_counter += 1
    if model.operator_counter > 20000
        @warn(
            "The addition operator has been used on JuMP expressions a large " *
            "number of times. This warning is safe to ignore but may " *
            "indicate that model generation is slower than necessary. For " *
            "performance reasons, you should not add expressions in a loop. " *
            "Instead of x += y, use add_to_expression!(x,y) to modify x in " *
            "place. If y is a single variable, you may also use " *
            "add_to_expression!(x, coef, y) for x += coef*y.",
            maxlog = 1
        )
    end
    return
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
function Base.:-(a::GenericAffExpr{S}, q::GenericQuadExpr{T}) where {S,T}
    result = -_copy_convert_coef(_MA.promote_operation(-, S, T), q)
    add_to_expression!(result.aff, a)
    return result
end

# GenericQuadExpr
Base.:+(lhs::GenericQuadExpr) = lhs
Base.:-(lhs::GenericQuadExpr) = map_coefficients(-, lhs)
# GenericQuadExpr--_Constant
# We don't do `+rhs` as `LinearAlgebra.UniformScaling` does not support unary `+`
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
# GenericQuadExpr--GenericAffExpr
function Base.:+(q::GenericQuadExpr, a::GenericAffExpr)
    return GenericQuadExpr(q.aff + a, copy(q.terms))
end
function Base.:-(q::GenericQuadExpr, a::GenericAffExpr)
    return GenericQuadExpr(q.aff - a, copy(q.terms))
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
LinearAlgebra.dot(lhs::_JuMPTypes, rhs::_JuMPTypes) = conj(lhs) * rhs
LinearAlgebra.dot(lhs::_JuMPTypes, rhs::_Constant) = conj(lhs) * rhs
LinearAlgebra.dot(lhs::_Constant, rhs::_JuMPTypes) = conj(lhs) * rhs

function Base.promote_rule(V::Type{<:AbstractVariableRef}, R::Type{<:Number})
    return GenericAffExpr{_complex_convert_type(value_type(V), R),V}
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

Base.conj(x::GenericVariableRef) = x
# Can remove the following code once == overloading is removed

function LinearAlgebra.issymmetric(x::Matrix{T}) where {T<:_JuMPTypes}
    (n = size(x, 1)) == size(x, 2) || return false
    for i in 1:n, j in (i+1):n
        isequal(x[i, j], x[j, i]) || return false
    end
    return true
end

function _throw_operator_error(
    ::Union{typeof(+),typeof(_MA.add_mul)},
    x::AbstractArray,
)
    msg =
        "Addition between an array and a JuMP scalar is not supported: " *
        "instead of `x + y`, do `x .+ y` for element-wise addition."
    if ndims(x) == 2 && size(x, 1) == size(x, 2)
        msg *=
            " If you are modifying the diagonal entries of a square matrix, " *
            "do `x + y * LinearAlgebra.I(n)`, where `n` is the side length."
    end
    return error(msg)
end

function _throw_operator_error(
    ::Union{typeof(-),typeof(_MA.sub_mul)},
    x::AbstractArray,
)
    msg =
        "Subtraction between an array and a JuMP scalar is not supported: " *
        "instead of `x - y`, do `x .- y` for element-wise subtraction."
    if ndims(x) == 2 && size(x, 1) == size(x, 2)
        msg *=
            " If you are modifying the diagonal entries of a square matrix, " *
            "do `x - y * LinearAlgebra.I(n)`, where `n` is the side length."
    end
    return error(msg)
end

Base.:+(::AbstractJuMPScalar, x::AbstractArray) = _throw_operator_error(+, x)
Base.:+(x::AbstractArray, ::AbstractJuMPScalar) = _throw_operator_error(+, x)
Base.:-(::AbstractJuMPScalar, x::AbstractArray) = _throw_operator_error(-, x)
Base.:-(x::AbstractArray, ::AbstractJuMPScalar) = _throw_operator_error(-, x)

function _MA.operate!!(
    op::Union{typeof(_MA.add_mul),typeof(_MA.sub_mul)},
    x::AbstractArray,
    ::AbstractJuMPScalar,
)
    return _throw_operator_error(op, x)
end

function _MA.operate!!(
    op::Union{typeof(_MA.add_mul),typeof(_MA.sub_mul)},
    ::AbstractJuMPScalar,
    x::AbstractArray,
)
    return _throw_operator_error(op, x)
end

_mult_upper(α, A) = parent(α * LinearAlgebra.UpperTriangular(parent(A)))
_mult_lower(α, A) = parent(α * LinearAlgebra.LowerTriangular(parent(A)))

function Base.:*(
    x::Union{
        GenericVariableRef{<:Real},
        GenericAffExpr{<:Real},
        GenericQuadExpr{<:Real},
    },
    A::LinearAlgebra.Hermitian,
)
    c = LinearAlgebra.sym_uplo(A.uplo)
    B = c == :U ? _mult_upper(x, A) : _mult_lower(x, A)
    # Intermediate conversion to `Matrix` is needed to work around
    # https://github.com/JuliaLang/julia/issues/52895
    return LinearAlgebra.Hermitian(Matrix(LinearAlgebra.Hermitian(B, c)), c)
end

function Base.:*(
    A::LinearAlgebra.Hermitian,
    x::Union{
        GenericVariableRef{<:Real},
        GenericAffExpr{<:Real},
        GenericQuadExpr{<:Real},
    },
)
    return x * A
end

function Base.complex(
    r::Union{
        GenericVariableRef{<:Real},
        GenericAffExpr{<:Real},
        GenericQuadExpr{<:Real},
    },
    i::Real,
)
    return r + im * i
end

function Base.complex(
    r::Real,
    i::Union{
        GenericVariableRef{<:Real},
        GenericAffExpr{<:Real},
        GenericQuadExpr{<:Real},
    },
)
    return r + im * i
end

function Base.complex(
    r::Union{
        GenericVariableRef{<:Real},
        GenericAffExpr{<:Real},
        GenericQuadExpr{<:Real},
    },
    i::Union{
        GenericVariableRef{<:Real},
        GenericAffExpr{<:Real},
        GenericQuadExpr{<:Real},
    },
)
    return r + im * i
end

# These methods exist in LinearAlgebra for subtypes of Real. Without them, we
# return a `Matrix` which looses the Hermitian information.
function Base.:+(
    A::LinearAlgebra.Symmetric{V,Matrix{V}},
    B::LinearAlgebra.Hermitian,
) where {
    V<:Union{
        GenericVariableRef{<:Real},
        GenericAffExpr{<:Real},
        GenericQuadExpr{<:Real},
    },
}
    return LinearAlgebra.Hermitian(A) + B
end

function Base.:+(
    A::LinearAlgebra.Hermitian,
    B::LinearAlgebra.Symmetric{V,Matrix{V}},
) where {
    V<:Union{
        GenericVariableRef{<:Real},
        GenericAffExpr{<:Real},
        GenericQuadExpr{<:Real},
    },
}
    return A + LinearAlgebra.Hermitian(B)
end

function Base.:-(
    A::LinearAlgebra.Symmetric{V,Matrix{V}},
    B::LinearAlgebra.Hermitian,
) where {
    V<:Union{
        GenericVariableRef{<:Real},
        GenericAffExpr{<:Real},
        GenericQuadExpr{<:Real},
    },
}
    return LinearAlgebra.Hermitian(A) - B
end

function Base.:-(
    A::LinearAlgebra.Hermitian,
    B::LinearAlgebra.Symmetric{V,Matrix{V}},
) where {
    V<:Union{
        GenericVariableRef{<:Real},
        GenericAffExpr{<:Real},
        GenericQuadExpr{<:Real},
    },
}
    return A - LinearAlgebra.Hermitian(B)
end

#  Copyright 2017, Iain Dunning, Joey Huchette, Miles Lubin, and contributors
#  This Source Code Form is subject to the terms of the Mozilla Public
#  License, v. 2.0. If a copy of the MPL was not distributed with this
#  file, You can obtain one at http://mozilla.org/MPL/2.0/.
#############################################################################
# JuMP
# An algebraic modeling language for Julia
# See http://github.com/JuliaOpt/JuMP.jl
#############################################################################

const JuMPTypes = Union{AbstractJuMPScalar, NonlinearExpression}

# Overloads
#
# Different objects that must all interact:
# 1. Number
# 2. AbstractVariableRef
# 4. GenericAffExpr
# 5. GenericQuadExpr

# Number
# Number--Number obviously already taken care of!
# Number--VariableRef
Base.:+(lhs::Number, rhs::AbstractVariableRef) = GenericAffExpr(convert(Float64, lhs), rhs => 1.0)
Base.:-(lhs::Number, rhs::AbstractVariableRef) = GenericAffExpr(convert(Float64, lhs), rhs => -1.0)
Base.:*(lhs::Number, rhs::AbstractVariableRef) = GenericAffExpr(0.0, rhs => convert(Float64,lhs))
# Number--GenericAffExpr
function Base.:+(lhs::Number, rhs::GenericAffExpr)
    result = copy(rhs)
    result.constant += lhs
    return result
end
function Base.:-(lhs::Number, rhs::GenericAffExpr)
    result = -rhs
    result.constant += lhs
    return result
end
Base.:*(lhs::Number, rhs::GenericAffExpr) = map_coefficients(c -> lhs * c, rhs)
# Number--QuadExpr
Base.:+(lhs::Number, rhs::GenericQuadExpr) = GenericQuadExpr(lhs+rhs.aff, copy(rhs.terms))
function Base.:-(lhs::Number, rhs::GenericQuadExpr)
    result = -rhs
    result.aff.constant += lhs
    return result
end
Base.:*(lhs::Number, rhs::GenericQuadExpr) = map_coefficients(c -> lhs * c, rhs)

# AbstractVariableRef (or, AbstractJuMPScalar)
# TODO: What is the role of AbstractJuMPScalar??
Base.:+(lhs::AbstractJuMPScalar) = lhs
Base.:-(lhs::AbstractVariableRef) = GenericAffExpr(0.0, lhs => -1.0)
Base.:*(lhs::AbstractJuMPScalar) = lhs # make this more generic so extensions don't have to define unary multiplication for our macros
# AbstractVariableRef--Number
Base.:+(lhs::AbstractVariableRef, rhs::Number) = (+)( rhs,lhs)
Base.:-(lhs::AbstractVariableRef, rhs::Number) = (+)(-rhs,lhs)
Base.:*(lhs::AbstractVariableRef, rhs::Number) = (*)(rhs,lhs)
Base.:/(lhs::AbstractVariableRef, rhs::Number) = (*)(1.0/rhs,lhs)
# AbstractVariableRef--AbstractVariableRef
Base.:+(lhs::V, rhs::V) where V <: AbstractVariableRef = GenericAffExpr(0.0, lhs => 1.0, rhs => 1.0)
Base.:-(lhs::V, rhs::V) where V <: AbstractVariableRef = GenericAffExpr(0.0, lhs => 1.0, rhs => -1.0)
function Base.:*(lhs::V, rhs::V) where V <: AbstractVariableRef
    GenericQuadExpr(GenericAffExpr{Float64, V}(), UnorderedPair(lhs, rhs) => 1.0)
end
# AbstractVariableRef--GenericAffExpr
function Base.:+(lhs::V, rhs::GenericAffExpr{C,V}) where {C, V <: AbstractVariableRef}
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

function Base.:*(lhs::V, rhs::GenericAffExpr{C,V}) where {C, V <: AbstractVariableRef}
    if !iszero(rhs.constant)
        result = GenericQuadExpr{C,V}(lhs*rhs.constant)
    else
        result = zero(GenericQuadExpr{C,V})
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
# GenericAffExpr--Number
Base.:+(lhs::GenericAffExpr, rhs::Number) = (+)(+rhs,lhs)
Base.:-(lhs::GenericAffExpr, rhs::Number) = (+)(-rhs,lhs)
Base.:*(lhs::GenericAffExpr, rhs::Number) = (*)(rhs,lhs)
Base.:/(lhs::GenericAffExpr, rhs::Number) = map_coefficients(c -> c/rhs, lhs)
function Base.:^(lhs::Union{AbstractVariableRef,GenericAffExpr}, rhs::Integer)
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
Base.:^(lhs::Union{AbstractVariableRef,GenericAffExpr}, rhs::Number) = error("Only exponents of 0, 1, or 2 are currently supported. Are you trying to build a nonlinear problem? Make sure you use @NLconstraint/@NLobjective.")
# GenericAffExpr--AbstractVariableRef
function Base.:+(lhs::GenericAffExpr{C,V}, rhs::V) where {C, V <: AbstractVariableRef}
    return add_to_expression!(copy(lhs), one(C), rhs)
end
function Base.:-(lhs::GenericAffExpr{C,V}, rhs::V) where {C, V <: AbstractVariableRef}
    return add_to_expression!(copy(lhs), -one(C), rhs)
end
# Don't fall back on AbstractVariableRef*GenericAffExpr to preserve lhs/rhs
# consistency (appears in printing).
function Base.:*(lhs::GenericAffExpr{C,V}, rhs::V) where {C, V <: AbstractVariableRef}
    if !iszero(lhs.constant)
        result = GenericQuadExpr{C,V}(lhs.constant*rhs)
    else
        result = zero(GenericQuadExpr{C,V})
    end
    for (coef, var) in linear_terms(lhs)
        add_to_expression!(result, coef, var, rhs)
    end
    return result
end
Base.:/(lhs::GenericAffExpr, rhs::AbstractVariableRef) = error("Cannot divide affine expression by a variable")
# AffExpr--AffExpr
function Base.:+(lhs::GenericAffExpr{C,V}, rhs::GenericAffExpr{C,V}) where {C,V<:JuMPTypes}
    if length(linear_terms(lhs)) > 50 || length(linear_terms(rhs)) > 50
        if length(linear_terms(lhs)) > 1
            operator_warn(owner_model(first(linear_terms(lhs))[2]))
        end
    end
    result_terms = copy(lhs.terms)
    # merge() returns a Dict(), so we need to call merge!() instead.
    # Note: merge!() doesn't appear to call sizehint!(). Is this important?
    merge!(+, result_terms, rhs.terms)
    return GenericAffExpr(lhs.constant + rhs.constant, result_terms)
end

function Base.:-(lhs::GenericAffExpr{C,V}, rhs::GenericAffExpr{C,V}) where {C,V<:JuMPTypes}
    result = copy(lhs)
    result.constant -= rhs.constant
    sizehint!(result, length(linear_terms(lhs)) + length(linear_terms(rhs)))
    for (coef, var) in linear_terms(rhs)
        add_to_expression!(result, -coef, var)
    end
    return result
end

function Base.:*(lhs::GenericAffExpr{C,V}, rhs::GenericAffExpr{C,V}) where {C,V<:JuMPTypes}
    result = zero(GenericQuadExpr{C,V})

    lhs_length = length(linear_terms(lhs))
    rhs_length = length(linear_terms(rhs))

    # Quadratic terms
    for (lhscoef, lhsvar) in linear_terms(lhs)
        for (rhscoef, rhsvar) in linear_terms(rhs)
            add_to_expression!(result, lhscoef*rhscoef, lhsvar, rhsvar)
        end
    end

    # Try to preallocate space for aff
    if !iszero(lhs.constant) && !iszero(rhs.constant)
        sizehint!(result.aff, lhs_length + rhs_length)
    elseif !iszero(lhs.constant)
        sizehint!(result.aff, rhs_length)
    elseif !iszero(rhs.constant)
        sizehint!(result.aff, lhs_length)
    end

    # [LHS constant] * [RHS linear terms]
    if !iszero(lhs.constant)
        c = lhs.constant
        for (rhscoef, rhsvar) in linear_terms(rhs)
            add_to_expression!(result.aff, c*rhscoef, rhsvar)
        end
    end

    # [RHS constant] * [LHS linear terms]
    if !iszero(rhs.constant)
        c = rhs.constant
        for (lhscoef, lhsvar) in linear_terms(lhs)
            add_to_expression!(result.aff, c*lhscoef, lhsvar)
        end
    end

    result.aff.constant = lhs.constant * rhs.constant

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
# GenericQuadExpr--Number
Base.:+(lhs::GenericQuadExpr, rhs::Number) = (+)(+rhs,lhs)
Base.:-(lhs::GenericQuadExpr, rhs::Number) = (+)(-rhs,lhs)
Base.:*(lhs::GenericQuadExpr, rhs::Number) = (*)(rhs,lhs)
Base.:/(lhs::GenericQuadExpr, rhs::Number) = (*)(inv(rhs),lhs)
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
    for (coef, var1, var2) in quadterms(q2)
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
    for (coef, var1, var2) in quadterms(q2)
        add_to_expression!(result, -coef, var1, var2)
    end
    for (coef, var) in linear_terms(q2)
        add_to_expression!(result, -coef, var)
    end
    result.aff.constant -= q2.aff.constant
    return result
end

Base.:(==)(lhs::GenericAffExpr,rhs::GenericAffExpr) = (lhs.terms == rhs.terms) && (lhs.constant == rhs.constant)
Base.:(==)(lhs::GenericQuadExpr,rhs::GenericQuadExpr) = (lhs.terms == rhs.terms) && (lhs.aff == rhs.aff)

#############################################################################
# Helpers to initialize memory for GenericAffExpr/GenericQuadExpr
#############################################################################

_sizehint_expr!(a::GenericAffExpr, n::Int) = sizehint!(a, n)

# TODO: Why do we allocate the same size for the quadratic and affine parts?
function _sizehint_expr!(q::GenericQuadExpr, n::Int)
        sizehint!(q.terms,  length(q.terms) + n)
        _sizehint_expr!(q.aff, n)
        nothing
end
_sizehint_expr!(q, n) = nothing

#############################################################################
# High-level operators
# Currently supported
#  - sum
#  - dot
#############################################################################

# TODO: specialize sum for Dict and JuMPArray of JuMP objects?
Base.sum(vars::Array{<:AbstractVariableRef}) = GenericAffExpr(0.0, [v => 1.0 for v in vars])
function Base.sum(array::AbstractArray{<:AbstractVariableRef})
    result_expression = zero(eltype(array))
    for variable in array
        add_to_expression!(result_expression, variable)
    end
    return result_expression
end
# TODO: Specialize for iterables.
function Base.sum(affs::AbstractArray{T}) where T<:GenericAffExpr
    new_aff = zero(T)
    for aff in affs
        add_to_expression!(new_aff, aff)
    end
    return new_aff
end

# Base Julia's generic fallback vecdot, aka Compat.dot, requires that dot, aka Compat.LinearAlgebra.dot, be defined
# for scalars, so instead of defining them one-by-one, we will
# fallback to the multiplication operator
Compat.LinearAlgebra.dot(lhs::JuMPTypes, rhs::JuMPTypes) = lhs*rhs
Compat.LinearAlgebra.dot(lhs::JuMPTypes, rhs::Number)    = lhs*rhs
Compat.LinearAlgebra.dot(lhs::Number,    rhs::JuMPTypes) = lhs*rhs

Compat.dot(lhs::AbstractVector{T},rhs::AbstractVector{S}) where {T<:JuMPTypes,S<:JuMPTypes} = _dot(lhs,rhs)
Compat.dot(lhs::AbstractVector{T},rhs::AbstractVector{S}) where {T<:JuMPTypes,S} = _dot(lhs,rhs)
Compat.dot(lhs::AbstractVector{T},rhs::AbstractVector{S}) where {T,S<:JuMPTypes} = _dot(lhs,rhs)

Compat.dot(lhs::AbstractArray{T,N},rhs::AbstractArray{S,N}) where {T<:JuMPTypes,S,N} = _dot(lhs,rhs)
Compat.dot(lhs::AbstractArray{T,N},rhs::AbstractArray{S,N}) where {T<:JuMPTypes,S<:JuMPTypes,N} = _dot(lhs,rhs)
Compat.dot(lhs::AbstractArray{T,N},rhs::AbstractArray{S,N}) where {T,S<:JuMPTypes,N} = _dot(lhs,rhs)

function _dot(lhs::AbstractArray{T}, rhs::AbstractArray{S}) where {T,S}
    size(lhs) == size(rhs) || error("Incompatible dimensions")
    ret = zero(one(T)*one(S))
    for (x,y) in zip(lhs,rhs)
        ret = destructive_add!(ret, x, y)
    end
    ret
end

###############################################################################
# A bunch of operator junk to make matrix multiplication and friends act
# reasonably sane with JuMP types

Base.promote_rule(V::Type{<:AbstractVariableRef},R::Type{<:Real}                 )                 = GenericAffExpr{Float64, V}
Base.promote_rule(V::Type{<:AbstractVariableRef}, ::Type{<:GenericAffExpr{T}}    ) where T         = GenericAffExpr{T, V}
Base.promote_rule(V::Type{<:AbstractVariableRef}, ::Type{<:GenericQuadExpr{T}}   ) where T         = GenericQuadExpr{T, V}
Base.promote_rule(::Type{GenericAffExpr{S, V}},  R::Type{<:Real}                 ) where {S, V}    = GenericAffExpr{promote_type(S, R), V}
Base.promote_rule(::Type{<:GenericAffExpr{S, V}}, ::Type{<:GenericQuadExpr{T, V}}) where {S, T, V} = GenericQuadExpr{promote_type(S, T), V}
Base.promote_rule(::Type{GenericQuadExpr{S, V}}, R::Type{<:Real}                 ) where {S, V}    = GenericQuadExpr{promote_type(S, R), V}

Base.transpose(x::AbstractJuMPScalar) = x

# Can remove the following code once == overloading is removed

function Compat.LinearAlgebra.issymmetric(x::Matrix{T}) where T<:JuMPTypes
    (n = size(x,1)) == size(x,2) || return false
    for i in 1:n, j in (i+1):n
        isequal(x[i,j], x[j,i]) || return false
    end
    true
end

# Special-case because the the base version wants to do fill!(::Array{AbstractVariableRef}, zero(GenericAffExpr{Float64,eltype(x)}))
one_indexed(A) = all(x -> isa(x, Base.OneTo), Compat.axes(A))
function Compat.LinearAlgebra.diagm(x::AbstractVector{<:AbstractVariableRef})
    @assert one_indexed(x) # Base.diagm doesn't work for non-one-indexed arrays in general.
    diagm(0=>copyto!(similar(x, GenericAffExpr{Float64,eltype(x)}), x))
end

###############################################################################
# Matrix/Vector Arithmetic with JuMP eltypes
###############################################################################

###############################################################################
# convenience/utility definitions

const GenericAffOrQuadExpr = Union{GenericAffExpr,GenericQuadExpr}

densify_with_jump_eltype(x::AbstractMatrix) = convert(Matrix, x)

function densify_with_jump_eltype(x::SparseMatrixCSC{V}) where V<:AbstractVariableRef
    return convert(Matrix{GenericAffExpr{Float64,V}}, x)
end

# see https://github.com/JuliaLang/julia/pull/18218
_A_mul_B_eltype(::Type{R}, ::Type{S}) where {R,S} = typeof(one(R) * one(S) + one(R) * one(S))

_A_mul_B_ret_dims(A::AbstractMatrix, B::AbstractVector) = (size(A, 1),)
_A_mul_B_ret_dims(A::AbstractMatrix, B::AbstractMatrix) = (size(A, 1), size(B, 2))

# don't do size checks here; defer that to `_A_mul_B(A, B)`
function _A_mul_B_ret(A, B, dims...)
    T = _A_mul_B_eltype(eltype(A), eltype(B))
    ret = Array{T}(undef, _A_mul_B_ret_dims(A, B))
    return fillwithzeros!(ret, T)
end

if VERSION < v"0.7-"
    _At_mul_B_ret_dims(A::AbstractMatrix, B::AbstractVector) = (size(A, 2),)
    _At_mul_B_ret_dims(A::AbstractMatrix, B::AbstractMatrix) = (size(A, 2), size(B, 2))
    function _At_mul_B_ret(A, B, dims...)
        T = _A_mul_B_eltype(eltype(A), eltype(B))
        ret = Array{T}(undef, _At_mul_B_ret_dims(A, B))
        return fillwithzeros!(ret, T)
    end
end

function fillwithzeros!(A, ::Type{T}) where {T}
    for I in eachindex(A)
        A[I] = zero(T)
    end
    return A
end

###############################################################################
# `_A_mul_B!(ret, A, B)` adds the result of `A*B` into the buffer `ret`. We roll our own
# matmul here (instead of using Julia's generic fallbacks) because doing so allows us to
# accumulate the expressions for the inner loops in-place. Additionally, Julia's generic
# fallbacks can be finnicky when your array elements aren't `<:Number`.
#
# No bounds/size checks are performed; it is expected that the caller has done this, has
# ensured that the eltype of `ret` is appropriate, and has zeroed the elements of `ret` (if
# desired).

function _A_mul_B!(ret::AbstractArray{T}, A, B) where {T<:JuMPTypes}
    for i ∈ 1:size(A, 1), j ∈ 1:size(B, 2)
        q = ret[i, j]
        _sizehint_expr!(q, size(A, 2))
        for k ∈ 1:size(A, 2)
            tmp = convert(T, A[i, k] * B[k, j])
            add_to_expression!(q, tmp)
        end
    end
    ret
end

function _A_mul_B!(ret::AbstractArray{<:GenericAffOrQuadExpr}, A::SparseMatrixCSC, B)
    nzv = nonzeros(A)
    rv  = rowvals(A)
    for col ∈ 1:size(A, 2)
        for k ∈ 1:size(ret, 2)
            for j ∈ nzrange(A, col)
                add_to_expression!(ret[rv[j], k], nzv[j] * B[col, k])
            end
        end
    end
    ret
end

function _A_mul_B!(ret::AbstractArray{<:GenericAffOrQuadExpr}, A::AbstractMatrix, B::SparseMatrixCSC)
    rowval = rowvals(B)
    nzval  = nonzeros(B)
    for multivec_row in 1:size(A, 1)
        for col ∈ 1:size(B, 2)
            idxset = nzrange(B, col)
            q = ret[multivec_row, col]
            _sizehint_expr!(q, length(idxset))
            for k ∈ idxset
                add_to_expression!(q, A[multivec_row, rowval[k]] * nzval[k])
            end
        end
    end
    ret
end

# TODO: implement sparse * sparse code as in base/sparse/linalg.jl (spmatmul)
_A_mul_B!(ret::AbstractArray{<:GenericAffOrQuadExpr}, A::SparseMatrixCSC, B::SparseMatrixCSC) = _A_mul_B!(ret, A, densify_with_jump_eltype(B))

function _A_mul_B(A, B)
    size(A, 2) == size(B, 1) || error("Incompatible sizes")
    ret = _A_mul_B_ret(A, B)
    _A_mul_B!(ret, A, B)
    ret
end

###############################################################################
# `_At_mul_B!(ret, A, B)` stores the result of `Aᵀ * B` into the buffer `ret`. We roll our
# own version here (instead of working with Julia's generic fallbacks) for the same reasons
# as above.

function _At_mul_B!(ret::AbstractArray{T}, A, B) where {T<:JuMPTypes}
    for i ∈ 1:size(A, 2), j ∈ 1:size(B, 2)
        q = ret[i, j]
        _sizehint_expr!(q, size(A, 1))
        for k ∈ 1:size(A, 1)
            tmp = convert(T, A[k, i] * B[k, j]) # transpose
            add_to_expression!(q, tmp)
        end
    end
    ret
end

if VERSION < v"0.7-"
    function _At_mul_B!(ret::AbstractArray{<:GenericAffOrQuadExpr}, A::SparseMatrixCSC, B)
        _A_mul_B!(ret, transpose(A), B) # TODO fully implement
    end
    function _At_mul_B!(ret::AbstractArray{<:GenericAffOrQuadExpr}, A::AbstractMatrix, B::SparseMatrixCSC)
        rowval = rowvals(A)
        nzval  = nonzeros(A)
        for multivec_row ∈ 1:size(B, 2) # transpose
            for col ∈ 1:size(A, 2)
                idxset = nzrange(A, col)
                q = ret[multivec_row, col]
                _sizehint_expr!(q, length(idxset))
                for k ∈ idxset
                    add_to_expression!(q, B[rowval[k], multivec_row] * nzval[k]) # transpose
                end
            end
        end
        ret
    end
    function _At_mul_B(A, B)
        size(A, 1) == size(B, 1) || error("Incompatible sizes")
        ret = _At_mul_B_ret(A, B)
        _At_mul_B!(ret, A, B)
        ret
    end
else
    function _At_mul_B!(ret::AbstractArray{<:GenericAffOrQuadExpr}, A::SparseMatrixCSC, B)
        _A_mul_B!(ret, copy(transpose(A)), B) # TODO fully implement
    end
    function _At_mul_B!(ret::AbstractArray{<:GenericAffOrQuadExpr}, A::AbstractMatrix, B::SparseMatrixCSC)
        _A_mul_B!(ret, transpose(A), B)
    end
    # This method of `_At_mul_B!` is adapted from upstream Julia. Note that we
    # confuse transpose with adjoint because they're the same for all JuMP types.
    #=
    > Copyright (c) 2009-2018: Jeff Bezanson, Stefan Karpinski, Viral B. Shah,
    > and other contributors:
    >
    > https://github.com/JuliaLang/julia/contributors
    >
    > Permission is hereby granted, free of charge, to any person obtaining
    > a copy of this software and associated documentation files (the
    > "Software"), to deal in the Software without restriction, including
    > without limitation the rights to use, copy, modify, merge, publish,
    > distribute, sublicense, and/or sell copies of the Software, and to
    > permit persons to whom the Software is furnished to do so, subject to
    > the following conditions:
    >
    > The above copyright notice and this permission notice shall be
    > included in all copies or substantial portions of the Software.
    >
    > THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
    > EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
    > MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
    > NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE
    > LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION
    > OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION
    > WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
    =#
    function _At_mul_B!(ret::StridedVecOrMat{<:GenericAffOrQuadExpr}, A::SparseMatrixCSC, B::StridedVecOrMat)
        A.n == size(ret, 1) || throw(DimensionMismatch())
        A.m == size(B, 1) || throw(DimensionMismatch())
        size(B, 2) == size(ret, 2) || throw(DimensionMismatch())
        nzv = A.nzval
        rv = A.rowval
        # ret is already filled with zeros by _return_arrayt.
        for k = 1:size(ret, 2)
            @inbounds for col = 1:A.n
                tmp = zero(eltype(ret))
                for j = A.colptr[col]:(A.colptr[col + 1] - 1)
                    tmp += adjoint(nzv[j]) * B[rv[j],k]
                end
                ret[col,k] += tmp
            end
        end
        ret
    end
    function _At_mul_B(A, B)
        size(A, 1) == size(B, 1) || error("Incompatible sizes")
        ret = _A_mul_B_ret(transpose(A), B)
        _At_mul_B!(ret, A, B)
        ret
    end
end

# TODO: implement sparse * sparse code as in base/sparse/linalg.jl (spmatmul)
_At_mul_B!(ret::AbstractArray{<:GenericAffOrQuadExpr}, A::SparseMatrixCSC, B::SparseMatrixCSC) = _At_mul_B!(ret, A, densify_with_jump_eltype(B))

###############################################################################
# Interception of Base's matrix/vector arithmetic machinery

# TODO: Intercepting "externally owned" method calls by dispatching on type parameters
# (rather than outermost wrapper type) is generally bad practice, but refactoring this code
# to use a different mechanism would be a lot of work. In the future, this interception code
# would be more easily/robustly replaced by using a tool like
# https://github.com/jrevels/Cassette.jl.

function Base.:*(A::Union{Matrix{<:JuMPTypes},SparseMatrixCSC{<:JuMPTypes}},
                 B::Union{Matrix,Vector,SparseMatrixCSC})
    return _A_mul_B(A, B)
end

function Base.:*(A::Union{Matrix{<:JuMPTypes},SparseMatrixCSC{<:JuMPTypes}},
                 B::Union{Matrix{<:JuMPTypes},Vector{<:JuMPTypes},SparseMatrixCSC{<:JuMPTypes}})
    return _A_mul_B(A, B)
end

function Base.:*(A::Union{Matrix,SparseMatrixCSC},
                 B::Union{Matrix{<:JuMPTypes},Vector{<:JuMPTypes},SparseMatrixCSC{<:JuMPTypes}})
    return _A_mul_B(A, B)
end

if VERSION < v"0.7-"
    # these methods are called when one does A.'*B
    Base.At_mul_B(A::Union{Matrix{T},SparseMatrixCSC{T}}, B::Union{Matrix, Vector, SparseMatrixCSC}) where {T<:JuMPTypes} = _At_mul_B(A, B)
    Base.At_mul_B(A::Union{Matrix{T},SparseMatrixCSC{T}}, B::Union{Matrix{R}, Vector{R}, SparseMatrixCSC{R}}) where {T<:JuMPTypes,R<:JuMPTypes} = _At_mul_B(A, B)
    Base.At_mul_B(A::Union{Matrix,SparseMatrixCSC}, B::Union{Matrix{T}, Vector{T}, SparseMatrixCSC{T}}) where {T<:JuMPTypes} = _At_mul_B(A, B)

    # these methods are called when one does A'*B
    # these methods are the same as above since complex numbers are not implemented in JuMP
    Base.Ac_mul_B(A::Union{Matrix{T},SparseMatrixCSC{T}}, B::Union{Matrix, Vector, SparseMatrixCSC}) where {T<:JuMPTypes} = _At_mul_B(A, B)
    Base.Ac_mul_B(A::Union{Matrix{T},SparseMatrixCSC{T}}, B::Union{Matrix{R}, Vector{R}, SparseMatrixCSC{R}}) where {T<:JuMPTypes,R<:JuMPTypes} = _At_mul_B(A, B)
    Base.Ac_mul_B(A::Union{Matrix,SparseMatrixCSC}, B::Union{Matrix{T}, Vector{T}, SparseMatrixCSC{T}}) where {T<:JuMPTypes} = _At_mul_B(A, B)
else
    # TODO: This is a stopgap solution to get (most) tests passing on Julia 0.7. A lot of
    # cases probably still don't work. (Like A * A where A is a sparse matrix of a JuMP
    # type). This code needs a big refactor.
    Base.:*(A::Adjoint{<:JuMPTypes,<:SparseMatrixCSC}, B::Vector) = _At_mul_B(parent(A), B)
    Base.:*(A::Adjoint{<:Any,<:SparseMatrixCSC}, B::Vector{<:JuMPTypes}) = _At_mul_B(parent(A), B)
    Base.:*(A::Adjoint{<:JuMPTypes,<:SparseMatrixCSC}, B::Vector{<:JuMPTypes}) = _At_mul_B(parent(A), B)

    Base.:*(A::Transpose{<:JuMPTypes,<:SparseMatrixCSC}, B::Vector) = _At_mul_B(parent(A), B)
    Base.:*(A::Transpose{<:Any,<:SparseMatrixCSC}, B::Vector{<:JuMPTypes}) = _At_mul_B(parent(A), B)
    Base.:*(A::Transpose{<:JuMPTypes,<:SparseMatrixCSC}, B::Vector{<:JuMPTypes}) = _At_mul_B(parent(A), B)

    Base.:*(A::Adjoint{<:JuMPTypes,<:SparseMatrixCSC}, B::Matrix) = _At_mul_B(parent(A), B)
    Base.:*(A::Adjoint{<:Any,<:SparseMatrixCSC}, B::Matrix{<:JuMPTypes}) = _At_mul_B(parent(A), B)
    Base.:*(A::Adjoint{<:JuMPTypes,<:SparseMatrixCSC}, B::Matrix{<:JuMPTypes}) = _At_mul_B(parent(A), B)

    Base.:*(A::Transpose{<:JuMPTypes,<:SparseMatrixCSC}, B::Matrix) = _At_mul_B(parent(A), B)
    Base.:*(A::Transpose{<:Any,<:SparseMatrixCSC}, B::Matrix{<:JuMPTypes}) = _At_mul_B(parent(A), B)
    Base.:*(A::Transpose{<:JuMPTypes,<:SparseMatrixCSC}, B::Matrix{<:JuMPTypes}) = _At_mul_B(parent(A), B)
end

Base.:*(A::Number, B::SparseMatrixCSC{T}) where {T<:JuMPTypes} = SparseMatrixCSC(B.m, B.n, copy(B.colptr), copy(B.rowval), A .* B.nzval)
Base.:*(A::JuMPTypes, B::SparseMatrixCSC) = SparseMatrixCSC(B.m, B.n, copy(B.colptr), copy(B.rowval), A .* B.nzval)
Base.:*(A::SparseMatrixCSC{T}, B::Number) where {T<:JuMPTypes} = SparseMatrixCSC(A.m, A.n, copy(A.colptr), copy(A.rowval), A.nzval .* B)
Base.:*(A::SparseMatrixCSC, B::JuMPTypes) = SparseMatrixCSC(A.m, A.n, copy(A.colptr), copy(A.rowval), A.nzval .* B)

Base.:*(x::AbstractArray{T}) where {T<:JuMPTypes} = x

Base.:/(A::SparseMatrixCSC{T}, B::Number) where {T<:JuMPTypes} = SparseMatrixCSC(A.m, A.n, copy(A.colptr), copy(A.rowval), A.nzval ./ B)

Base.:+(x::AbstractArray{T}) where {T<:JuMPTypes} = x

function Base.:-(x::AbstractArray{T}) where T<:JuMPTypes
    ret = similar(x, typeof(-one(T)))
    for I in eachindex(ret)
        ret[I] = -x[I]
    end
    ret
end

# TODO This will interact poorly with other packages that overload broadcast on Julia >v0.7;
# we should be using the new broadcast interface instead.
for (op, opsymbol) in [(+, :+), (-, :-), (*, :*), (/, :/)]
    @eval begin
        Base.broadcast(::typeof($op), x::Number, y::JuMPTypes) = $opsymbol(x, y)
        Base.broadcast(::typeof($op), x::JuMPTypes, y::Number) = $opsymbol(x, y)
    end
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

Base.:*(::T,::S) where {T<:GenericQuadExpr,S<:Union{AbstractVariableRef,GenericAffExpr,GenericQuadExpr}} =
    error( "*(::$T,::$S) is not defined. $op_hint")
Base.:*(lhs::GenericQuadExpr, rhs::GenericQuadExpr) =
    error( "*(::GenericQuadExpr,::GenericQuadExpr) is not defined. $op_hint")
Base.:*(::S,::T) where {T<:GenericQuadExpr,S<:Union{AbstractVariableRef,GenericAffExpr,GenericQuadExpr}} =
    error( "*(::$S,::$T) is not defined. $op_hint")
Base.:/(::S,::T) where {S<:Union{Number,AbstractVariableRef,GenericAffExpr,GenericQuadExpr},T<:Union{AbstractVariableRef,GenericAffExpr,GenericQuadExpr}} =
    error( "/(::$S,::$T) is not defined. $op_hint")

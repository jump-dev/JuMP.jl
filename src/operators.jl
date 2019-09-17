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
const Constant = Union{Number, UniformScaling}
_float(x::Number) = convert(Float64, x)
_float(J::UniformScaling) = _float(J.λ)

# Overloads
#
# Different objects that must all interact:
# 1. Constant
# 2. AbstractVariableRef
# 4. GenericAffExpr
# 5. GenericQuadExpr

# Constant
# Constant--Constant obviously already taken care of!
# Constant--VariableRef
Base.:+(lhs::Constant, rhs::AbstractVariableRef) = GenericAffExpr(_float(lhs), rhs => 1.0)
Base.:-(lhs::Constant, rhs::AbstractVariableRef) = GenericAffExpr(_float(lhs), rhs => -1.0)
Base.:*(lhs::Constant, rhs::AbstractVariableRef) = GenericAffExpr(0.0, rhs => _float(lhs))
# Constant--GenericAffExpr
function Base.:+(lhs::Constant, rhs::GenericAffExpr)
    result = copy(rhs)
    result.constant += lhs
    return result
end
function Base.:-(lhs::Constant, rhs::GenericAffExpr)
    result = -rhs
    result.constant += lhs
    return result
end
function Base.:*(lhs::Constant, rhs::GenericAffExpr)
    f = _float(lhs)
    return map_coefficients(c -> f * c, rhs)
end
# Constant--QuadExpr
Base.:+(lhs::Constant, rhs::GenericQuadExpr) = GenericQuadExpr(lhs+rhs.aff, copy(rhs.terms))
function Base.:-(lhs::Constant, rhs::GenericQuadExpr)
    result = -rhs
    result.aff.constant += lhs
    return result
end
Base.:*(lhs::Constant, rhs::GenericQuadExpr) = map_coefficients(c -> lhs * c, rhs)

# AbstractVariableRef (or, AbstractJuMPScalar)
# TODO: What is the role of AbstractJuMPScalar??
Base.:+(lhs::AbstractJuMPScalar) = lhs
Base.:-(lhs::AbstractVariableRef) = GenericAffExpr(0.0, lhs => -1.0)
Base.:*(lhs::AbstractJuMPScalar) = lhs # make this more generic so extensions don't have to define unary multiplication for our macros
# AbstractVariableRef--Constant
Base.:+(lhs::AbstractVariableRef, rhs::Constant) = (+)( rhs,lhs)
Base.:-(lhs::AbstractVariableRef, rhs::Constant) = (+)(-rhs,lhs)
Base.:*(lhs::AbstractVariableRef, rhs::Constant) = (*)(rhs,lhs)
Base.:/(lhs::AbstractVariableRef, rhs::Constant) = (*)(1.0/rhs,lhs)
# AbstractVariableRef--AbstractVariableRef
Base.:+(lhs::V, rhs::V) where {V <: AbstractVariableRef} = GenericAffExpr(0.0, lhs => 1.0, rhs => 1.0)
Base.:-(lhs::V, rhs::V) where {V <: AbstractVariableRef} = GenericAffExpr(0.0, lhs => 1.0, rhs => -1.0)
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
# GenericAffExpr--Constant
Base.:+(lhs::GenericAffExpr, rhs::Constant) = (+)(rhs,lhs)
Base.:-(lhs::GenericAffExpr, rhs::Constant) = (+)(-rhs,lhs)
Base.:*(lhs::GenericAffExpr, rhs::Constant) = (*)(rhs,lhs)
Base.:/(lhs::GenericAffExpr, rhs::Constant) = map_coefficients(c -> c/rhs, lhs)
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
Base.:^(lhs::Union{AbstractVariableRef, GenericAffExpr}, rhs::Constant) = error("Only exponents of 0, 1, or 2 are currently supported. Are you trying to build a nonlinear problem? Make sure you use @NLconstraint/@NLobjective.")
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
# GenericQuadExpr--Constant
Base.:+(lhs::GenericQuadExpr, rhs::Constant) = (+)(+rhs,lhs)
Base.:-(lhs::GenericQuadExpr, rhs::Constant) = (+)(-rhs,lhs)
Base.:*(lhs::GenericQuadExpr, rhs::Constant) = (*)(rhs,lhs)
Base.:/(lhs::GenericQuadExpr, rhs::Constant) = (*)(inv(rhs),lhs)
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

# TODO: specialize sum for DenseAxisArray and SparseAxisArray of JuMP objects?
function Base.sum(array::AbstractArray{<:AbstractVariableRef})
    result_expression = zero(eltype(array))
    for variable in array
        add_to_expression!(result_expression, variable)
    end
    return result_expression
end
# TODO: Specialize for iterables.
function Base.sum(affs::AbstractArray{T}) where {T <: AbstractJuMPScalar}
    new_aff = zero(T)
    for aff in affs
        add_to_expression!(new_aff, aff)
    end
    return new_aff
end

# Base Julia's generic fallback vecdot, aka dot, requires that dot, aka LinearAlgebra.dot, be defined
# for scalars, so instead of defining them one-by-one, we will
# fallback to the multiplication operator
LinearAlgebra.dot(lhs::_JuMPTypes, rhs::_JuMPTypes) = lhs*rhs
LinearAlgebra.dot(lhs::_JuMPTypes, rhs::Constant) = lhs*rhs
LinearAlgebra.dot(lhs::Constant, rhs::_JuMPTypes) = lhs*rhs

LinearAlgebra.dot(lhs::AbstractVector{T}, rhs::AbstractVector{S}) where {T <: _JuMPTypes, S <: _JuMPTypes} = _dot(lhs,rhs)
LinearAlgebra.dot(lhs::AbstractVector{T}, rhs::AbstractVector{S}) where {T <: _JuMPTypes, S} = _dot(lhs,rhs)
LinearAlgebra.dot(lhs::AbstractVector{T}, rhs::AbstractVector{S}) where {T, S <: _JuMPTypes} = _dot(lhs,rhs)

LinearAlgebra.dot(lhs::AbstractArray{T, N}, rhs::AbstractArray{S, N}) where {T <: _JuMPTypes, S, N} = _dot(lhs,rhs)
LinearAlgebra.dot(lhs::AbstractArray{T, N}, rhs::AbstractArray{S, N}) where {T <: _JuMPTypes, S <: _JuMPTypes, N} = _dot(lhs,rhs)
LinearAlgebra.dot(lhs::AbstractArray{T, N}, rhs::AbstractArray{S, N}) where {T, S <: _JuMPTypes, N} = _dot(lhs,rhs)

function _dot(lhs::AbstractArray{T}, rhs::AbstractArray{S}) where {T, S}
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

Base.promote_rule(V::Type{<:AbstractVariableRef}, R::Type{<:Real}) = GenericAffExpr{Float64, V}
Base.promote_rule(V::Type{<:AbstractVariableRef}, ::Type{<:GenericAffExpr{T}}) where {T} = GenericAffExpr{T, V}
Base.promote_rule(V::Type{<:AbstractVariableRef}, ::Type{<:GenericQuadExpr{T}}) where {T} = GenericQuadExpr{T, V}
Base.promote_rule(::Type{GenericAffExpr{S, V}}, R::Type{<:Real}) where {S, V} = GenericAffExpr{promote_type(S, R), V}
Base.promote_rule(::Type{<:GenericAffExpr{S, V}}, ::Type{<:GenericQuadExpr{T, V}}) where {S, T, V} = GenericQuadExpr{promote_type(S, T), V}
Base.promote_rule(::Type{GenericQuadExpr{S, V}}, R::Type{<:Real}) where {S, V} = GenericQuadExpr{promote_type(S, R), V}

Base.transpose(x::AbstractJuMPScalar) = x

# Can remove the following code once == overloading is removed

function LinearAlgebra.issymmetric(x::Matrix{T}) where {T <: _JuMPTypes}
    (n = size(x,1)) == size(x,2) || return false
    for i in 1:n, j in (i+1):n
        isequal(x[i,j], x[j,i]) || return false
    end
    true
end

# Special-case because the the base version wants to do fill!(::Array{AbstractVariableRef}, zero(GenericAffExpr{Float64,eltype(x)}))
_one_indexed(A) = all(x -> isa(x, Base.OneTo), axes(A))
function LinearAlgebra.diagm(x::AbstractVector{<:AbstractVariableRef})
    @assert _one_indexed(x) # Base.diagm doesn't work for non-one-indexed arrays in general.
    diagm(0=>copyto!(similar(x, GenericAffExpr{Float64, eltype(x)}), x))
end

###############################################################################
# Matrix/Vector Arithmetic with JuMP eltypes
###############################################################################

###############################################################################
# convenience/utility definitions

function _fill_with_zeros!(A)
    for I in eachindex(A)
        A[I] = zero(eltype(A))
    end
    return A
end

################################################################################
# `_mul!(ret, A, B)` sets the result of `A*B` into the buffer `ret`. We roll our
# own matmul here (instead of using Julia's generic fallbacks) because doing so
# allows us to accumulate the expressions for the inner loops in-place.
# Additionally, Julia's generic fallbacks can be finnicky when your array
# elements aren't `<:Number`.

const GenericAffOrQuadExpr{C, V} = Union{GenericAffExpr{C, V}, GenericQuadExpr{C, V}}

function _mul!(ret::AbstractVecOrMat{<:GenericAffOrQuadExpr}, A, B)
    size(A, 2) == size(B, 1) || throw(DimensionMismatch())
    size(A, 1) == size(ret, 1) || throw(DimensionMismatch())
    size(B, 2) == size(ret, 2) || throw(DimensionMismatch())
    _fill_with_zeros!(ret)
    for i ∈ 1:size(A, 1), j ∈ 1:size(B, 2)
        q = zero(eltype(ret))
        _sizehint_expr!(q, size(A, 2))
        for k ∈ 1:size(A, 2)
            add_to_expression!(q, A[i, k], B[k, j])
        end
        ret[i, j] = q
    end
    return ret
end

# This method of `mul!` is adapted from upstream Julia. Note that we
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

const TransposeOrAdjoint{T<:Union{GenericAffOrQuadExpr, Real}, MT} = Union{Transpose{T, MT}, Adjoint{T, MT}}
_mirror_transpose_or_adjoint(x, ::Transpose) = transpose(x)
_mirror_transpose_or_adjoint(x, ::Adjoint) = adjoint(x)

function _mul!(ret::AbstractVecOrMat{<:GenericAffOrQuadExpr},
               adjA::TransposeOrAdjoint{<:Any, <:SparseMatrixCSC},
               B::AbstractVecOrMat,
               α_expr=one(eltype(ret)),
               β=zero(eltype(ret)))
    A = parent(adjA)
    size(A, 2) == size(ret, 1) || throw(DimensionMismatch())
    size(A, 1) == size(B, 1) || throw(DimensionMismatch())
    size(B, 2) == size(ret, 2) || throw(DimensionMismatch())
    A_nonzeros = nonzeros(A)
    A_rowvals = rowvals(A)
    # See SparseArrays/src/linalg.jl
    if !isone(β)
        !iszero(β) ? rmul!(ret, β) : _fill_with_zeros!(ret)
    end
    α = convert(Float64, α_expr)
    for k ∈ 1:size(ret, 2)
        @inbounds for col ∈ 1:A.n
            tmp = zero(eltype(ret))
            for j ∈ A.colptr[col]:(A.colptr[col + 1] - 1)
                A_val = _mirror_transpose_or_adjoint(A_nonzeros[j], adjA)
                add_to_expression!(tmp, A_val, B[A_rowvals[j], k])
            end
            add_to_expression!(ret[col, k], α, tmp)
        end
    end
    return ret
end
function _mul!(ret::AbstractVecOrMat{<:GenericAffOrQuadExpr},
               A::SparseMatrixCSC, B::AbstractVecOrMat,
               α_expr=one(eltype(ret)),
               β=zero(eltype(ret)))
    size(A, 2) == size(B, 1) || throw(DimensionMismatch())
    size(A, 1) == size(ret, 1) || throw(DimensionMismatch())
    size(B, 2) == size(ret, 2) || throw(DimensionMismatch())
    A_nonzeros = nonzeros(A)
    A_rowvals = rowvals(A)
    # See SparseArrays/src/linalg.jl
    if !isone(β)
        !iszero(β) ? rmul!(ret, β) : _fill_with_zeros!(ret)
    end
    α = convert(Float64, α_expr)
    for col ∈ 1:size(A, 2)
        for k ∈ 1:size(ret, 2)
            αxj = α * B[col,k]
            for j ∈ nzrange(A, col)
                add_to_expression!(ret[A_rowvals[j], k], A_nonzeros[j], αxj)
            end
        end
    end
    return ret
end

function _mul!(ret::AbstractMatrix{<:GenericAffOrQuadExpr},
               A::AbstractMatrix, B::SparseMatrixCSC,
               α_expr=one(eltype(ret)),
               β=zero(eltype(ret)))
    rowval = rowvals(B)
    B_nonzeros = nonzeros(B)
    # See SparseArrays/src/linalg.jl
    if !isone(β)
        !iszero(β) ? rmul!(ret, β) : _fill_with_zeros!(ret)
    end
    α = convert(Float64, α_expr)
    for multivec_row in 1:size(A, 1)
        for col ∈ 1:size(B, 2)
            idxset = nzrange(B, col)
            tmp = zero(eltype(ret))
            _sizehint_expr!(tmp, length(idxset))
            for k ∈ idxset
                add_to_expression!(tmp, A[multivec_row, rowval[k]], B_nonzeros[k])
            end
            add_to_expression!(ret[multivec_row, col], α, tmp)
        end
    end
    return ret
end
function _mul!(ret::AbstractMatrix{<:GenericAffOrQuadExpr},
               A::AbstractMatrix,
               adjB::TransposeOrAdjoint{<:Any, <:SparseMatrixCSC},
               α_expr=one(eltype(ret)),
               β=zero(eltype(ret)))
    B = parent(adjB)
    B_rowvals = rowvals(B)
    B_nonzeros = nonzeros(B)
    # See SparseArrays/src/linalg.jl
    if !isone(β)
        !iszero(β) ? rmul!(ret, β) : _fill_with_zeros!(ret)
    end
    α = convert(Float64, α_expr)
    for B_col ∈ 1:size(B, 2), k ∈ nzrange(B, B_col)
        B_row = B_rowvals[k]
        B_val = _mirror_transpose_or_adjoint(B_nonzeros[k], adjB)
        αB_val = α * B_val
        for A_row in 1:size(A, 1)
            add_to_expression!(ret[A_row, B_row], A[A_row, B_col], αB_val)
        end
    end
    return ret
end


function _densify_with_jump_eltype(x::SparseMatrixCSC{V}) where {V <: AbstractVariableRef}
    return convert(Matrix{GenericAffExpr{Float64, V}}, x)
end
_densify_with_jump_eltype(x::AbstractMatrix) = convert(Matrix, x)

# TODO: Implement sparse * sparse code as in base/sparse/linalg.jl (spmatmul).
function _mul!(ret::AbstractMatrix{<:GenericAffOrQuadExpr},
               A::SparseMatrixCSC,
               B::SparseMatrixCSC)
    return mul!(ret, A, _densify_with_jump_eltype(B))
end

# TODO: Implement sparse * sparse code as in base/sparse/linalg.jl (spmatmul).
function _mul!(ret::AbstractMatrix{<:GenericAffOrQuadExpr},
               A::TransposeOrAdjoint{<:Any, <:SparseMatrixCSC},
               B::SparseMatrixCSC)
    return mul!(ret, A, _densify_with_jump_eltype(B))
end

# TODO: Implement sparse * sparse code as in base/sparse/linalg.jl (spmatmul).
function _mul!(ret::AbstractMatrix{<:GenericAffOrQuadExpr},
               A::SparseMatrixCSC,
               B::TransposeOrAdjoint{<:Any, <:SparseMatrixCSC})
    return mul!(ret, _densify_with_jump_eltype(A), B)
end

###############################################################################
# Interception of Base's matrix/vector arithmetic machinery

# Redirect calls with `eltype(ret) <: GenericAffOrQuadExpr` to `_mul!` to
# replace it with an implementation more efficient than `generic_matmatmul!` and
# `generic_matvecmul!` since it takes into account the mutability of the affine
# and quadratic JuMP expression.
# We need `args...` because SparseArrays` also gives `α` and `β` arguments.

function LinearAlgebra.mul!(ret::AbstractMatrix{<:GenericAffOrQuadExpr},
                            A::AbstractVecOrMat, B::AbstractVecOrMat, args...)
    return _mul!(ret, A, B, args...)
end
function LinearAlgebra.mul!(ret::AbstractVector{<:GenericAffOrQuadExpr},
                            A::AbstractVecOrMat, B::AbstractVector, args...)
    return _mul!(ret, A, B, args...)
end

function LinearAlgebra.mul!(ret::AbstractVector{<:GenericAffOrQuadExpr},
                            A::Transpose{<:Any,<:AbstractVecOrMat},
                            B::AbstractVector, args...)
    return _mul!(ret, A, B, args...)
end
function LinearAlgebra.mul!(ret::AbstractVector{<:GenericAffOrQuadExpr},
                            A::Adjoint{<:Any,<:AbstractVecOrMat},
                            B::AbstractVector, args...)
    return _mul!(ret, A, B, args...)
end
function LinearAlgebra.mul!(ret::AbstractMatrix{<:GenericAffOrQuadExpr},
                            A::Transpose{<:Any,<:AbstractVecOrMat},
                            B::AbstractMatrix, args...)
    return _mul!(ret, A, B, args...)
end
function LinearAlgebra.mul!(ret::AbstractMatrix{<:GenericAffOrQuadExpr},
                            A::Adjoint{<:Any,<:AbstractVecOrMat},
                            B::AbstractMatrix, args...)
    return _mul!(ret, A, B, args...)
end
function LinearAlgebra.mul!(ret::AbstractMatrix{<:GenericAffOrQuadExpr},
                            A::AbstractMatrix,
                            B::Transpose{<:Any,<:AbstractVecOrMat}, args...)
    return _mul!(ret, A, B, args...)
end
function LinearAlgebra.mul!(ret::AbstractMatrix{<:GenericAffOrQuadExpr},
                            A::AbstractMatrix,
                            B::Adjoint{<:Any,<:AbstractVecOrMat}, args...)
    return _mul!(ret, A, B, args...)
end

# SparseArrays promotes the element types of `A` and `B` to the same type
# which always produce quadratic expressions even if only one of them was
# affine and the other one constant. Moreover, it does not always go through
# `mul!` which prevents us from using mutability of JuMP affine and quadratic
# expression. For this reason we intercept the calls and redirect them to
# `_mul` which does that `LinearAlgebra/src/matmul.jl` does for abstract
# matrices and vector, i.e., use `matprod` to estimate the resulting element
# type, allocate the resulting array and redirect to `mul!`.

_A_mul_B_eltype(::Type{R}, ::Type{S}) where {R, S} = typeof(one(R) * one(S) + one(R) * one(S))
_A_mul_B_ret_dims(A::AbstractMatrix, B::AbstractVector) = (size(A, 1),)
_A_mul_B_ret_dims(A::AbstractMatrix, B::AbstractMatrix) = (size(A, 1), size(B, 2))
function _mul(A::AbstractVecOrMat, B::AbstractVecOrMat)
    T = _A_mul_B_eltype(eltype(A), eltype(B))
    ret = Array{T}(undef, _A_mul_B_ret_dims(A, B))
    return mul!(ret, A, B)
end

# A few are overwritten below but many more need to be redirected to `_mul` in
# `linalg.jl`.
Base.:*(A::SparseMatrixCSC{<:AbstractJuMPScalar}, B::SparseMatrixCSC{<:AbstractJuMPScalar}) = _mul(A, B)
Base.:*(A::SparseMatrixCSC{<:Any}, B::SparseMatrixCSC{<:AbstractJuMPScalar}) = _mul(A, B)
Base.:*(A::SparseMatrixCSC{<:AbstractJuMPScalar}, B::SparseMatrixCSC{<:Any}) = _mul(A, B)

Base.:*(A::SparseMatrixCSC{<:AbstractJuMPScalar}, B::Adjoint{<:AbstractJuMPScalar, <:SparseMatrixCSC}) = _mul(A, B)
Base.:*(A::SparseMatrixCSC{<:Any}, B::Adjoint{<:AbstractJuMPScalar, <:SparseMatrixCSC}) = _mul(A, B)
Base.:*(A::SparseMatrixCSC{<:AbstractJuMPScalar}, B::Adjoint{<:Any, <:SparseMatrixCSC}) = _mul(A, B)

Base.:*(A::Adjoint{<:AbstractJuMPScalar, <:SparseMatrixCSC}, B::SparseMatrixCSC{<:AbstractJuMPScalar}) = _mul(A, B)
Base.:*(A::Adjoint{<:Any, <:SparseMatrixCSC}, B::SparseMatrixCSC{<:AbstractJuMPScalar}) = _mul(A, B)
Base.:*(A::Adjoint{<:AbstractJuMPScalar, <:SparseMatrixCSC}, B::SparseMatrixCSC{<:Any}) = _mul(A, B)

Base.:*(A::StridedMatrix{<:AbstractJuMPScalar}, B::SparseMatrixCSC{<:AbstractJuMPScalar}) = _mul(A, B)
Base.:*(A::StridedMatrix{<:Any}, B::SparseMatrixCSC{<:AbstractJuMPScalar}) = _mul(A, B)
Base.:*(A::StridedMatrix{<:AbstractJuMPScalar}, B::SparseMatrixCSC{<:Any}) = _mul(A, B)

Base.:*(A::SparseMatrixCSC{<:AbstractJuMPScalar}, B::StridedMatrix{<:AbstractJuMPScalar}) = _mul(A, B)
Base.:*(A::SparseMatrixCSC{<:Any}, B::StridedMatrix{<:AbstractJuMPScalar}) = _mul(A, B)
Base.:*(A::SparseMatrixCSC{<:AbstractJuMPScalar}, B::StridedMatrix{<:Any}) = _mul(A, B)

Base.:*(A::Adjoint{<:AbstractJuMPScalar, <:SparseMatrixCSC}, B::StridedMatrix{<:AbstractJuMPScalar}) = _mul(A, B)
Base.:*(A::Adjoint{<:Any, <:SparseMatrixCSC}, B::StridedMatrix{<:AbstractJuMPScalar}) = _mul(A, B)
Base.:*(A::Adjoint{<:AbstractJuMPScalar, <:SparseMatrixCSC}, B::StridedMatrix{<:Any}) = _mul(A, B)

# TODO: Intercepting "externally owned" method calls by dispatching on type parameters
# (rather than outermost wrapper type) is generally bad practice, but refactoring this code
# to use a different mechanism would be a lot of work. In the future, this interception code
# would be more easily/robustly replaced by using a tool like
# https://github.com/jrevels/Cassette.jl.

# Base doesn't define efficient fallbacks for sparse array arithmetic involving
# non-`<:Number` scalar elements, so we define some of these for `<:JuMPType` scalar
# elements here.

function Base.:*(A::Number, B::SparseMatrixCSC{T}) where {T <: _JuMPTypes}
    return SparseMatrixCSC(B.m, B.n, copy(B.colptr), copy(rowvals(B)), A .* nonzeros(B))
end

function Base.:*(A::SparseMatrixCSC{T}, B::Number) where {T <: _JuMPTypes}
    return SparseMatrixCSC(A.m, A.n, copy(A.colptr), copy(rowvals(A)), nonzeros(A) .* B)
end

function Base.:*(A::_JuMPTypes, B::SparseMatrixCSC)
    return SparseMatrixCSC(B.m, B.n, copy(B.colptr), copy(rowvals(B)), A .* nonzeros(B))
end

function Base.:*(A::SparseMatrixCSC, B::_JuMPTypes)
    return SparseMatrixCSC(A.m, A.n, copy(A.colptr), copy(rowvals(A)), nonzeros(A) .* B)
end

function Base.:/(A::SparseMatrixCSC{T}, B::Number) where {T <: _JuMPTypes}
    return SparseMatrixCSC(A.m, A.n, copy(A.colptr), copy(rowvals(A)), nonzeros(A) ./ B)
end

Base.:*(x::AbstractArray{T}) where {T <: _JuMPTypes} = x

Base.:+(x::AbstractArray{T}) where {T <: _JuMPTypes} = x

function Base.:-(x::AbstractArray{T}) where {T <: _JuMPTypes}
    ret = similar(x, typeof(-one(T)))
    for I in eachindex(ret)
        ret[I] = -x[I]
    end
    return ret
end

# Fix https://github.com/JuliaLang/julia/issues/32374 as done in
# https://github.com/JuliaLang/julia/pull/32375. This hack should
# be removed once we drop Julia v1.0.
function Base.:-(A::Symmetric{<:JuMP.AbstractVariableRef})
    return Symmetric(-A.data, LinearAlgebra.sym_uplo(A.uplo))
end
function Base.:-(A::Hermitian{<:JuMP.AbstractVariableRef})
    return Hermitian(-A.data, LinearAlgebra.sym_uplo(A.uplo))
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
Base.:/(::S, ::T) where {S <: Union{Constant, AbstractVariableRef, GenericAffExpr, GenericQuadExpr},
                         T <: Union{AbstractVariableRef, GenericAffExpr, GenericQuadExpr}} =
    error( "/(::$S,::$T) is not defined. $op_hint")

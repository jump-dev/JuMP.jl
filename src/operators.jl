#  Copyright 2017, Iain Dunning, Joey Huchette, Miles Lubin, and contributors
#  This Source Code Form is subject to the terms of the Mozilla Public
#  License, v. 2.0. If a copy of the MPL was not distributed with this
#  file, You can obtain one at http://mozilla.org/MPL/2.0/.
#############################################################################
# JuMP
# An algebraic modeling language for Julia
# See http://github.com/JuliaOpt/JuMP.jl
#############################################################################

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
    sizehint!(result, length(linearterms(rhs)) + 1)
    add_to_expression!(result, one(C), lhs)
    for (coef, var) in linearterms(rhs)
        add_to_expression!(result, coef, var)
    end
    return result
end

function Base.:-(lhs::V, rhs::GenericAffExpr{C,V}) where {C,V <: AbstractVariableRef}
    # For the variables to have the proper order in the result, we need to add the lhs first.
    result = zero(rhs)
    result.constant = -rhs.constant
    sizehint!(result, length(linearterms(rhs)) + 1)
    add_to_expression!(result, one(C), lhs)
    for (coef, var) in linearterms(rhs)
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
    for (coef, var) in linearterms(rhs)
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
        return convert(GenericQuadExpr{Float64, variablereftype(lhs)}, lhs)
    elseif rhs == 0
        return one(GenericQuadExpr{Float64, variablereftype(lhs)})
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
    for (coef, var) in linearterms(lhs)
        add_to_expression!(result, coef, var, rhs)
    end
    return result
end
Base.:/(lhs::GenericAffExpr, rhs::AbstractVariableRef) = error("Cannot divide affine expression by a variable")
# AffExpr--AffExpr
function Base.:+(lhs::GenericAffExpr{C,V}, rhs::GenericAffExpr{C,V}) where {C,V<:JuMPTypes}
    operator_warn(lhs,rhs)
    result_terms = copy(lhs.terms)
    # merge() returns a Dict(), so we need to call merge!() instead.
    # Note: merge!() doesn't appear to call sizehint!(). Is this important?
    merge!(+, result_terms, rhs.terms)
    return GenericAffExpr(lhs.constant + rhs.constant, result_terms)
end

function Base.:-(lhs::GenericAffExpr{C,V}, rhs::GenericAffExpr{C,V}) where {C,V<:JuMPTypes}
    result = copy(lhs)
    result.constant -= rhs.constant
    sizehint!(result, length(linearterms(lhs)) + length(linearterms(rhs)))
    for (coef, var) in linearterms(rhs)
        add_to_expression!(result, -coef, var)
    end
    return result
end

function Base.:*(lhs::GenericAffExpr{C,V}, rhs::GenericAffExpr{C,V}) where {C,V<:JuMPTypes}
    result = zero(GenericQuadExpr{C,V})

    lhs_length = length(linearterms(lhs))
    rhs_length = length(linearterms(rhs))

    # Quadratic terms
    for (lhscoef, lhsvar) in linearterms(lhs)
        for (rhscoef, rhsvar) in linearterms(rhs)
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
        for (rhscoef, rhsvar) in linearterms(rhs)
            add_to_expression!(result.aff, c*rhscoef, rhsvar)
        end
    end

    # [RHS constant] * [LHS linear terms]
    if !iszero(rhs.constant)
        c = rhs.constant
        for (lhscoef, lhsvar) in linearterms(lhs)
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
    for (coef, var) in linearterms(q2)
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
    for (coef, var) in linearterms(q2)
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

# Base Julia's generic fallback vecdot requires that dot be defined
# for scalars, so instead of defining them one-by-one, we will
# fallback to the multiplication operator
Compat.LinearAlgebra.dot(lhs::JuMPTypes, rhs::JuMPTypes) = lhs*rhs
Compat.LinearAlgebra.dot(lhs::JuMPTypes, rhs::Number)    = lhs*rhs
Compat.LinearAlgebra.dot(lhs::Number,    rhs::JuMPTypes) = lhs*rhs

Compat.LinearAlgebra.dot(lhs::AbstractVector{T},rhs::AbstractVector{S}) where {T<:JuMPTypes,S<:JuMPTypes} = _dot(lhs,rhs)
Compat.LinearAlgebra.dot(lhs::AbstractVector{T},rhs::AbstractVector{S}) where {T<:JuMPTypes,S} = _dot(lhs,rhs)
Compat.LinearAlgebra.dot(lhs::AbstractVector{T},rhs::AbstractVector{S}) where {T,S<:JuMPTypes} = _dot(lhs,rhs)

Compat.LinearAlgebra.vecdot(lhs::AbstractArray{T,N},rhs::AbstractArray{S,N}) where {T<:JuMPTypes,S,N} = _dot(lhs,rhs)
Compat.LinearAlgebra.vecdot(lhs::AbstractArray{T,N},rhs::AbstractArray{S,N}) where {T<:JuMPTypes,S<:JuMPTypes,N} = _dot(lhs,rhs)
Compat.LinearAlgebra.vecdot(lhs::AbstractArray{T,N},rhs::AbstractArray{S,N}) where {T,S<:JuMPTypes,N} = _dot(lhs,rhs)

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
function Compat.LinearAlgebra.diagm(x::AbstractVector{<:AbstractVariableRef})
    @assert one_indexed(x) # Base.diagm doesn't work for non-one-indexed arrays in general.
    diagm(copy!(similar(x, GenericAffExpr{Float64,eltype(x)}), x))
end

###############
# The _multiply!(buf,y,z) adds the results of y*z into the buffer buf. No bounds/size
# checks are performed; it is expected that the caller has done this, has ensured
# that the eltype of buf is appropriate, and has zeroed the elements of buf (if desired).

function _multiply!(ret::Array{T}, lhs::Array, rhs::Array) where T<:JuMPTypes
    m, n = size(lhs,1), size(lhs,2)
    r, s = size(rhs,1), size(rhs,2)
    for i ∈ 1:m, j ∈ 1:s
        q = ret[i,j]
        _sizehint_expr!(q, n)
        for k ∈ 1:n
            tmp = convert(T, lhs[i,k]*rhs[k,j])
            add_to_expression!(q, tmp)
        end
    end
    ret
end

# this computes lhs.'*rhs and places it in ret
function _multiplyt!(ret::Array{T}, lhs::Array, rhs::Array) where T<:JuMPTypes
    m, n = size(lhs,2), size(lhs,1) # transpose
    r, s = size(rhs,1), size(rhs,2)
    for i ∈ 1:m, j ∈ 1:s
        q = ret[i,j]
        _sizehint_expr!(q, n)
        for k ∈ 1:n
            tmp = convert(T, lhs[k,i]*rhs[k,j]) # transpose
            add_to_expression!(q, tmp)
        end
    end
    ret
end

function _multiply!(ret::Array{T}, lhs::SparseMatrixCSC, rhs::Array) where T<:Union{GenericAffExpr,GenericQuadExpr}
    nzv = nonzeros(lhs)
    rv  = rowvals(lhs)
    for col ∈ 1:lhs.n
        for k ∈ 1:size(ret, 2)
            for j ∈ nzrange(lhs, col)
                add_to_expression!(ret[rv[j], k], nzv[j] * rhs[col,k])
            end
        end
    end
    ret
end

# this computes lhs.'*rhs and places it in ret
function _multiplyt!(ret::Array{T}, lhs::SparseMatrixCSC, rhs::Array) where T<:Union{GenericAffExpr,GenericQuadExpr}
    _multiply!(ret, transpose(lhs), rhs) # TODO fully implement
end

function _multiply!(ret::Array{T}, lhs::Matrix, rhs::SparseMatrixCSC) where T<:Union{GenericAffExpr,GenericQuadExpr}
    rowval = rowvals(rhs)
    nzval  = nonzeros(rhs)
    for multivec_row in 1:size(lhs,1)
        for col ∈ 1:rhs.n
            idxset = nzrange(rhs, col)
            q = ret[multivec_row, col]
            _sizehint_expr!(q, length(idxset))
            for k ∈ idxset
                add_to_expression!(q, lhs[multivec_row, rowval[k]] * nzval[k])
            end
        end
    end
    ret
end

# this computes lhs.'*rhs and places it in ret
function _multiplyt!(ret::Array{T}, lhs::Matrix, rhs::SparseMatrixCSC) where T<:Union{GenericAffExpr,GenericQuadExpr}
    rowval = rowvals(rhs)
    nzval  = nonzeros(rhs)
    for multivec_row ∈ 1:size(lhs,2) # transpose
        for col ∈ 1:rhs.n
            idxset = nzrange(rhs, col)
            q = ret[multivec_row, col]
            _sizehint_expr!(q, length(idxset))
            for k ∈ idxset
                add_to_expression!(q, lhs[rowval[k], multivec_row] * nzval[k]) # transpose
            end
        end
    end
    ret
end

# See https://github.com/JuliaLang/julia/issues/27015
function Base.Matrix(S::SparseMatrixCSC{V}) where V<:AbstractVariableRef
    A = zeros(GenericAffExpr{Float64, V}, S.m, S.n)
    for Sj in 1:S.n
        for Sk in nzrange(S, Sj)
            Si = S.rowval[Sk]
            Sv = S.nzval[Sk]
            A[Si, Sj] = Sv
        end
    end
    return A
end
# TODO: implement sparse * sparse code as in base/sparse/linalg.jl (spmatmul)
_multiply!(ret::AbstractArray{T}, lhs::SparseMatrixCSC, rhs::SparseMatrixCSC) where {T<:JuMPTypes} = _multiply!(ret, lhs, Matrix(rhs))
_multiplyt!(ret::AbstractArray{T}, lhs::SparseMatrixCSC, rhs::SparseMatrixCSC) where {T<:JuMPTypes} = _multiplyt!(ret, lhs, Matrix(rhs))

function Base.:*(A::Union{Matrix{T}, SparseMatrixCSC{T}},
                 x::Union{Matrix, Vector, SparseMatrixCSC}
                 ) where {T <: JuMPTypes}
    return _matmul(A, x)
end
function Base.:*(A::Union{Matrix{T}, SparseMatrixCSC{T}},
                 x::Union{Matrix{R}, Vector{R}, SparseMatrixCSC{R}}
                 ) where {T <: JuMPTypes, R <: JuMPTypes}
    return _matmul(A, x)
end
function Base.:*(A::Union{Matrix, SparseMatrixCSC},
                 x::Union{Matrix{T}, Vector{T}, SparseMatrixCSC{T}}
                 ) where {T <: JuMPTypes}
    return _matmul(A, x)
end

if VERSION >= v"0.7-"
    # This is a stopgap solution to get (most) tests passing on Julia 0.7. A lot
    # of cases probably still don't work. (Like A * A where A is a sparse matrix
    # of a JuMP type). This code needs a big refactor.
    function Base.:*(adjA::Adjoint{<:JuMPTypes,<:SparseMatrixCSC},
                     x::Vector)
        A = adjA.parent
        ret = _return_arrayt(A, x)
        return mul!(ret, adjoint(A), x)
    end
    function Base.:*(adjA::Adjoint{<:Any,<:SparseMatrixCSC},
                     x::Vector{<:JuMPTypes})
        A = adjA.parent
        ret = _return_arrayt(A, x)
        return mul!(ret, adjoint(A), x)
    end
    function Base.:*(adjA::Adjoint{<:JuMPTypes,<:SparseMatrixCSC},
                     x::Vector{<:JuMPTypes})
        A = adjA.parent
        ret = _return_arrayt(A, x)
        return mul!(ret, adjoint(A), x)
    end
    function Base.:*(adjA::Transpose{<:JuMPTypes,<:SparseMatrixCSC},
                     x::Vector)
        A = adjA.parent
        ret = _return_arrayt(A, x)
        return mul!(ret, transpose(A), x)
    end
    function Base.:*(adjA::Transpose{<:Any,<:SparseMatrixCSC},
                     x::Vector{<:JuMPTypes})
        A = adjA.parent
        ret = _return_arrayt(A, x)
        return mul!(ret, transpose(A), x)
    end
    function Base.:*(adjA::Transpose{<:JuMPTypes,<:SparseMatrixCSC},
                     x::Vector{<:JuMPTypes})
        A = adjA.parent
        ret = _return_arrayt(A, x)
        return mul!(ret, transpose(A), x)
    end
    # Matrix versions.
    function Base.:*(adjA::Adjoint{<:JuMPTypes,<:SparseMatrixCSC},
                     x::Matrix)
        A = adjA.parent
        ret = _return_arrayt(A, x)
        return mul!(ret, adjoint(A), x)
    end
    function Base.:*(adjA::Adjoint{<:Any,<:SparseMatrixCSC},
                     x::Matrix{<:JuMPTypes})
        A = adjA.parent
        ret = _return_arrayt(A, x)
        return mul!(ret, adjoint(A), x)
    end
    function Base.:*(adjA::Adjoint{<:JuMPTypes,<:SparseMatrixCSC},
                     x::Matrix{<:JuMPTypes})
        A = adjA.parent
        ret = _return_arrayt(A, x)
        return mul!(ret, adjoint(A), x)
    end
    function Base.:*(adjA::Transpose{<:JuMPTypes,<:SparseMatrixCSC},
                     x::Matrix)
        A = adjA.parent
        ret = _return_arrayt(A, x)
        return mul!(ret, transpose(A), x)
    end
    function Base.:*(adjA::Transpose{<:Any,<:SparseMatrixCSC},
                     x::Matrix{<:JuMPTypes})
        A = adjA.parent
        ret = _return_arrayt(A, x)
        return mul!(ret, transpose(A), x)
    end
    function Base.:*(adjA::Transpose{<:JuMPTypes,<:SparseMatrixCSC},
                     x::Matrix{<:JuMPTypes})
        A = adjA.parent
        ret = _return_arrayt(A, x)
        return mul!(ret, transpose(A), x)
    end
    # mul! is adapted from upstream Julia.
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
    # We confuse transpose with adjoint because they're the same for all JuMP
    # types. Note this doesn't extend the LinearAlgebra version. It assumes
    # that C is already filled with zeros.
    function mul!(C::StridedVecOrMat,
                                       adjA::Union{Adjoint{<:Any,<:SparseMatrixCSC},
                                                   Transpose{<:Any,<:SparseMatrixCSC}},
                                       B::StridedVecOrMat)
        A = adjA.parent
        A.n == size(C, 1) || throw(DimensionMismatch())
        A.m == size(B, 1) || throw(DimensionMismatch())
        size(B, 2) == size(C, 2) || throw(DimensionMismatch())
        nzv = A.nzval
        rv = A.rowval
        # C is already filled with zeros by _return_arrayt.
        for k = 1:size(C, 2)
            @inbounds for col = 1:A.n
                tmp = zero(eltype(C))
                for j = A.colptr[col]:(A.colptr[col + 1] - 1)
                    tmp += adjoint(nzv[j])*B[rv[j],k]
                end
                C[col,k] += tmp
            end
        end
        C
    end
else
    _multiply!(ret, lhs, rhs) = A_mul_B!(ret, lhs, rhs)
    _multiplyt!(ret, lhs, rhs) = At_mul_B!(ret, lhs, rhs)

    # these methods are called when one does A.'*v or A'*v respectively
    Base.At_mul_B(A::Union{Matrix{T},SparseMatrixCSC{T}}, x::Union{Matrix, Vector, SparseMatrixCSC}) where {T<:JuMPTypes} = _matmult(A, x)
    Base.At_mul_B(A::Union{Matrix{T},SparseMatrixCSC{T}}, x::Union{Matrix{R}, Vector{R}, SparseMatrixCSC{R}}) where {T<:JuMPTypes,R<:JuMPTypes} = _matmult(A, x)
    Base.At_mul_B(A::Union{Matrix,SparseMatrixCSC}, x::Union{Matrix{T}, Vector{T}, SparseMatrixCSC{T}}) where {T<:JuMPTypes} = _matmult(A, x)
    # these methods are the same as above since complex numbers are not implemented in JuMP
    Base.Ac_mul_B(A::Union{Matrix{T},SparseMatrixCSC{T}}, x::Union{Matrix, Vector, SparseMatrixCSC}) where {T<:JuMPTypes} = _matmult(A, x)
    Base.Ac_mul_B(A::Union{Matrix{T},SparseMatrixCSC{T}}, x::Union{Matrix{R}, Vector{R}, SparseMatrixCSC{R}}) where {T<:JuMPTypes,R<:JuMPTypes} = _matmult(A, x)
    Base.Ac_mul_B(A::Union{Matrix,SparseMatrixCSC}, x::Union{Matrix{T}, Vector{T}, SparseMatrixCSC{T}}) where {T<:JuMPTypes} = _matmult(A, x)
end

function _matmul(A, x)
    m, n = size(A,1), size(A,2)
    r, s = size(x,1), size(x,2)
    n == r || error("Incompatible sizes")
    ret = _return_array(A, x)
    _multiply!(ret, A, x)
    ret
end

function _matmult(A, x)
    m, n = size(A,2), size(A,1) # transpose
    r, s = size(x,1), size(x,2)
    n == r || error("Incompatible sizes")
    ret = _return_arrayt(A, x)
    _multiplyt!(ret, A, x)
    ret
end

# See https://github.com/JuliaLang/julia/pull/18218
_matprod_type(R, S) = typeof(one(R) * one(S) + one(R) * one(S))
# Don't do size checks here in _return_array, defer that to (*)
_return_array(A::AbstractMatrix{R}, x::AbstractVector{S}) where {R,S} = _fillwithzeros(Array{_matprod_type(R,S)}(undef,size(A,1)))
_return_array(A::AbstractMatrix{R}, x::AbstractMatrix{S}) where {R,S} = _fillwithzeros(Array{_matprod_type(R,S)}(undef,size(A,1), size(x,2)))
# these are for transpose return matrices
_return_arrayt(A::AbstractMatrix{R}, x::AbstractVector{S}) where {R,S} = _fillwithzeros(Array{_matprod_type(R,S)}(undef,size(A,2)))
_return_arrayt(A::AbstractMatrix{R}, x::AbstractMatrix{S}) where {R,S} = _fillwithzeros(Array{_matprod_type(R,S)}(undef,size(A,2), size(x, 2)))

# helper so we don't fill the buffer array with the same object
function _fillwithzeros(arr::AbstractArray{T}) where T
    for I in eachindex(arr)
        arr[I] = zero(T)
    end
    arr
end

# Special-case sparse matrix scalar multiplication/division
Base.:*(lhs::Number, rhs::SparseMatrixCSC{T}) where {T<:JuMPTypes} =
    SparseMatrixCSC(rhs.m, rhs.n, copy(rhs.colptr), copy(rhs.rowval), lhs .* rhs.nzval)
Base.:*(lhs::JuMPTypes, rhs::SparseMatrixCSC) =
    SparseMatrixCSC(rhs.m, rhs.n, copy(rhs.colptr), copy(rhs.rowval), lhs .* rhs.nzval)
Base.:*(lhs::SparseMatrixCSC{T}, rhs::Number) where {T<:JuMPTypes} =
    SparseMatrixCSC(lhs.m, lhs.n, copy(lhs.colptr), copy(lhs.rowval), lhs.nzval .* rhs)
Base.:*(lhs::SparseMatrixCSC, rhs::JuMPTypes) =
    SparseMatrixCSC(lhs.m, lhs.n, copy(lhs.colptr), copy(lhs.rowval), lhs.nzval .* rhs)
Base.:/(lhs::SparseMatrixCSC{T}, rhs::Number) where {T<:JuMPTypes} =
    SparseMatrixCSC(lhs.m, lhs.n, copy(lhs.colptr), copy(lhs.rowval), lhs.nzval ./ rhs)


for (op,opsymbol) in [(+,:+), (-,:-), (*,:*), (/,:/)]
    @eval begin
        Base.broadcast(::typeof($op),lhs::Number,rhs::JuMPTypes) = $opsymbol(lhs,rhs)
        Base.broadcast(::typeof($op),lhs::JuMPTypes,rhs::Number) = $opsymbol(lhs,rhs)
    end
end


Base.:+(x::AbstractArray{T}) where {T<:JuMPTypes} = x
function Base.:-(x::AbstractArray{T}) where T<:JuMPTypes
    ret = similar(x, typeof(-one(T)))
    for I in eachindex(ret)
        ret[I] = -x[I]
    end
    ret
end
Base.:*(x::AbstractArray{T}) where {T<:JuMPTypes} = x

###############################################################################
# Add nonlinear function fallbacks for JuMP built-in types
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

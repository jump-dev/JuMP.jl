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
# 2. Variable
# 4. [Generic]AffExpr
# 5. QuadExpr

# Number
# Number--Number obviously already taken care of!
# Number--Variable
(+)(lhs::Number, rhs::Variable) = AffExpr([rhs],[+1.],convert(Float64,lhs))
(-)(lhs::Number, rhs::Variable) = AffExpr([rhs],[-1.],convert(Float64,lhs))
(*)(lhs::Number, rhs::Variable) = AffExpr([rhs],[convert(Float64,lhs)], 0.)
# Number--GenericAffExpr
(+)(lhs::Number, rhs::GenericAffExpr) = GenericAffExpr(copy(rhs.vars),copy(rhs.coeffs),lhs+rhs.constant)
(-)(lhs::Number, rhs::GenericAffExpr) = GenericAffExpr(copy(rhs.vars),    -rhs.coeffs ,lhs-rhs.constant)
(*)(lhs::Number, rhs::GenericAffExpr) = GenericAffExpr(copy(rhs.vars),[lhs*rhs.coeffs[i] for i=1:length(rhs.coeffs)],lhs*rhs.constant)
# Number--QuadExpr
(+)(lhs::Number, rhs::QuadExpr) = QuadExpr(copy(rhs.qvars1),copy(rhs.qvars2),copy(rhs.qcoeffs),lhs+rhs.aff)
(-)(lhs::Number, rhs::QuadExpr) = QuadExpr(copy(rhs.qvars1),copy(rhs.qvars2),    -rhs.qcoeffs ,lhs-rhs.aff)
(*)(lhs::Number, rhs::QuadExpr) = QuadExpr(copy(rhs.qvars1),copy(rhs.qvars2), lhs*rhs.qcoeffs ,lhs*rhs.aff)

# Variable (or, AbstractJuMPScalar)
(+)(lhs::AbstractJuMPScalar) = lhs
(-)(lhs::Variable) = AffExpr([lhs],[-1.0],0.0)
(*)(lhs::AbstractJuMPScalar) = lhs # make this more generic so extensions don't have to define unary multiplication for our macros
# Variable--Number
(+)(lhs::Variable, rhs::Number) = (+)( rhs,lhs)
(-)(lhs::Variable, rhs::Number) = (+)(-rhs,lhs)
(*)(lhs::Variable, rhs::Number) = (*)(rhs,lhs)
(/)(lhs::Variable, rhs::Number) = (*)(1./rhs,lhs)
# Variable--Variable
(+)(lhs::Variable, rhs::Variable) = AffExpr([lhs,rhs], [1.,+1.], 0.)
(-)(lhs::Variable, rhs::Variable) = AffExpr([lhs,rhs], [1.,-1.], 0.)
(*)(lhs::Variable, rhs::Variable) = QuadExpr([lhs],[rhs],[1.],AffExpr(Variable[],Float64[],0.))
# Variable--AffExpr
(+)(lhs::V, rhs::GenericAffExpr{C,V}) where {C,V<:JuMPTypes} =
    GenericAffExpr{C,V}(vcat(lhs,rhs.vars),vcat(one(C),rhs.coeffs), rhs.constant)
(-)(lhs::V, rhs::GenericAffExpr{C,V}) where {C,V<:JuMPTypes} =
    GenericAffExpr{C,V}(vcat(lhs,rhs.vars),vcat(one(C),-rhs.coeffs),-rhs.constant)
function (*)(lhs::Variable, rhs::AffExpr)
    n = length(rhs.vars)
    if rhs.constant != 0.
        ret = QuadExpr([lhs for i=1:n],copy(rhs.vars),copy(rhs.coeffs),AffExpr([lhs], [rhs.constant], 0.))
    else
        ret = QuadExpr([lhs for i=1:n],copy(rhs.vars),copy(rhs.coeffs),zero(AffExpr))
    end
end
(/)(lhs::Variable, rhs::AffExpr) = error("Cannot divide a variable by an affine expression")
# Variable--QuadExpr
(+)(v::Variable, q::QuadExpr) = QuadExpr(copy(q.qvars1),copy(q.qvars2),copy(q.qcoeffs),v+q.aff)
(-)(v::Variable, q::QuadExpr) = QuadExpr(copy(q.qvars1),copy(q.qvars2),    -q.qcoeffs ,v-q.aff)

# GenericAffExpr
(+)(lhs::GenericAffExpr) = lhs
(-)(lhs::GenericAffExpr) = GenericAffExpr(lhs.vars, -lhs.coeffs, -lhs.constant)
(*)(lhs::GenericAffExpr) = lhs
# AffExpr--Number
(+)(lhs::GenericAffExpr, rhs::Number) = (+)(+rhs,lhs)
(-)(lhs::GenericAffExpr, rhs::Number) = (+)(-rhs,lhs)
(*)(lhs::GenericAffExpr, rhs::Number) = (*)(rhs,lhs)
(/)(lhs::GenericAffExpr, rhs::Number) = (*)(1.0/rhs,lhs)
function (^)(lhs::Union{Variable,AffExpr}, rhs::Integer)
    if rhs == 2
        return lhs*lhs
    elseif rhs == 1
        return QuadExpr(lhs)
    elseif rhs == 0
        return QuadExpr(1)
    else
        error("Only exponents of 0, 1, or 2 are currently supported. Are you trying to build a nonlinear problem? Make sure you use @NLconstraint/@NLobjective.")
    end
end
(^)(lhs::Union{Variable,AffExpr}, rhs::Number) = error("Only exponents of 0, 1, or 2 are currently supported. Are you trying to build a nonlinear problem? Make sure you use @NLconstraint/@NLobjective.")
# AffExpr--Variable
(+)(lhs::GenericAffExpr{C,V}, rhs::V) where {C,V<:JuMPTypes} = GenericAffExpr{C,V}(vcat(lhs.vars,rhs),vcat(lhs.coeffs,one(C)), lhs.constant)
(-)(lhs::GenericAffExpr{C,V}, rhs::V) where {C,V<:JuMPTypes} = GenericAffExpr{C,V}(vcat(lhs.vars,rhs),vcat(lhs.coeffs,-one(C)),lhs.constant)
# Don't fall back on Variable*AffExpr to preserve lhs/rhs consistency (appears in printing)
function (*)(lhs::AffExpr, rhs::Variable)
    n = length(lhs.vars)
    if lhs.constant != 0.
        ret = QuadExpr(copy(lhs.vars),[rhs for i=1:n],copy(lhs.coeffs),AffExpr([rhs], [lhs.constant], 0.))
    else
        ret = QuadExpr(copy(lhs.vars),[rhs for i=1:n],copy(lhs.coeffs),zero(AffExpr))
    end
end
(/)(lhs::AffExpr, rhs::Variable) = error("Cannot divide affine expression by a variable")
# AffExpr--AffExpr
(+)(lhs::GenericAffExpr{C,V}, rhs::GenericAffExpr{C,V}) where {C,V<:JuMPTypes} = (operator_warn(lhs,rhs); GenericAffExpr(vcat(lhs.vars,rhs.vars),vcat(lhs.coeffs, rhs.coeffs),lhs.constant+rhs.constant))
(-)(lhs::GenericAffExpr{C,V}, rhs::GenericAffExpr{C,V}) where {C,V<:JuMPTypes} = GenericAffExpr(vcat(lhs.vars,rhs.vars),vcat(lhs.coeffs,-rhs.coeffs),lhs.constant-rhs.constant)
function (*)(lhs::AffExpr, rhs::AffExpr)
    ret = QuadExpr(Variable[],Variable[],Float64[],AffExpr(Variable[],Float64[],0.))

    # Quadratic terms
    n = length(lhs.coeffs)
    m = length(rhs.coeffs)
    sizehint!(ret.qvars1, n*m)
    sizehint!(ret.qvars2, n*m)
    sizehint!(ret.qcoeffs, n*m)
    for i = 1:n
        for j = 1:m
            push!(ret.qvars1,  lhs.vars[i])
            push!(ret.qvars2,  rhs.vars[j])
            push!(ret.qcoeffs, lhs.coeffs[i]*rhs.coeffs[j])
        end
    end

    # Try to preallocate space for aff
    if lhs.constant != 0 && rhs.constant != 0
        sizehint!(ret.aff.vars, n+m)
        sizehint!(ret.aff.coeffs, n+m)
    elseif lhs.constant != 0
        sizehint!(ret.aff.vars, n)
        sizehint!(ret.aff.coeffs, n)
    elseif rhs.constant != 0
        sizehint!(ret.aff.vars, m)
        sizehint!(ret.aff.coeffs, m)
    end

    # [LHS constant] * RHS
    if lhs.constant != 0
        c = lhs.constant
        for j = 1:m
            push!(ret.aff.vars,   rhs.vars[j])
            push!(ret.aff.coeffs, rhs.coeffs[j] * c)
        end
        ret.aff.constant += c * rhs.constant
    end

    # [RHS constant] * LHS
    if rhs.constant != 0
        c = rhs.constant
        for i = 1:n
            push!(ret.aff.vars,   lhs.vars[i])
            push!(ret.aff.coeffs, lhs.coeffs[i] * c)
        end
        # Don't do the following line
        #ret.aff.constant += c * lhs.constant
        # If lhs.constant is 0, its a waste of time
        # If lhs.constant is non-zero, its already done
    end

    return ret
end
# AffExpr--QuadExpr
(+)(a::AffExpr, q::QuadExpr) = QuadExpr(copy(q.qvars1),copy(q.qvars2),copy(q.qcoeffs),a+q.aff)
(-)(a::AffExpr, q::QuadExpr) = QuadExpr(copy(q.qvars1),copy(q.qvars2),    -q.qcoeffs ,a-q.aff)

# QuadExpr
(+)(lhs::QuadExpr) = lhs
(-)(lhs::QuadExpr) = 0.0-lhs
(*)(lhs::QuadExpr) = lhs
# QuadExpr--Number
(+)(lhs::QuadExpr, rhs::Number) = (+)(+rhs,lhs)
(-)(lhs::QuadExpr, rhs::Number) = (+)(-rhs,lhs)
(*)(lhs::QuadExpr, rhs::Number) = (*)(rhs,lhs)
(/)(lhs::QuadExpr, rhs::Number) = (*)(1.0/rhs,lhs)
# QuadExpr--Variable
(+)(q::QuadExpr, v::Variable) = (+)(v,q)
(-)(q::QuadExpr, v::Variable) = QuadExpr(copy(q.qvars1),copy(q.qvars2),copy(q.qcoeffs),q.aff-v)
(*)(q::QuadExpr, v::Variable) = error("Cannot multiply a quadratic expression by a variable")
(/)(q::QuadExpr, v::Variable) = error("Cannot divide a quadratic expression by a variable")
# QuadExpr--AffExpr
(+)(q::QuadExpr, a::AffExpr) = QuadExpr(copy(q.qvars1),copy(q.qvars2),copy(q.qcoeffs),q.aff+a)
(-)(q::QuadExpr, a::AffExpr) = QuadExpr(copy(q.qvars1),copy(q.qvars2),copy(q.qcoeffs),q.aff-a)
(*)(q::QuadExpr, a::AffExpr) = error("Cannot multiply a quadratic expression by an aff. expression")
(/)(q::QuadExpr, a::AffExpr) = error("Cannot divide a quadratic expression by an aff. expression")
# QuadExpr--QuadExpr
(+)(q1::QuadExpr, q2::QuadExpr) = QuadExpr( vcat(q1.qvars1, q2.qvars1),     vcat(q1.qvars2, q2.qvars2),
                                            vcat(q1.qcoeffs, q2.qcoeffs),   q1.aff + q2.aff)
(-)(q1::QuadExpr, q2::QuadExpr) = QuadExpr( vcat(q1.qvars1, q2.qvars1),     vcat(q1.qvars2, q2.qvars2),
                                            vcat(q1.qcoeffs, -q2.qcoeffs),  q1.aff - q2.aff)

(==)(lhs::AffExpr,rhs::AffExpr) = (lhs.vars == rhs.vars) && (lhs.coeffs == rhs.coeffs) && (lhs.constant == rhs.constant)
(==)(lhs::QuadExpr,rhs::QuadExpr) = (lhs.qvars1 == rhs.qvars1) && (lhs.qvars2 == rhs.qvars2) && (lhs.qcoeffs == rhs.qcoeffs) && (lhs.aff == rhs.aff)

#############################################################################
# Helpers to initialize memory for AffExpr/QuadExpr
#############################################################################

_sizehint_expr!(q::GenericAffExpr, n::Int) = begin
        sizehint!(q.vars,   length(q.vars)   + n)
        sizehint!(q.coeffs, length(q.coeffs) + n)
        nothing
end

_sizehint_expr!(q::GenericQuadExpr, n::Int) = begin
        sizehint!(q.qvars1,  length(q.qvars1)  + n)
        sizehint!(q.qvars2,  length(q.qvars2)  + n)
        sizehint!(q.qcoeffs, length(q.qcoeffs) + n)
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
Base.sum(j::Array{Variable}) = AffExpr(vec(j), ones(length(j)), 0.0)
Base.sum(j::AbstractArray{Variable}) = sum([j[i] for i in eachindex(j)]) # to handle non-one-indexed arrays.
function Base.sum(affs::AbstractArray{T}) where T<:GenericAffExpr
    new_aff = zero(T)
    for aff in affs
        append!(new_aff, aff)
    end
    return new_aff
end

import Base.vecdot

# Base Julia's generic fallback vecdot requires that dot be defined
# for scalars, so instead of defining them one-by-one, we will
# fallback to the multiplication operator
Base.dot(lhs::JuMPTypes, rhs::JuMPTypes) = lhs*rhs
Base.dot(lhs::JuMPTypes, rhs::Number)    = lhs*rhs
Base.dot(lhs::Number,    rhs::JuMPTypes) = lhs*rhs

Base.dot(lhs::AbstractVector{T},rhs::AbstractVector{S}) where {T<:JuMPTypes,S<:JuMPTypes} = _dot(lhs,rhs)
Base.dot(lhs::AbstractVector{T},rhs::AbstractVector{S}) where {T<:JuMPTypes,S} = _dot(lhs,rhs)
Base.dot(lhs::AbstractVector{T},rhs::AbstractVector{S}) where {T,S<:JuMPTypes} = _dot(lhs,rhs)

# TODO: qualify Base.vecdot once v0.3 support is dropped
vecdot(lhs::AbstractArray{T,N},rhs::AbstractArray{S,N}) where {T<:JuMPTypes,S,N} = _dot(lhs,rhs)
vecdot(lhs::AbstractArray{T,N},rhs::AbstractArray{S,N}) where {T<:JuMPTypes,S<:JuMPTypes,N} = _dot(lhs,rhs)
vecdot(lhs::AbstractArray{T,N},rhs::AbstractArray{S,N}) where {T,S<:JuMPTypes,N} = _dot(lhs,rhs)

function _dot(lhs::AbstractArray{T}, rhs::AbstractArray{S}) where {T,S}
    size(lhs) == size(rhs) || error("Incompatible dimensions")
    ret = zero(one(T)*one(S))
    for (x,y) in zip(lhs,rhs)
        ret = addtoexpr(ret, x, y)
    end
    ret
end

###############################################################################
# A bunch of operator junk to make matrix multiplication and friends act
# reasonably sane with JuMP types

Base.promote_rule(::Type{Variable},::Type{R}       ) where {R<:Real} = AffExpr
Base.promote_rule(         ::Type{Variable},::Type{AffExpr} ) = AffExpr
Base.promote_rule(         ::Type{Variable},::Type{QuadExpr}) = QuadExpr
Base.promote_rule(::Type{AffExpr}, ::Type{R}       ) where {R<:Real} = AffExpr
Base.promote_rule(         ::Type{AffExpr}, ::Type{QuadExpr}) = QuadExpr
Base.promote_rule(::Type{QuadExpr},::Type{R}       ) where {R<:Real} = QuadExpr

Base.transpose(x::AbstractJuMPScalar) = x

# Can remove the following code once == overloading is removed

function Base.issymmetric(x::Matrix{T}) where T<:JuMPTypes
    (n = size(x,1)) == size(x,2) || return false
    for i in 1:n, j in (i+1):n
        isequal(x[i,j], x[j,i]) || return false
    end
    true
end

# Special-case because the the base version wants to do fill!(::Array{Variable}, zero(AffExpr))
function Base.diagm(x::AbstractVector{Variable})
    @assert one_indexed(x) # Base.diagm doesn't work for non-one-indexed arrays in general.
    diagm(copy!(similar(x, AffExpr), x))
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
            append!(q, tmp)
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
            append!(q, tmp)
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
                append!(ret[rv[j], k], nzv[j] * rhs[col,k])
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
                append!(q, lhs[multivec_row, rowval[k]] * nzval[k])
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
                append!(q, lhs[rowval[k], multivec_row] * nzval[k]) # transpose
            end
        end
    end
    ret
end


# TODO: implement sparse * sparse code as in base/sparse/linalg.jl (spmatmul)
_multiply!(ret::AbstractArray{T}, lhs::SparseMatrixCSC, rhs::SparseMatrixCSC) where {T<:JuMPTypes} = _multiply!(ret, lhs, full(rhs))
_multiplyt!(ret::AbstractArray{T}, lhs::SparseMatrixCSC, rhs::SparseMatrixCSC) where {T<:JuMPTypes} = _multiplyt!(ret, lhs, full(rhs))

_multiply!(ret, lhs, rhs) = A_mul_B!(ret, lhs, rhs)
_multiplyt!(ret, lhs, rhs) = At_mul_B!(ret, lhs, rhs)

(*)(             A::Union{Matrix{T},SparseMatrixCSC{T}}, x::Union{Matrix,   Vector,   SparseMatrixCSC}) where {T<:JuMPTypes}    = _matmul(A, x)
(*)(A::Union{Matrix{T},SparseMatrixCSC{T}}, x::Union{Matrix{R},Vector{R},SparseMatrixCSC{R}}) where {T<:JuMPTypes,R<:JuMPTypes} = _matmul(A, x)
(*)(             A::Union{Matrix,   SparseMatrixCSC},    x::Union{Matrix{T},Vector{T},SparseMatrixCSC{T}}) where {T<:JuMPTypes} = _matmul(A, x)

import Base.At_mul_B
import Base.Ac_mul_B
# these methods are called when one does A.'*v or A'*v respectively
At_mul_B(A::Union{Matrix{T},SparseMatrixCSC{T}}, x::Union{Matrix, Vector, SparseMatrixCSC}) where {T<:JuMPTypes} = _matmult(A, x)
At_mul_B(A::Union{Matrix{T},SparseMatrixCSC{T}}, x::Union{Matrix{R}, Vector{R}, SparseMatrixCSC{R}}) where {T<:JuMPTypes,R<:JuMPTypes} = _matmult(A, x)
At_mul_B(A::Union{Matrix,SparseMatrixCSC}, x::Union{Matrix{T}, Vector{T}, SparseMatrixCSC{T}}) where {T<:JuMPTypes} = _matmult(A, x)
# these methods are the same as above since complex numbers are not implemented in JuMP
Ac_mul_B(A::Union{Matrix{T},SparseMatrixCSC{T}}, x::Union{Matrix, Vector, SparseMatrixCSC}) where {T<:JuMPTypes} = _matmult(A, x)
Ac_mul_B(A::Union{Matrix{T},SparseMatrixCSC{T}}, x::Union{Matrix{R}, Vector{R}, SparseMatrixCSC{R}}) where {T<:JuMPTypes,R<:JuMPTypes} = _matmult(A, x)
Ac_mul_B(A::Union{Matrix,SparseMatrixCSC}, x::Union{Matrix{T}, Vector{T}, SparseMatrixCSC{T}}) where {T<:JuMPTypes} = _matmult(A, x)

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
_return_array(A::AbstractMatrix{R}, x::AbstractVector{S}) where {R,S} = _fillwithzeros(Array{_matprod_type(R,S)}(size(A,1)))
_return_array(A::AbstractMatrix{R}, x::AbstractMatrix{S}) where {R,S} = _fillwithzeros(Array{_matprod_type(R,S)}(size(A,1), size(x,2)))
# these are for transpose return matrices
_return_arrayt(A::AbstractMatrix{R}, x::AbstractVector{S}) where {R,S} = _fillwithzeros(Array{_matprod_type(R,S)}(size(A,2)))
_return_arrayt(A::AbstractMatrix{R}, x::AbstractMatrix{S}) where {R,S} = _fillwithzeros(Array{_matprod_type(R,S)}(size(A,2), size(x, 2)))

# helper so we don't fill the buffer array with the same object
function _fillwithzeros(arr::AbstractArray{T}) where T
    for I in eachindex(arr)
        arr[I] = zero(T)
    end
    arr
end

# Special-case sparse matrix scalar multiplication/division
(*)(lhs::Number, rhs::SparseMatrixCSC{T}) where {T<:JuMPTypes} =
    SparseMatrixCSC(rhs.m, rhs.n, copy(rhs.colptr), copy(rhs.rowval), lhs .* rhs.nzval)
(*)(lhs::JuMPTypes, rhs::SparseMatrixCSC) =
    SparseMatrixCSC(rhs.m, rhs.n, copy(rhs.colptr), copy(rhs.rowval), lhs .* rhs.nzval)
(*)(lhs::SparseMatrixCSC{T}, rhs::Number) where {T<:JuMPTypes} =
    SparseMatrixCSC(lhs.m, lhs.n, copy(lhs.colptr), copy(lhs.rowval), lhs.nzval .* rhs)
(*)(lhs::SparseMatrixCSC, rhs::JuMPTypes) =
    SparseMatrixCSC(lhs.m, lhs.n, copy(lhs.colptr), copy(lhs.rowval), lhs.nzval .* rhs)
(/)(lhs::SparseMatrixCSC{T}, rhs::Number) where {T<:JuMPTypes} =
    SparseMatrixCSC(lhs.m, lhs.n, copy(lhs.colptr), copy(lhs.rowval), lhs.nzval ./ rhs)


for (op,opsymbol) in [(+,:+), (-,:-), (*,:*), (/,:/)]
    @eval begin
        Base.broadcast(::typeof($op),lhs::Number,rhs::JuMPTypes) = $opsymbol(lhs,rhs)
        Base.broadcast(::typeof($op),lhs::JuMPTypes,rhs::Number) = $opsymbol(lhs,rhs)
    end
end


(+)(x::AbstractArray{T}) where {T<:JuMPTypes} = x
function (-)(x::AbstractArray{T}) where T<:JuMPTypes
    ret = similar(x, typeof(-one(T)))
    for I in eachindex(ret)
        ret[I] = -x[I]
    end
    ret
end
(*)(x::AbstractArray{T}) where {T<:JuMPTypes} = x

###############################################################################
# Add nonlinear function fallbacks for JuMP built-in types
const op_hint = "Are you trying to build a nonlinear problem? Make sure you use @NLconstraint/@NLobjective."
for (func,_) in Calculus.symbolic_derivatives_1arg(), typ in [:Variable,:AffExpr,:QuadExpr]
    errstr = "$func is not defined for type $typ. $op_hint"
    if isdefined(Base, func)
        @eval Base.$(func)(::$typ) = error($errstr)
    end
end

*(::T,::S) where {T<:QuadExpr,S<:Union{Variable,AffExpr,QuadExpr}} =
    error( "*(::$T,::$S) is not defined. $op_hint")
(*)(lhs::QuadExpr, rhs::QuadExpr) =
    error( "*(::QuadExpr,::QuadExpr) is not defined. $op_hint")
*(::S,::T) where {T<:QuadExpr,S<:Union{Variable,AffExpr,QuadExpr}} =
    error( "*(::$S,::$T) is not defined. $op_hint")
/(::S,::T) where {S<:Union{Number,Variable,AffExpr,QuadExpr},T<:Union{Variable,AffExpr,QuadExpr}} =
    error( "/(::$S,::$T) is not defined. $op_hint")

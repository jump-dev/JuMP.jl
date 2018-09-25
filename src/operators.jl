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
# 3. [Generic]Norm
# 4. [Generic]AffExpr
# 5. QuadExpr
# 6. [Generic]NormExpr

# Number
# Number--Number obviously already taken care of!
# Number--Variable
Base.:+(lhs::Number, rhs::Variable) = AffExpr([rhs],[+1.],convert(Float64,lhs))
Base.:-(lhs::Number, rhs::Variable) = AffExpr([rhs],[-1.],convert(Float64,lhs))
Base.:*(lhs::Number, rhs::Variable) = AffExpr([rhs],[convert(Float64,lhs)], 0.)
# Number--GenericNorm
Base.:+(lhs::Number, rhs::GenericNorm{P,C,V}) where {P,C,V} = GenericNormExpr{P,C,V}(copy(rhs),  one(C), convert(C,lhs))
Base.:-(lhs::Number, rhs::GenericNorm{P,C,V}) where {P,C,V} = GenericNormExpr{P,C,V}(copy(rhs), -one(C), convert(C,lhs))
Base.:*(lhs::Number, rhs::GenericNorm{P,C,V}) where {P,C,V} = GenericNormExpr{P,C,V}(copy(rhs),     lhs, zero(GenericAffExpr{C,V}))
# Number--GenericAffExpr
Base.:+(lhs::Number, rhs::GenericAffExpr) = GenericAffExpr(copy(rhs.vars),copy(rhs.coeffs),lhs+rhs.constant)
Base.:-(lhs::Number, rhs::GenericAffExpr) = GenericAffExpr(copy(rhs.vars),    -rhs.coeffs ,lhs-rhs.constant)
Base.:*(lhs::Number, rhs::GenericAffExpr) = GenericAffExpr(copy(rhs.vars),[lhs*rhs.coeffs[i] for i=1:length(rhs.coeffs)],lhs*rhs.constant)
# Number--QuadExpr
Base.:+(lhs::Number, rhs::QuadExpr) = QuadExpr(copy(rhs.qvars1),copy(rhs.qvars2),copy(rhs.qcoeffs),lhs+rhs.aff)
Base.:-(lhs::Number, rhs::QuadExpr) = QuadExpr(copy(rhs.qvars1),copy(rhs.qvars2),    -rhs.qcoeffs ,lhs-rhs.aff)
Base.:*(lhs::Number, rhs::QuadExpr) = QuadExpr(copy(rhs.qvars1),copy(rhs.qvars2), lhs*rhs.qcoeffs ,lhs*rhs.aff)
# Number--GenericNormExpr
Base.:+(lhs::Number, rhs::T) where {T<:GenericNormExpr} = T(copy(rhs.norm),     rhs.coeff, lhs+rhs.aff)
Base.:-(lhs::Number, rhs::T) where {T<:GenericNormExpr} = T(copy(rhs.norm),    -rhs.coeff, lhs-rhs.aff)
Base.:*(lhs::Number, rhs::T) where {T<:GenericNormExpr} = T(copy(rhs.norm), lhs*rhs.coeff, lhs*rhs.aff)

# Variable (or, AbstractJuMPScalar)
Base.:+(lhs::AbstractJuMPScalar) = lhs
Base.:-(lhs::Variable) = AffExpr([lhs],[-1.0],0.0)
Base.:*(lhs::AbstractJuMPScalar) = lhs # make this more generic so extensions don't have to define unary multiplication for our macros
# Variable--Number
Base.:+(lhs::Variable, rhs::Number) = (+)( rhs,lhs)
Base.:-(lhs::Variable, rhs::Number) = (+)(-rhs,lhs)
Base.:*(lhs::Variable, rhs::Number) = (*)(rhs,lhs)
Base.:/(lhs::Variable, rhs::Number) = (*)(1 ./ rhs,lhs)
# Variable--Variable
Base.:+(lhs::Variable, rhs::Variable) = AffExpr([lhs,rhs], [1.,+1.], 0.)
Base.:-(lhs::Variable, rhs::Variable) = AffExpr([lhs,rhs], [1.,-1.], 0.)
Base.:*(lhs::Variable, rhs::Variable) = QuadExpr([lhs],[rhs],[1.],AffExpr(Variable[],Float64[],0.))
# Variable--Norm
Base.:+(lhs::Variable, rhs::GenericNorm{2,C,V}) where {C,V<:Variable} = GenericSOCExpr{C,V}(copy(rhs),  one(C), GenericAffExpr{C,V}(lhs))
Base.:-(lhs::Variable, rhs::GenericNorm{2,C,V}) where {C,V<:Variable} = GenericSOCExpr{C,V}(copy(rhs), -one(C), GenericAffExpr{C,V}(lhs))
# Variable--AffExpr
Base.:+(lhs::V, rhs::GenericAffExpr{C,V}) where {C,V<:JuMPTypes} =
    GenericAffExpr{C,V}(vcat(rhs.vars,lhs),vcat(rhs.coeffs,one(C)), rhs.constant)
Base.:-(lhs::V, rhs::GenericAffExpr{C,V}) where {C,V<:JuMPTypes} =
    GenericAffExpr{C,V}(vcat(rhs.vars,lhs),vcat(-rhs.coeffs,one(C)),-rhs.constant)
function Base.:*(lhs::Variable, rhs::AffExpr)
    n = length(rhs.vars)
    if rhs.constant != 0.
        ret = QuadExpr([lhs for i=1:n],copy(rhs.vars),copy(rhs.coeffs),AffExpr([lhs], [rhs.constant], 0.))
    else
        ret = QuadExpr([lhs for i=1:n],copy(rhs.vars),copy(rhs.coeffs),zero(AffExpr))
    end
end
Base.:/(lhs::Variable, rhs::AffExpr) = error("Cannot divide a variable by an affine expression")
# Variable--QuadExpr
Base.:+(v::Variable, q::QuadExpr) = QuadExpr(copy(q.qvars1),copy(q.qvars2),copy(q.qcoeffs),v+q.aff)
Base.:-(v::Variable, q::QuadExpr) = QuadExpr(copy(q.qvars1),copy(q.qvars2),    -q.qcoeffs ,v-q.aff)
# Variable--GenericSOCExpr
Base.:+(lhs::Variable, rhs::GenericSOCExpr{C,V}) where {C,V<:Variable} = GenericSOCExpr{C,V}(copy(rhs.norm),  rhs.coeff, lhs+rhs.aff)
Base.:-(lhs::Variable, rhs::GenericSOCExpr{C,V}) where {C,V<:Variable} = GenericSOCExpr{C,V}(copy(rhs.norm), -rhs.coeff, lhs-rhs.aff)

# GenericNorm
Base.:+(lhs::GenericNorm) = lhs
Base.:-(lhs::GenericNorm{P,C,V}) where {P,C,V} = GenericNormExpr{P,C,V}(copy(lhs), -one(C), zero(GenericAffExpr{C,V}))
Base.:*(lhs::GenericNorm) = lhs
# GenericNorm--Number
Base.:+(lhs::GenericNorm{P,C,V},rhs::Number) where {P,C,V} =
    GenericNormExpr{P,C,V}(copy(lhs), one(C), GenericAffExpr{C,V}( rhs))
Base.:-(lhs::GenericNorm{P,C,V},rhs::Number) where {P,C,V} =
    GenericNormExpr{P,C,V}(copy(lhs), one(C), GenericAffExpr{C,V}(-rhs))
Base.:*(lhs::GenericNorm{P,C,V},rhs::Number) where {P,C,V} =
    GenericNormExpr{P,C,V}(copy(lhs), rhs, zero(GenericAffExpr{C,V}))
Base.:/(lhs::GenericNorm{P,C,V},rhs::Number) where {P,C,V} =
    GenericNormExpr{P,C,V}(copy(lhs), 1/rhs, zero(GenericAffExpr{C,V}))
# Norm--Variable
Base.:+(lhs::Norm,rhs::Variable) = SOCExpr(copy(lhs), 1.0, AffExpr( rhs))
Base.:-(lhs::Norm,rhs::Variable) = SOCExpr(copy(lhs), 1.0, AffExpr(-rhs))
# GenericNorm--GenericNorm
# GenericNorm--GenericAffExpr
Base.:+(lhs::GenericNorm{P,C,V},rhs::GenericAffExpr{C,V}) where {P,C,V} =
    GenericNormExpr{P,C,V}(copy(lhs), one(C), copy(rhs))
Base.:-(lhs::GenericNorm{P,C,V},rhs::GenericAffExpr{C,V}) where {P,C,V} =
    GenericNormExpr{P,C,V}(copy(lhs), one(C),     -rhs)
# Norm--QuadExpr
# Norm--SOCExpr

# GenericAffExpr
Base.:+(lhs::GenericAffExpr) = lhs
Base.:-(lhs::GenericAffExpr) = GenericAffExpr(lhs.vars, -lhs.coeffs, -lhs.constant)
Base.:*(lhs::GenericAffExpr) = lhs
# AffExpr--Number
Base.:+(lhs::GenericAffExpr, rhs::Number) = (+)(+rhs,lhs)
Base.:-(lhs::GenericAffExpr, rhs::Number) = (+)(-rhs,lhs)
Base.:*(lhs::GenericAffExpr, rhs::Number) = (*)(rhs,lhs)
Base.:/(lhs::GenericAffExpr, rhs::Number) = (*)(1.0/rhs,lhs)
function Base.:^(lhs::Union{Variable,AffExpr}, rhs::Integer)
    if rhs == 2
        return lhs*lhs
    elseif rhs == 1
        return convert(QuadExpr, lhs)
    elseif rhs == 0
        return convert(QuadExpr, 1)
    else
        error("Only exponents of 0, 1, or 2 are currently supported. Are you trying to build a nonlinear problem? Make sure you use @NLconstraint/@NLobjective.")
    end
end
Base.:^(lhs::Union{Variable,AffExpr}, rhs::Number) = error("Only exponents of 0, 1, or 2 are currently supported. Are you trying to build a nonlinear problem? Make sure you use @NLconstraint/@NLobjective.")
# AffExpr--Variable
Base.:+(lhs::GenericAffExpr{C,V}, rhs::V) where {C,V<:JuMPTypes} = (+)(rhs,lhs)
Base.:-(lhs::GenericAffExpr{C,V}, rhs::V) where {C,V<:JuMPTypes} = GenericAffExpr{C,V}(vcat(lhs.vars,rhs),vcat(lhs.coeffs,-one(C)),lhs.constant)
Base.:*(lhs::AffExpr, rhs::Variable) = (*)(rhs,lhs)
Base.:/(lhs::AffExpr, rhs::Variable) = error("Cannot divide affine expression by a variable")
# AffExpr--Norm
Base.:+(lhs::GenericAffExpr{C,V}, rhs::GenericNorm{P,C,V}) where {P,C,V} = GenericNormExpr{P,C,V}(copy(rhs),  one(C), copy(lhs))
Base.:-(lhs::GenericAffExpr{C,V}, rhs::GenericNorm{P,C,V}) where {P,C,V} = GenericNormExpr{P,C,V}(copy(rhs), -one(C), copy(lhs))
# AffExpr--AffExpr
Base.:+(lhs::GenericAffExpr{C,V}, rhs::GenericAffExpr{C,V}) where {C,V<:JuMPTypes} = (operator_warn(lhs,rhs); GenericAffExpr(vcat(lhs.vars,rhs.vars),vcat(lhs.coeffs, rhs.coeffs),lhs.constant+rhs.constant))
Base.:-(lhs::GenericAffExpr{C,V}, rhs::GenericAffExpr{C,V}) where {C,V<:JuMPTypes} = GenericAffExpr(vcat(lhs.vars,rhs.vars),vcat(lhs.coeffs,-rhs.coeffs),lhs.constant-rhs.constant)
function Base.:*(lhs::AffExpr, rhs::AffExpr)
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
Base.:+(a::AffExpr, q::QuadExpr) = QuadExpr(copy(q.qvars1),copy(q.qvars2),copy(q.qcoeffs),a+q.aff)
Base.:-(a::AffExpr, q::QuadExpr) = QuadExpr(copy(q.qvars1),copy(q.qvars2),    -q.qcoeffs ,a-q.aff)
# AffExpr--SOCExpr
Base.:+(lhs::GenericAffExpr{C,V}, rhs::GenericSOCExpr{C,V}) where {C,V} = GenericSOCExpr{C,V}(copy(rhs.norm),  copy(rhs.coeff), lhs+rhs.aff)
Base.:-(lhs::GenericAffExpr{C,V}, rhs::GenericSOCExpr{C,V}) where {C,V} = GenericSOCExpr{C,V}(copy(rhs.norm), -copy(rhs.coeff), lhs-rhs.aff)

# QuadExpr
Base.:+(lhs::QuadExpr) = lhs
Base.:-(lhs::QuadExpr) = 0.0-lhs
Base.:*(lhs::QuadExpr) = lhs
# QuadExpr--Number
Base.:+(lhs::QuadExpr, rhs::Number) = (+)(+rhs,lhs)
Base.:-(lhs::QuadExpr, rhs::Number) = (+)(-rhs,lhs)
Base.:*(lhs::QuadExpr, rhs::Number) = (*)(rhs,lhs)
Base.:/(lhs::QuadExpr, rhs::Number) = (*)(1.0/rhs,lhs)
# QuadExpr--Variable
Base.:+(q::QuadExpr, v::Variable) = (+)(v,q)
Base.:-(q::QuadExpr, v::Variable) = QuadExpr(copy(q.qvars1),copy(q.qvars2),copy(q.qcoeffs),q.aff-v)
Base.:*(q::QuadExpr, v::Variable) = error("Cannot multiply a quadratic expression by a variable")
Base.:/(q::QuadExpr, v::Variable) = error("Cannot divide a quadratic expression by a variable")
# QuadExpr--AffExpr
Base.:+(q::QuadExpr, a::AffExpr) = QuadExpr(copy(q.qvars1),copy(q.qvars2),copy(q.qcoeffs),q.aff+a)
Base.:-(q::QuadExpr, a::AffExpr) = QuadExpr(copy(q.qvars1),copy(q.qvars2),copy(q.qcoeffs),q.aff-a)
Base.:*(q::QuadExpr, a::AffExpr) = error("Cannot multiply a quadratic expression by an aff. expression")
Base.:/(q::QuadExpr, a::AffExpr) = error("Cannot divide a quadratic expression by an aff. expression")
# QuadExpr--QuadExpr
Base.:+(q1::QuadExpr, q2::QuadExpr) = QuadExpr( vcat(q1.qvars1, q2.qvars1),     vcat(q1.qvars2, q2.qvars2),
                                            vcat(q1.qcoeffs, q2.qcoeffs),   q1.aff + q2.aff)
Base.:-(q1::QuadExpr, q2::QuadExpr) = QuadExpr( vcat(q1.qvars1, q2.qvars1),     vcat(q1.qvars2, q2.qvars2),
                                            vcat(q1.qcoeffs, -q2.qcoeffs),  q1.aff - q2.aff)

# GenericNormExpr
Base.:+(lhs::GenericNormExpr) = lhs
Base.:-(lhs::T) where {T<:GenericNormExpr} = T(copy(lhs.norm), -lhs.coeff, -lhs.aff)
Base.:*(lhs::GenericNormExpr) = lhs
# GenericNormExpr--Number
Base.:+(lhs::T,rhs::Number) where {T<:GenericNormExpr} = T(copy(lhs.norm), lhs.coeff, lhs.aff+rhs)
Base.:-(lhs::T,rhs::Number) where {T<:GenericNormExpr} = T(copy(lhs.norm), lhs.coeff, lhs.aff-rhs)
Base.:*(lhs::T,rhs::Number) where {T<:GenericNormExpr} = T(copy(lhs.norm), lhs.coeff*rhs, lhs.aff*rhs)
Base.:/(lhs::T,rhs::Number) where {T<:GenericNormExpr} = T(copy(lhs.norm), lhs.coeff/rhs, lhs.aff/rhs)
# SOCExpr--Variable
Base.:+(lhs::SOCExpr,rhs::Variable) = SOCExpr(copy(lhs.norm), lhs.coeff, lhs.aff+rhs)
Base.:-(lhs::SOCExpr,rhs::Variable) = SOCExpr(copy(lhs.norm), lhs.coeff, lhs.aff-rhs)
# GenericNormExpr--GenericNorm
# GenericNormExpr--GenericAffExpr
Base.:+(lhs::GenericNormExpr{P,C,V},rhs::GenericAffExpr{C,V}) where {P,C,V} =
    GenericNormExpr{P,C,V}(copy(lhs.norm), lhs.coeff, lhs.aff+rhs)
Base.:-(lhs::GenericNormExpr{P,C,V},rhs::GenericAffExpr{C,V}) where {P,C,V} =
    GenericNormExpr{P,C,V}(copy(lhs.norm), lhs.coeff, lhs.aff-rhs)
# SOCExpr--QuadExpr
# SOCExpr--SOCExpr

Base.:(==)(lhs::AffExpr,rhs::AffExpr) = (lhs.vars == rhs.vars) && (lhs.coeffs == rhs.coeffs) && (lhs.constant == rhs.constant)
Base.:(==)(lhs::QuadExpr,rhs::QuadExpr) = (lhs.qvars1 == rhs.qvars1) && (lhs.qvars2 == rhs.qvars2) && (lhs.qcoeffs == rhs.qcoeffs) && (lhs.aff == rhs.aff)

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

Base.sum(j::JuMPArray) = sum(j.innerArray)
Base.sum(j::JuMPDict)  = sum(values(j.tupledict))
Base.sum(j::JuMPArray{Variable}) = AffExpr(vec(j.innerArray), ones(length(j.innerArray)), 0.0)
Base.sum(j::JuMPDict{Variable})  = AffExpr(collect(values(j.tupledict)), ones(length(j.tupledict)), 0.0)
Base.sum(j::Array{Variable}) = AffExpr(vec(j), ones(length(j)), 0.0)
Base.sum(j::AbstractArray{Variable}) = sum([j[i] for i in eachindex(j)]) # to handle non-one-indexed arrays.
function Base.sum(affs::AbstractArray{T}) where T<:GenericAffExpr
    new_aff = zero(T)
    for aff in affs
        append!(new_aff, aff)
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

Compat.LinearAlgebra.dot(lhs::AbstractArray{T,N},rhs::AbstractArray{S,N}) where {T<:JuMPTypes,S<:JuMPTypes,N} = _dot(lhs,rhs)
Compat.LinearAlgebra.dot(lhs::AbstractArray{T,N},rhs::AbstractArray{S,N}) where {T<:JuMPTypes,S,N} = _dot(lhs,rhs)
Compat.LinearAlgebra.dot(lhs::AbstractArray{T,N},rhs::AbstractArray{S,N}) where {T,S<:JuMPTypes,N} = _dot(lhs,rhs)

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

_throw_transpose_error() = error("Transpose not currently implemented for JuMPArrays with arbitrary index sets.")
Compat.LinearAlgebra.transpose(x::AbstractJuMPScalar) = x
Compat.LinearAlgebra.transpose( x::JuMPArray) = _throw_transpose_error()
Compat.adjoint(x::AbstractJuMPScalar) = x
Compat.adjoint(x::JuMPArray) = _throw_transpose_error()

# Can remove the following code once == overloading is removed

function Compat.LinearAlgebra.issymmetric(x::Matrix{T}) where T<:JuMPTypes
    (n = size(x,1)) == size(x,2) || return false
    for i in 1:n, j in (i+1):n
        isequal(x[i,j], x[j,i]) || return false
    end
    true
end

# Special-case because the the base version wants to do fill!(::Array{Variable}, zero(AffExpr))
function Compat.LinearAlgebra.diagm(x::AbstractVector{Variable})
    @assert one_indexed(x) # Base.diagm doesn't work for non-one-indexed arrays in general.
    diagm(copyto!(similar(x, AffExpr), x))
end

###############
# The _multiply!(buf,y,z) adds the results of y*z into the buffer buf. No bounds/size
# checks are performed; it is expected that the caller has done this, has ensured
# that the eltype of buf is appropriate, and has zeroed the elements of buf (if desired).

# this is a generic fallback
function _multiply!(ret::Array{T}, lhs::Array, rhs::Array) where T
    #@warn("A terrible fallback is being called!")  # occasionally check for this in testing
    m, n = size(lhs,1), size(lhs,2)
    r, s = size(rhs,1), size(rhs,2)
    for i ∈ 1:m, j ∈ 1:s
        for k ∈ 1:n
            tmp = convert(T, lhs[i,k]*rhs[k,j])
            ret[i,j] += tmp
        end
    end
    ret
end


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

# See https://github.com/JuliaLang/julia/issues/27015
function Base.Matrix(S::SparseMatrixCSC{V}) where V<:Variable
    A = zeros(GenericAffExpr{Float64, V}, S.m, S.n)
    for Sj ∈ 1:S.n
        for Sk ∈ nzrange(S, Sj)
            Si = S.rowval[Sk]
            Sv = S.nzval[Sk]
            A[Si, Sj] = Sv
        end
    end
    A
end


# TODO: implement sparse * sparse code as in base/sparse/linalg.jl (spmatmul)
_multiply!(ret::AbstractArray{T}, lhs::SparseMatrixCSC, rhs::SparseMatrixCSC) where {T<:JuMPTypes} = _multiply!(ret, lhs, Matrix(rhs))
_multiplyt!(ret::AbstractArray{T}, lhs::SparseMatrixCSC, rhs::SparseMatrixCSC) where {T<:JuMPTypes} = _multiplyt!(ret, lhs, Matrix(rhs))

function Base.:*(A::Union{Matrix{T},SparseMatrixCSC{T}},
                 x::Union{Matrix,Vector,SparseMatrixCSC}) where {T<:JuMPTypes}
    _matmul(A, x)
end
function Base.:*(A::Union{Matrix{T},SparseMatrixCSC{T}},
                 x::Union{Matrix{R},Vector{R},SparseMatrixCSC{R}}) where {T<:JuMPTypes,R<:JuMPTypes}
    _matmul(A, x)
end
function Base.:*(A::Union{Matrix,SparseMatrixCSC},
                 x::Union{Matrix{T},Vector{T},SparseMatrixCSC{T}}) where {T<:JuMPTypes}
    _matmul(A, x)
end


for op in [:+, :-]; @eval begin
    function Base.$op(lhs::Number,rhs::AbstractArray{T}) where T<:JuMPTypes
        ret = similar(rhs, typeof($op(lhs, zero(T))))
        for I in eachindex(ret)
            ret[I] = $op(lhs, rhs[I])
        end
        ret
    end
    function Base.$op(lhs::AbstractArray{T},rhs::Number) where T<:JuMPTypes
        ret = similar(lhs, typeof($op(zero(T), rhs)))
        for I in eachindex(ret)
            ret[I] = $op(lhs[I], rhs)
        end
        ret
    end
    function Base.$op(lhs::T,rhs::AbstractArray{S}) where {T<:JuMPTypes,S}
        ret = similar(rhs, typeof($op(lhs, zero(S))))
        for I in eachindex(ret)
            ret[I] = $op(lhs, rhs[I])
        end
        ret
    end
    function Base.$op(lhs::AbstractArray{S},rhs::T) where {T<:JuMPTypes,S}
        ret = similar(lhs, typeof($op(zero(S), rhs)))
        for I in eachindex(ret)
            ret[I] = $op(lhs[I], rhs)
        end
        ret
    end
end; end

for op in [:*, :/]; @eval begin
    function Base.$op(lhs::Number,rhs::AbstractArray{T}) where T<:JuMPTypes
        ret = similar(rhs, typeof($op(lhs, zero(T))))
        for I in eachindex(ret)
            ret[I] = $op(lhs, rhs[I])
        end
        ret
    end
    function Base.$op(lhs::AbstractArray{T},rhs::Number) where T<:JuMPTypes
        ret = similar(lhs, typeof($op(zero(T), rhs)))
        for I in eachindex(ret)
            ret[I] = $op(lhs[I], rhs)
        end
        ret
    end
    function Base.$op(lhs::T,rhs::AbstractArray{S}) where {T<:JuMPTypes,S}
        ret = similar(rhs, typeof($op(lhs, zero(S))))
        for I in eachindex(ret)
            ret[I] = $op(lhs, rhs[I])
        end
        ret
    end
    function Base.$op(lhs::AbstractArray{S},rhs::T) where {T<:JuMPTypes,S}
        ret = similar(lhs, typeof($op(zero(S), rhs)))
        for I in eachindex(ret)
            ret[I] = $op(lhs[I], rhs)
        end
        ret
    end
end; end

if VERSION ≥ v"0.7-"
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
    # WARNING these are inefficient fallbacks
    function Base.:*(adjA::Adjoint{<:JuMPTypes,<:SparseMatrixCSC},
                     x::SparseMatrixCSC)
        adjA * Matrix(x)
    end
    function Base.:*(adjA::Transpose{<:JuMPTypes,<:SparseMatrixCSC},
                     x::SparseMatrixCSC)
        adjA * Matrix(x)
    end
    #function Base.:*(adjA::Adjoint{<:JuMPTypes,<:SparseMatrixCSC},
    #                 x::SparseMatrixCSC)
    #    # WARNING this is an inefficient hack
    #    adjA * Array(x)
    #end
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
_return_array(A::AbstractMatrix{R}, x::AbstractVector{S}) where {R,S} = _fillwithzeros(Array{_matprod_type(R,S)}(undef, size(A,1)))
_return_array(A::AbstractMatrix{R}, x::AbstractMatrix{S}) where {R,S} = _fillwithzeros(Array{_matprod_type(R,S)}(undef, size(A,1), size(x,2)))
# these are for transpose return matrices
_return_arrayt(A::AbstractMatrix{R}, x::AbstractVector{S}) where {R,S} = _fillwithzeros(Array{_matprod_type(R,S)}(undef, size(A,2)))
_return_arrayt(A::AbstractMatrix{R}, x::AbstractMatrix{S}) where {R,S} = _fillwithzeros(Array{_matprod_type(R,S)}(undef, size(A,2), size(x, 2)))

# helper so we don't fill the buffer array with the same object
function _fillwithzeros(arr::AbstractArray{T}) where T
    for I in eachindex(arr)
        arr[I] = zero(T)
    end
    arr
end

# Special-case sparse matrix scalar multiplication/division
function Base.:*(lhs::Number, rhs::SparseMatrixCSC{T}) where {T<:JuMPTypes}
    return map(x -> lhs * x, rhs)
end
function Base.:*(lhs::JuMPTypes, rhs::SparseMatrixCSC)
    return map(x -> lhs * x, rhs)
end
function Base.:*(lhs::SparseMatrixCSC{T}, rhs::Number) where {T<:JuMPTypes}
    return map(x -> x * rhs, lhs)
end
function Base.:*(lhs::SparseMatrixCSC, rhs::JuMPTypes)
    return map(x -> x * rhs, lhs)
end
function Base.:/(lhs::SparseMatrixCSC{T}, rhs::Number) where {T<:JuMPTypes}
    return map(x -> x / rhs, lhs)
end

if VERSION ≥ v"0.7-"
    Base.BroadcastStyle(::Type{<:JuMPTypes}) = Base.Broadcast.DefaultArrayStyle{0}()
    Base.broadcastable(x::JuMPTypes) = fill(x, ())
else
    for (op,opsymbol) in [(+,:+), (-,:-), (*,:*), (/,:/)]
        @eval begin
            Base.broadcast(::typeof($op),lhs::Number,rhs::JuMPTypes) = $opsymbol(lhs,rhs)
            Base.broadcast(::typeof($op),lhs::JuMPTypes,rhs::Number) = $opsymbol(lhs,rhs)
        end
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
for (func,_) in Calculus.symbolic_derivatives_1arg(), typ in [:Variable,:AffExpr,:QuadExpr]
    errstr = "$func is not defined for type $typ. $op_hint"
    if isdefined(Base, func)
        @eval Base.$(func)(::$typ) = error($errstr)
    end
end

Base.:*(::T,::S) where {T<:QuadExpr,S<:Union{Variable,AffExpr,QuadExpr}} =
    error( "*(::$T,::$S) is not defined. $op_hint")
Base.:*(lhs::QuadExpr, rhs::QuadExpr) =
    error( "*(::QuadExpr,::QuadExpr) is not defined. $op_hint")
Base.:*(::S,::T) where {T<:QuadExpr,S<:Union{Variable,AffExpr,QuadExpr}} =
    error( "*(::$S,::$T) is not defined. $op_hint")
Base.:/(::S,::T) where {S<:Union{Number,Variable,AffExpr,QuadExpr},T<:Union{Variable,AffExpr,QuadExpr}} =
    error( "/(::$S,::$T) is not defined. $op_hint")

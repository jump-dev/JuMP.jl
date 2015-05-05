#  Copyright 2015, Iain Dunning, Joey Huchette, Miles Lubin, and contributors
#  This Source Code Form is subject to the terms of the Mozilla Public
#  License, v. 2.0. If a copy of the MPL was not distributed with this
#  file, You can obtain one at http://mozilla.org/MPL/2.0/.
#############################################################################
# JuMP
# An algebraic modelling langauge for Julia
# See http://github.com/JuliaOpt/JuMP.jl
#############################################################################

# Overloads
#
# Different objects that must all interact:
# 1. Number
# 2. Variable
# 3. AffExpr
# 4. QuadExpr

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

# This is not well defined if variable types are different, but needed to avoid ambiguities
(+)(lhs::GenericAffExpr, rhs::GenericAffExpr) = (+)(promote(lhs,rhs)...)
(-)(lhs::GenericAffExpr, rhs::GenericAffExpr) = (-)(promote(lhs,rhs)...)

# Variable
(+)(lhs::Variable) = lhs
(-)(lhs::Variable) = AffExpr([lhs],[-1.0],0.0)
(*)(lhs::Variable) = lhs
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
(+){CoefType,VarType}(lhs::VarType, rhs::GenericAffExpr{CoefType,VarType}) =
    GenericAffExpr{CoefType,VarType}(vcat(rhs.vars,lhs),vcat( rhs.coeffs,one(CoefType)), rhs.constant)
(-){CoefType,VarType}(lhs::VarType, rhs::GenericAffExpr{CoefType,VarType}) =
    GenericAffExpr{CoefType,VarType}(vcat(rhs.vars,lhs),vcat(-rhs.coeffs,one(CoefType)),-rhs.constant)
function (*)(lhs::Variable, rhs::AffExpr)
    n = length(rhs.vars)
    if rhs.constant != 0.
        ret = QuadExpr([lhs for i=1:n],copy(rhs.vars),copy(rhs.coeffs),AffExpr([lhs], [rhs.constant], 0.))
    else
        ret = QuadExpr([lhs for i=1:n],copy(rhs.vars),copy(rhs.coeffs),AffExpr())
    end
end
(/)(lhs::Variable, rhs::AffExpr) = error("Cannot divide a variable by an affine expression")
# Variable--QuadExpr
(+)(v::Variable, q::QuadExpr) = QuadExpr(copy(q.qvars1),copy(q.qvars2),copy(q.qcoeffs),v+q.aff)
(-)(v::Variable, q::QuadExpr) = QuadExpr(copy(q.qvars1),copy(q.qvars2),    -q.qcoeffs ,v-q.aff)

# AffExpr (GenericAffExpr)
(+)(lhs::GenericAffExpr) = lhs
(-)(lhs::GenericAffExpr) = GenericAffExpr(lhs.vars, -lhs.coeffs, -lhs.constant)
(*)(lhs::GenericAffExpr) = lhs
# AffExpr--Number
(+)(lhs::GenericAffExpr, rhs::Number) = (+)(+rhs,lhs)
(-)(lhs::GenericAffExpr, rhs::Number) = (+)(-rhs,lhs)
(*)(lhs::GenericAffExpr, rhs::Number) = (*)(rhs,lhs)
(/)(lhs::GenericAffExpr, rhs::Number) = (*)(1.0/rhs,lhs)
function (^)(lhs::Union(Variable,AffExpr), rhs::Integer)
    rhs == 2 || error("Only exponents of 2 are currently supported. Are you trying to build a nonlinear problem? Make sure you use @addNLConstraint/@setNLObjective.")
    return lhs*lhs
end
(^)(lhs::Union(Variable,AffExpr), rhs::Number) = error("Only exponents of 2 are currently supported. Are you trying to build a nonlinear problem? Make sure you use @addNLConstraint/@setNLObjective.")
# AffExpr--Variable
(+){CoefType,VarType}(lhs::GenericAffExpr{CoefType,VarType}, rhs::VarType) = (+)(rhs,lhs)
(-){CoefType,VarType}(lhs::GenericAffExpr{CoefType,VarType}, rhs::VarType) = GenericAffExpr{CoefType,VarType}(vcat(lhs.vars,rhs),vcat(lhs.coeffs,-one(CoefType)),lhs.constant)
(*)(lhs::AffExpr, rhs::Variable) = (*)(rhs,lhs)
(/)(lhs::AffExpr, rhs::Variable) = error("Cannot divide affine expression by a variable")
# AffExpr--AffExpr
(+){T<:GenericAffExpr}(lhs::T, rhs::T) = GenericAffExpr(vcat(lhs.vars,rhs.vars),vcat(lhs.coeffs, rhs.coeffs),lhs.constant+rhs.constant)
(-){T<:GenericAffExpr}(lhs::T, rhs::T) = GenericAffExpr(vcat(lhs.vars,rhs.vars),vcat(lhs.coeffs,-rhs.coeffs),lhs.constant-rhs.constant)
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

_deprecate_comparisons(sgn) =
    Base.warn_once("The comparison operator $sgn has been deprecated for constructing constraints. Use the macro form @addConstraint instead.")

# LinearConstraint
# Number--???
for (sgn, osgn) in ( (:<=,:>=), (:(==),:(==)), (:>=,:<=) )
    for typ in (:Variable, :AffExpr, :QuadExpr)
        @eval $(sgn)(lhs::Number, rhs::$(typ)) = $(osgn)(rhs, lhs)
    end
    # Variable--???
    for typ in (:Number, :Variable, :AffExpr, :QuadExpr)
        @eval $(sgn)(lhs::Variable, rhs::$(typ)) = $(sgn)(lhs-rhs, 0.0)
    end
end
# AffExpr--???
function (<=)(lhs::AffExpr, rhs::Number)
    _deprecate_comparisons(:(<=))
    LinearConstraint(lhs,-Inf,rhs-lhs.constant)
end
function (==)(lhs::AffExpr, rhs::Number)
    _deprecate_comparisons(:(==))
    LinearConstraint(lhs,rhs-lhs.constant,rhs-lhs.constant)
end
function (>=)(lhs::AffExpr, rhs::Number)
    _deprecate_comparisons(:(>=))
    LinearConstraint(lhs,rhs-lhs.constant,Inf)
end
for sgn in (:<=, :(==), :>=)
    for typ in (:Variable, :AffExpr, :QuadExpr)
        @eval $(sgn)(lhs::AffExpr, rhs::$(typ)) = $(sgn)(lhs-rhs, 0.0)
    end
end
# There's no easy way to allow operator overloads for range constraints.
# Use macros instead.

# QuadConstraint
# QuadConstraint--Number
for sgn in (:<=, :(==), :>=)
    @eval begin
        function $(sgn)(lhs::QuadExpr, rhs::Number)
            _deprecate_comparisons($sgn)
            QuadConstraint( QuadExpr(copy(lhs.qvars1), copy(lhs.qvars2), lhs.qcoeffs,lhs.aff - rhs), $(quot(sgn)))
        end
    end
    for typ in (:Variable, :AffExpr, :QuadExpr)
        @eval $(sgn)(lhs::QuadExpr, rhs::$(typ)) = $(sgn)(lhs-rhs, 0)
    end
end

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
function Base.sum{S,T}(affs::Array{GenericAffExpr{S,T}})
    new_aff = GenericAffExpr{S,T}()
    for aff in affs
        append!(new_aff, aff)
    end
    return new_aff
end

function Base.dot{T,S,N}(lhs::Array{T,N}, rhs::JuMPArray{S,N})
    size(lhs) == size(rhs.innerArray) || error("Incompatible dimensions")
    dot(lhs,rhs.innerArray)
end
Base.dot{S,T,N}(lhs::JuMPArray{S,N},rhs::Array{T,N}) = dot(rhs,lhs)

Base.dot{T<:Real,N}(lhs::Array{T,N}, rhs::JuMPArray{Float64,N}) = dot(vec(lhs), vec(rhs.innerArray))
Base.dot{T<:Real,N}(lhs::JuMPArray{Float64,N}, rhs::Array{T,N}) = dot(vec(rhs), vec(lhs.innerArray))
Base.dot{T<:Real,N}(lhs::Array{T,N}, rhs::Array{Variable,N})   = AffExpr(vec(rhs), vec(float(lhs)), 0.0)
Base.dot{T<:Real,N}(rhs::Array{Variable,N}, lhs::Array{T,N})   = AffExpr(vec(rhs), vec(float(lhs)), 0.0)

function Base.dot{N}(lhs::JuMPArray{Variable,N},rhs::JuMPArray{Variable,N})
    size(lhs.innerArray) == size(rhs.innerArray) || error("Incompatible dimensions")
    return QuadExpr(vec(lhs.innerArray), vec(rhs.innerArray), ones(length(lhs.innerArray)), AffExpr())
end

function Base.dot{N}(lhs::JuMPArray{Float64,N},rhs::JuMPArray{Float64,N})
    size(lhs.innerArray) == size(rhs.innerArray) || error("Incompatible dimensions")
    return sum(lhs.innerArray .* rhs.innerArray)
end

###############################################################################
# Add nonlinear function fallbacks for JuMP built-in types
const op_hint = "Are you trying to build a nonlinear problem? Make sure you use @addNLConstraint/@setNLObjective."
for (func,_) in Calculus.symbolic_derivatives_1arg(), typ in [:Variable,:AffExpr,:QuadExpr]
    errstr = "$func is not defined for type $typ. $op_hint"
    @eval Base.($(quot(func)))(::$typ) = error($errstr)
end

*{T<:QuadExpr,S<:Union(Variable,AffExpr,QuadExpr)}(::T,::S) =
    error( "*(::$T,::$S) is not defined. $op_hint")
(*)(lhs::QuadExpr, rhs::QuadExpr) =
    error( "*(::QuadExpr,::QuadExpr) is not defined. $op_hint")
*{T<:QuadExpr,S<:Union(Variable,AffExpr,QuadExpr)}(::S,::T) =
    error( "*(::$S,::$T) is not defined. $op_hint")
/{S<:Union(Number,Variable,AffExpr,QuadExpr),T<:Union(Variable,AffExpr,QuadExpr)}(::S,::T) =
    error( "/(::$S,::$T) is not defined. $op_hint")

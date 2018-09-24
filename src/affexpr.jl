#  Copyright 2017, Iain Dunning, Joey Huchette, Miles Lubin, and contributors
#  This Source Code Form is subject to the terms of the Mozilla Public
#  License, v. 2.0. If a copy of the MPL was not distributed with this
#  file, You can obtain one at http://mozilla.org/MPL/2.0/.
#############################################################################
# JuMP
# An algebraic modeling language for Julia
# See http://github.com/JuliaOpt/JuMP.jl
#############################################################################
# src/affexpr.jl
# Defines all types relating to affine expressions
# - GenericAffExpr              ∑ aᵢ xᵢ  +  c
#   - AffExpr                   Alias for (Float64, Variable)
# - GenericRangeConstraint      l ≤ ∑ aᵢ xᵢ ≤ u
#   - LinearConstraint          Alias for AffExpr
# Operator overloads in src/operators.jl
#############################################################################

#############################################################################
# GenericAffExpr
# ∑ aᵢ xᵢ  +  c
mutable struct GenericAffExpr{CoefType,VarType} <: AbstractJuMPScalar
    vars::Vector{VarType}
    coeffs::Vector{CoefType}
    constant::CoefType
end
GenericAffExpr(a::GenericAffExpr) = GenericAffExpr(a.vars, a.coeffs, a.constant)

coeftype(::GenericAffExpr{C,V}) where {C,V} = C

# variables must be ∈ ℝ but coeffs can be ∈ ℂ
Compat.adjoint(a::GenericAffExpr{<:Real}) = a
function Compat.adjoint(a::GenericAffExpr)
    GenericAffExpr(a.vars, adjoint.(a.coeffs), a.constant')
end

Base.zero(::Type{GenericAffExpr{C,V}}) where {C,V} = GenericAffExpr{C,V}(V[],C[],zero(C))
Base.one(::Type{GenericAffExpr{C,V}}) where { C,V} = GenericAffExpr{C,V}(V[],C[], one(C))
Base.zero(a::GenericAffExpr) = zero(typeof(a))
Base.one( a::GenericAffExpr) =  one(typeof(a))
Base.copy(a::GenericAffExpr) = GenericAffExpr(copy(a.vars),copy(a.coeffs),copy(a.constant))

Base.convert(::Type{GenericAffExpr{C,V}}, a::GenericAffExpr{C,V}) where {C,V} = a

# Old iterator protocol - iterates over tuples (aᵢ,xᵢ)
struct LinearTermIterator{GAE<:GenericAffExpr}
    aff::GAE
end

linearterms(aff::GenericAffExpr) = LinearTermIterator(aff)

if VERSION < v"0.7-"
    Base.start(lti::LinearTermIterator) = 1
    Base.done( lti::LinearTermIterator, state::Int) = state > length(lti.aff.vars)
    Base.next( lti::LinearTermIterator, state::Int) = ((lti.aff.coeffs[state], lti.aff.vars[state]), state+1)
else
    function Base.iterate(lti::LinearTermIterator)
        if length(lti.aff.vars) ≥ 1
            ((lti.aff.coeffs[1], lti.aff.vars[1]), 2)
        else
            nothing
        end
    end
    function Base.iterate(lti::LinearTermIterator, state::Int)
        if state ≤ length(lti.aff.vars)
            ((lti.aff.coeffs[state], lti.aff.vars[state]), state+1)
        else
            nothing
        end
    end
end

# More efficient ways to grow an affine expression
# Add a single term to an affine expression
function Base.push!(aff::GenericAffExpr{C,V}, new_coeff::C, new_var::V) where {C,V}
    push!(aff.coeffs, new_coeff)
    push!(aff.vars, new_var)
    aff
end
# Add an affine expression to an existing affine expression
function Base.append!(aff::GenericAffExpr{C,V}, other::GenericAffExpr{C,V}) where {C,V}
    append!(aff.vars, other.vars)
    append!(aff.coeffs, other.coeffs)
    aff.constant += other.constant
    aff
end
# For consistency, allow appending constants and individual variables
Base.append!(aff::GenericAffExpr{C,C}, other::C) where {C} = error() # for ambiguity
function Base.append!(aff::GenericAffExpr{C,V}, other::C) where {C,V}
    aff.constant += other
    aff
end
function Base.append!(aff::GenericAffExpr{C,V}, other::Real) where {C,V}
    aff.constant += other
    aff
end
Base.append!(aff::GenericAffExpr{C,V}, other::V) where {C,V} = push!(aff,one(C),other)

function Base.isequal(aff::GenericAffExpr{C,V},other::GenericAffExpr{C,V}) where {C,V}
    isequal(aff.constant, other.constant)  || return false
    length(aff.vars) == length(other.vars) || return false
    for i in 1:length(aff.vars)
        isequal(aff.vars[i],   other.vars[i])   || return false
        isequal(aff.coeffs[i], other.coeffs[i]) || return false
    end
    return true
end


# Alias for (Float64, Variable), the specific GenericAffExpr used by JuMP
const AffExpr = GenericAffExpr{Float64,Variable}
AffExpr() = zero(AffExpr)
AffExpr(v::Variable) = AffExpr([v], [1.], 0.)
AffExpr(v::Real) = AffExpr(Variable[], Float64[], v)
AffExpr(a::AffExpr) = AffExpr(a.vars, a.coeffs, a.constant)

Base.isempty(a::AffExpr) = (length(a.vars) == 0 && a.constant == 0.)
Base.convert(::Type{AffExpr}, v::Union{Variable,<:Real}) = AffExpr(v)

# Check all coefficients are finite, i.e. not NaN, not Inf, not -Inf
function assert_isfinite(a::AffExpr)
    coeffs = a.coeffs
    for i in 1:length(a.vars)
        isfinite(coeffs[i]) || error("Invalid coefficient $(coeffs[i]) on variable $(a.vars[i])")
    end
end

setobjective(m::Model, sense::Symbol, x::Variable) = setobjective(m, sense, convert(AffExpr,x))
function setobjective(m::Model, sense::Symbol, a::AffExpr)
    if isa(m.internalModel, MathProgBase.AbstractNonlinearModel)
        # Give the correct answer when changing objectives in an NLP.
        # A better approach would be to update and reuse the evaluator
        m.internalModelLoaded = false
    end
    if length(m.obj.qvars1) != 0
        # Go through the quadratic path so that we properly clear
        # current quadratic terms.
        setobjective(m, sense, convert(QuadExpr,a))
    else
        setobjectivesense(m, sense)
        m.obj = convert(QuadExpr,a)
    end
end

# Copy an affine expression to a new model by converting all the
# variables to the new model's variables
function Base.copy(a::AffExpr, new_model::Model)
    AffExpr(copy(a.vars, new_model), copy(a.coeffs), a.constant)
end

function getvalue(a::AffExpr)
    ret = a.constant
    for it in 1:length(a.vars)
        ret += a.coeffs[it] * getvalue(a.vars[it])
    end
    ret
end

##########################################################################
# GenericRangeConstraint
# l ≤ ∑ aᵢ xᵢ ≤ u
# The constant part of the internal expression is assumed to be zero
mutable struct GenericRangeConstraint{TermsType} <: AbstractConstraint
    terms::TermsType
    lb::Float64
    ub::Float64
end

#  b ≤ expr ≤ b   →   ==
# -∞ ≤ expr ≤ u   →   <=
#  l ≤ expr ≤ ∞   →   >=
#  l ≤ expr ≤ u   →   range
function sense(c::GenericRangeConstraint)
    if c.lb != -Inf
        if c.ub != Inf
            if c.ub == c.lb
                return :(==)
            else
                return :range
            end
        else
                return :(>=)
        end
    else #if c.lb == -Inf
        c.ub == Inf && error("'Free' constraint sense not supported")
        return :(<=)
    end
end

function rhs(c::GenericRangeConstraint)
    s = sense(c)
    s == :range && error("Range constraints do not have a well-defined RHS")
    s == :(<=) ? c.ub : c.lb
end

# Alias for AffExpr
const LinearConstraint = GenericRangeConstraint{AffExpr}

function Base.copy(c::LinearConstraint, new_model::Model)
    return LinearConstraint(copy(c.terms, new_model), c.lb, c.ub)
end

function addconstraint(m::Model, c::LinearConstraint)
    push!(m.linconstr,c)
    if m.internalModelLoaded
        if hasmethod(MathProgBase.addconstr!, (typeof(m.internalModel),Vector{Int},Vector{Float64},Float64,Float64))
            assert_isfinite(c.terms)
            indices, coeffs = merge_duplicates(Cint, c.terms, m.indexedVector, m)
            MathProgBase.addconstr!(m.internalModel,indices,coeffs,c.lb,c.ub)
        else
            warn_once("Solver does not appear to support adding constraints to an existing model. JuMP's internal model will be discarded.")
            m.internalModelLoaded = false
        end
    end
    return LinConstrRef(m,length(m.linconstr))
end
addconstraint(m::Model, c::Array{LinearConstraint}) =
    error("The operators <=, >=, and == can only be used to specify scalar constraints. If you are trying to add a vectorized constraint, use the element-wise dot comparison operators (.<=, .>=, or .==) instead")

function addVectorizedConstraint(m::Model, v::Array{LinearConstraint})
    ret = Array{LinConstrRef}(undef, size(v))
    for I in eachindex(v)
        ret[I] = addconstraint(m, v[I])
    end
    ret
end

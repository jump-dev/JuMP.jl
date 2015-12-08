#  Copyright 2015, Iain Dunning, Joey Huchette, Miles Lubin, and contributors
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
type GenericAffExpr{CoefType,VarType}
    vars::Vector{VarType}
    coeffs::Vector{CoefType}
    constant::CoefType
end
coeftype{C,V}(::GenericAffExpr{C,V}) = C

Base.zero{C,V}(::Type{GenericAffExpr{C,V}}) = GenericAffExpr{C,V}(V[],C[],zero(C))
Base.one{ C,V}(::Type{GenericAffExpr{C,V}}) = GenericAffExpr{C,V}(V[],C[], one(C))
Base.zero(a::GenericAffExpr) = zero(typeof(a))
Base.one( a::GenericAffExpr) =  one(typeof(a))
Base.copy(a::GenericAffExpr) = GenericAffExpr(copy(a.vars),copy(a.coeffs),copy(a.constant))

# Iterator protocol - iterates over tuples (aᵢ,xᵢ)
Base.start(aff::GenericAffExpr) = 1
Base.done( aff::GenericAffExpr, state::Int) = state > length(aff.vars)
Base.next( aff::GenericAffExpr, state::Int) = ((aff.coeffs[state], aff.vars[state]), state+1)
function Base.in{C,V}(x::V, aff::GenericAffExpr{C,V})
    acc = zero(C)
    for (coef,term) in aff
        if isequal(x, term)
            acc += coef
        end
    end
    return !(acc == zero(C))
end


# More efficient ways to grow an affine expression
# Add a single term to an affine expression
function Base.push!{C,V}(aff::GenericAffExpr{C,V}, new_coeff::C, new_var::V)
    push!(aff.coeffs, new_coeff)
    push!(aff.vars, new_var)
    aff
end
# Add an affine expression to an existing affine expression
function Base.append!{C,V}(aff::GenericAffExpr{C,V}, other::GenericAffExpr{C,V})
    append!(aff.vars, other.vars)
    append!(aff.coeffs, other.coeffs)
    aff.constant += other.constant
    aff
end
# For consistency, allow appending constants and individual variables
Base.append!{C}(aff::GenericAffExpr{C,C}, other::C) = error() # for ambiguity
function Base.append!{C,V}(aff::GenericAffExpr{C,V}, other::C)
    aff.constant += other
    aff
end
Base.append!{C,V}(aff::GenericAffExpr{C,V}, other::V) = push!(aff,one(C),other)

function Base.isequal{C,V}(aff::GenericAffExpr{C,V},other::GenericAffExpr{C,V})
    isequal(aff.constant, other.constant)  || return false
    length(aff.vars) == length(other.vars) || return false
    for i in 1:length(aff.vars)
        isequal(aff.vars[i],   other.vars[i])   || return false
        isequal(aff.coeffs[i], other.coeffs[i]) || return false
    end
    return true
end


# Alias for (Float64, Variable), the specific GenericAffExpr used by JuMP
typealias AffExpr GenericAffExpr{Float64,Variable}
AffExpr() = zero(AffExpr)

Base.isempty(a::AffExpr) = (length(a.vars) == 0 && a.constant == 0.)
Base.convert(::Type{AffExpr}, v::Variable) = AffExpr([v], [1.], 0.)
Base.convert(::Type{AffExpr}, v::Real) = AffExpr(Variable[], Float64[], v)

# Check all coefficients are finite, i.e. not NaN, not Inf, not -Inf
function assert_isfinite(a::AffExpr)
    coeffs = a.coeffs
    for i in 1:length(a.vars)
        isfinite(coeffs[i]) || error("Invalid coefficient $(coeffs[i]) on variable $(a.vars[i])")
    end
end

setObjective(m::Model, sense::Symbol, x::Variable) = setObjective(m, sense, convert(AffExpr,x))
function setObjective(m::Model, sense::Symbol, a::AffExpr)
    setObjectiveSense(m, sense)
    m.obj = convert(QuadExpr,a)
end

# Copy an affine expression to a new model by converting all the
# variables to the new model's variables
function Base.copy(a::AffExpr, new_model::Model)
    AffExpr(copy(a.vars, new_model), copy(a.coeffs), a.constant)
end

function getValue(a::AffExpr)
    ret = a.constant
    for it in 1:length(a.vars)
        ret += a.coeffs[it] * getValue(a.vars[it])
    end
    ret
end
function getValue(arr::Array{AffExpr})
    ret = similar(arr, Float64)
    for I in eachindex(arr)
        ret[I] = getValue(arr[I])
    end
    ret
end

##########################################################################
# GenericRangeConstraint
# l ≤ ∑ aᵢ xᵢ ≤ u
# The constant part of the internal expression is assumed to be zero
type GenericRangeConstraint{TermsType} <: JuMPConstraint
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
typealias LinearConstraint GenericRangeConstraint{AffExpr}

function Base.copy(c::LinearConstraint, new_model::Model)
    return LinearConstraint(copy(c.terms, new_model), c.lb, c.ub)
end

function addConstraint(m::Model, c::LinearConstraint)
    push!(m.linconstr,c)
    if m.internalModelLoaded
        if method_exists(MathProgBase.addconstr!, (typeof(m.internalModel),Vector{Int},Vector{Float64},Float64,Float64))
            assert_isfinite(c.terms)
            indices, coeffs = merge_duplicates(Cint, c.terms, m.indexedVector, m)
            MathProgBase.addconstr!(m.internalModel,indices,coeffs,c.lb,c.ub)
        else
            Base.warn_once("Solver does not appear to support adding constraints to an existing model. Hot-start is disabled.")
            m.internalModelLoaded = false
        end
    end
    return ConstraintRef{LinearConstraint}(m,length(m.linconstr))
end
addConstraint(m::Model, c::Array{LinearConstraint}) =
    error("The operators <=, >=, and == can only be used to specify scalar constraints. If you are trying to add a vectorized constraint, use the element-wise dot comparison operators (.<=, .>=, or .==) instead")

function addVectorizedConstraint(m::Model, v::Array{LinearConstraint})
    ret = Array(ConstraintRef{LinearConstraint}, size(v))
    for I in eachindex(v)
        ret[I] = addConstraint(m, v[I])
    end
    ret
end

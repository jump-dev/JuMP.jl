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
type GenericAffExpr{CoefType,VarType} <: AbstractJuMPScalar
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

# Old iterator protocol - iterates over tuples (aᵢ,xᵢ)
struct LinearTermIterator{GAE<:GenericAffExpr}
    aff::GAE
end

"""
    linearterms(aff::GenericAffExpr)

Provides an iterator over the `(a_i::C,x_i::V)` terms in affine expression ``\\sum_i a_i x_i + b``.
"""
linearterms(aff::GenericAffExpr) = LinearTermIterator(aff)

Base.start(lti::LinearTermIterator) = 1
Base.done( lti::LinearTermIterator, state::Int) = state > length(lti.aff.vars)
Base.next( lti::LinearTermIterator, state::Int) = ((lti.aff.coeffs[state], lti.aff.vars[state]), state+1)

"""
    Base.push!{C,V}(aff::GenericAffExpr{C,V}, new_coeff::C, new_var::V)

An efficient way to grow an affine expression by one term. For example, to add `5x` to an existing expression `aff`, use `push!(aff, 5.0, x)`. This is significantly more efficient than `aff += 5.0*x`.
"""
function Base.push!{C,V}(aff::GenericAffExpr{C,V}, new_coeff::C, new_var::V)
    push!(aff.coeffs, new_coeff)
    push!(aff.vars, new_var)
    aff
end

# Add an affine expression to an existing affine expression
"""
    Base.append!{C,V}(aff::GenericAffExpr{C,V}, other::GenericAffExpr{C,V})

Efficiently append the terms of an affine expression to an existing affine expression. For example, given `aff = 5.0*x` and `other = 7.0*y + 3.0*z`, we can grow `aff` using `append!(aff, other)` which results in `aff` equaling `5x + 7y + 3z`. This is significantly more efficient than using `aff += other`.
"""
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
function Base.append!{C,V}(aff::GenericAffExpr{C,V}, other::Real)
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

# Check if two AffExprs are equal regardless of the order, and after merging duplicates
# Mostly useful for testing.
function isequal_canonical(aff::GenericAffExpr{C,V}, other::GenericAffExpr{C,V}) where {C,V}
    function canonicalize(a)
        d = Dict{V,C}()
        @assert length(a.vars) == length(a.coeffs)
        for k in 1:length(a.vars)
            d[a.vars[k]] = a.coeffs[k] + get(d, a.vars[k], zero(C))
        end
        return d
    end
    d1 = canonicalize(aff)
    d2 = canonicalize(other)
    return isequal(d1,d1) && aff.constant == other.constant
end


# Alias for (Float64, Variable), the specific GenericAffExpr used by JuMP
const AffExpr = GenericAffExpr{Float64,Variable}
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

function MOI.ScalarAffineFunction(a::AffExpr)
    return MOI.ScalarAffineFunction(instancereference.(a.vars), a.coeffs, a.constant)
end

function AffExpr(m::Model, f::MOI.ScalarAffineFunction)
    return AffExpr(Variable.(m,f.variables), f.coefficients, f.constant)
end

"""
    _fillvaf!(outputindex, variables, coefficients, offset::Int, oi::Int, coeff, aff::AffExpr)

Fills the vectors outputindex, variables, coefficients at indices starting at `offset+1` with the terms of `aff`.
The output index for all terms is `oi`.
"""
function _fillvaf!(outputindex, variables, coefficients, offset::Int, oi::Int, aff::AffExpr)
    for i in 1:length(aff.vars)
        outputindex[offset+i] = oi
        variables[offset+i] = instancereference(aff.vars[i])
        coefficients[offset+i] = aff.coeffs[i]
    end
    offset + length(aff.vars)
end

function MOI.VectorAffineFunction(affs::Vector{AffExpr})
    len = sum(aff -> length(aff.vars), affs)
    outputindex = Vector{Int}(len)
    variables = Vector{MOI.VariableReference}(len)
    coefficients = Vector{Float64}(len)
    constant = Vector{Float64}(length(affs))
    offset = 0
    for (i, aff) in enumerate(affs)
        constant[i] = aff.constant
        offset = _fillvaf!(outputindex, variables, coefficients, offset, i, aff)
    end
    MOI.VectorAffineFunction(outputindex, variables, coefficients, constant)
end

# TODO this could be interpreted as a SingleVariable objective, but that should require explict syntax
setobjective(m::Model, sense::Symbol, x::Variable) = setobjective(m, sense, convert(AffExpr,x))


function setobjective(m::Model, sense::Symbol, a::AffExpr)
    @assert !m.solverinstanceattached # TODO
    if sense == :Min
        moisense = MOI.MinSense
    else
        @assert sense == :Max
        moisense = MOI.MaxSense
    end
    MOI.setobjective!(m.instance, moisense, MOI.ScalarAffineFunction(a))
    nothing
end

"""
    objectivefunction(m::Model, ::Type{AffExpr})

Return an `AffExpr` object representing the objective function.
Error if the objective is not linear.
"""
function objectivefunction(m::Model, ::Type{AffExpr})
    f = MOI.getattribute(m.instance, MOI.ObjectiveFunction())::MOI.ScalarAffineFunction
    return AffExpr(m, f)
end


# Copy an affine expression to a new model by converting all the
# variables to the new model's variables
function Base.copy(a::AffExpr, new_model::Model)
    AffExpr(copy(a.vars, new_model), copy(a.coeffs), a.constant)
end

# """
#     getvalue(a::AffExpr)
#
# Evaluate an `AffExpr` given the current solution values.
# """
# function getvalue(a::AffExpr)
#     ret = a.constant
#     for it in 1:length(a.vars)
#         ret += a.coeffs[it] * getvalue(a.vars[it])
#     end
#     ret
# end

# TODO GenericLinearConstraint

struct LinearConstraint{S <: MOI.AbstractScalarSet} <: AbstractConstraint
    func::AffExpr
    set::S
end

# function Base.copy(c::LinearConstraint, new_model::Model)
#     return LinearConstraint(copy(c.terms, new_model), c.lb, c.ub)
# end

"""
    addconstraint(m::Model, c::LinearConstraint)

Add a linear constraint to `Model m`.
"""
function addconstraint(m::Model, c::LinearConstraint)
    @assert !m.solverinstanceattached # TODO
    cref = MOI.addconstraint!(m.instance, MOI.ScalarAffineFunction(c.func), c.set)
    return ConstraintRef(m, cref)
end
addconstraint(m::Model, c::Array{LinearConstraint}) =
    error("The operators <=, >=, and == can only be used to specify scalar constraints. If you are trying to add a vectorized constraint, use the element-wise dot comparison operators (.<=, .>=, or .==) instead")

#  Copyright 2017, Iain Dunning, Joey Huchette, Miles Lubin, and contributors
#  This Source Code Form is subject to the terms of the Mozilla Public
#  License, v. 2.0. If a copy of the MPL was not distributed with this
#  file, You can obtain one at http://mozilla.org/MPL/2.0/.
#############################################################################
# JuMP
# An algebraic modeling language for Julia
# See http://github.com/JuliaOpt/JuMP.jl
#############################################################################
# src/quadexpr.jl
# Defines all types relating to expressions with a quadratic and affine part
# - GenericQuadExpr             ∑qᵢⱼ xᵢⱼ  +  ∑ aᵢ xᵢ  +  c
#   - QuadExpr                  Alias for (Float64, Variable)
# - QuadExprConstraint       ∑qᵢⱼ xᵢⱼ  +  ∑ aᵢ xᵢ  +  c  in set
# Operator overloads in src/operators.jl
#############################################################################


#############################################################################
# GenericQuadExpr
# ∑qᵢⱼ xᵢⱼ  +  ∑ aᵢ xᵢ  +  c
mutable struct GenericQuadExpr{CoefType,VarType} <: AbstractJuMPScalar
    qvars1::Vector{VarType}
    qvars2::Vector{VarType}
    qcoeffs::Vector{CoefType}
    aff::GenericAffExpr{CoefType,VarType}
end

Base.isempty(q::GenericQuadExpr) = (length(q.qvars1) == 0 && isempty(q.aff))
Base.zero(::Type{GenericQuadExpr{C,V}}) where {C,V} = GenericQuadExpr(V[], V[], C[], zero(GenericAffExpr{C,V}))
Base.one(::Type{GenericQuadExpr{C,V}}) where {C,V}  = GenericQuadExpr(V[], V[], C[],  one(GenericAffExpr{C,V}))
Base.zero(q::GenericQuadExpr) = zero(typeof(q))
Base.one(q::GenericQuadExpr)  =  one(typeof(q))
Base.copy(q::GenericQuadExpr) = GenericQuadExpr(copy(q.qvars1),copy(q.qvars2),copy(q.qcoeffs),copy(q.aff))

function Base.append!(q::GenericQuadExpr{T,S}, other::GenericQuadExpr{T,S}) where {T,S}
    append!(q.qvars1, other.qvars1)
    append!(q.qvars2, other.qvars2)
    append!(q.qcoeffs, other.qcoeffs)
    append!(q.aff, other.aff)
    q
end

function assert_isfinite(q::GenericQuadExpr)
    assert_isfinite(q.aff)
    for i in 1:length(q.qcoeffs)
        isfinite(q.qcoeffs[i]) || error("Invalid coefficient $(q.qcoeffs[i]) on quadratic term $(q.qvars1[i])*$(q.qvars2[i])")
    end
end

function Base.isequal(q::GenericQuadExpr{T,S},other::GenericQuadExpr{T,S}) where {T,S}
    isequal(q.aff,other.aff)   || return false
    length(q.qvars1) == length(other.qvars1) || return false
    for i in 1:length(q.qvars1)
        isequal(q.qvars1[i],  other.qvars1[i]) || return false
        isequal(q.qvars2[i],  other.qvars2[i]) || return false
        isequal(q.qcoeffs[i],other.qcoeffs[i]) || return false
    end
    return true
end

# Check if two QuadExprs are equal regardless of the order, and after merging duplicates
# Mostly useful for testing.
function isequal_canonical(quad::GenericQuadExpr{CoefType,VarType}, other::GenericQuadExpr{CoefType,VarType}) where {CoefType,VarType}
    function canonicalize(q)
        d = Dict{Set{VarType},CoefType}()
        @assert length(q.qvars1) == length(q.qvars2) == length(q.qcoeffs)
        for (v1,v2,c) in zip(q.qvars1, q.qvars2, q.qcoeffs)
            vset = Set((v1,v2))
            d[vset] = c + get(d, vset, zero(CoefType))
        end
        return d
    end
    d1 = canonicalize(quad)
    d2 = canonicalize(other)
    return isequal(d1,d2) && isequal_canonical(quad.aff, other.aff)
end

# Alias for (Float64, Variable)
const QuadExpr = GenericQuadExpr{Float64,Variable}
Base.convert(::Type{QuadExpr}, v::Union{Real,Variable,AffExpr}) = QuadExpr(Variable[], Variable[], Float64[], AffExpr(v))
QuadExpr() = zero(QuadExpr)

function MOI.ScalarQuadraticFunction(q::QuadExpr)
    scaledcoef = [(v1 == v2) ? 2coef : coef for (v1,v2,coef) in zip(q.qvars1, q.qvars2, q.qcoeffs)]
    return MOI.ScalarQuadraticFunction(index.(q.aff.vars), q.aff.coeffs, index.(q.qvars1), index.(q.qvars2), scaledcoef, q.aff.constant)
end

function QuadExpr(m::Model, f::MOI.ScalarQuadraticFunction)
    scaledcoef = [(v1 == v2) ? coef/2 : coef for (v1, v2, coef) in zip(f.quadratic_rowvariables, f.quadratic_colvariables, f.quadratic_coefficients)]
    return QuadExpr(Variable.(m,f.quadratic_rowvariables), Variable.(m, f.quadratic_colvariables), scaledcoef, AffExpr(Variable.(m,f.affine_variables), f.affine_coefficients, f.constant))
end

function setobjective(m::Model, sense::Symbol, a::QuadExpr)
    if sense == :Min
        moisense = MOI.MinSense
    else
        @assert sense == :Max
        moisense = MOI.MaxSense
    end
    MOI.set!(m.moibackend, MOI.ObjectiveSense(), moisense)
    MOI.set!(m.moibackend, MOI.ObjectiveFunction{MOI.ScalarQuadraticFunction{Float64}}(), MOI.ScalarQuadraticFunction(a))
    nothing
end

"""
    objectivefunction(m::Model, ::Type{QuadExpr})

Return a `QuadExpr` object representing the objective function.
Error if the objective is not quadratic.
"""
function objectivefunction(m::Model, ::Type{QuadExpr})
    f = MOI.get(m.moibackend, MOI.ObjectiveFunction{MOI.ScalarQuadraticFunction{Float64}}())::MOI.ScalarQuadraticFunction
    return QuadExpr(m, f)
end

# Copy a quadratic expression to a new model by converting all the
# variables to the new model's variables
function Base.copy(q::QuadExpr, new_model::Model)
    QuadExpr(copy(q.qvars1, new_model), copy(q.qvars2, new_model),
                copy(q.qcoeffs), copy(q.aff, new_model))
end

# TODO: resultvalue for QuadExpr

##########################################################################
# TODO: GenericQuadExprConstraint

struct QuadExprConstraint{S <: MOI.AbstractScalarSet} <: AbstractConstraint
    func::QuadExpr
    set::S
end

moi_function_and_set(c::QuadExprConstraint) = (MOI.ScalarQuadraticFunction(c.func), c.set)

function constraintobject(cr::ConstraintRef{Model}, ::Type{QuadExpr}, ::Type{SetType}) where {SetType <: MOI.AbstractScalarSet}
    m = cr.m
    f = MOI.get(m.moibackend, MOI.ConstraintFunction(), index(cr))::MOI.ScalarQuadraticFunction
    s = MOI.get(m.moibackend, MOI.ConstraintSet(), index(cr))::SetType
    return QuadExprConstraint(QuadExpr(m, f), s)
end

# TODO: VectorQuadExprConstraint

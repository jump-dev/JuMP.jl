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
type GenericQuadExpr{CoefType,VarType} <: AbstractJuMPScalar
    qvars1::Vector{VarType}
    qvars2::Vector{VarType}
    qcoeffs::Vector{CoefType}
    aff::GenericAffExpr{CoefType,VarType}
end
coeftype{C,V}(::GenericQuadExpr{C,V}) = C

Base.isempty(q::GenericQuadExpr) = (length(q.qvars1) == 0 && isempty(q.aff))
Base.zero{C,V}(::Type{GenericQuadExpr{C,V}}) = GenericQuadExpr(V[], V[], C[], zero(GenericAffExpr{C,V}))
Base.one{C,V}(::Type{GenericQuadExpr{C,V}})  = GenericQuadExpr(V[], V[], C[],  one(GenericAffExpr{C,V}))
Base.zero(q::GenericQuadExpr) = zero(typeof(q))
Base.one(q::GenericQuadExpr)  =  one(typeof(q))
Base.copy(q::GenericQuadExpr) = GenericQuadExpr(copy(q.qvars1),copy(q.qvars2),copy(q.qcoeffs),copy(q.aff))

function Base.append!{T,S}(q::GenericQuadExpr{T,S}, other::GenericQuadExpr{T,S})
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

function Base.isequal{T,S}(q::GenericQuadExpr{T,S},other::GenericQuadExpr{T,S})
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
    return isequal(d1,d1) && isequal_canonical(quad.aff, other.aff)
end

# Alias for (Float64, Variable)
const QuadExpr = GenericQuadExpr{Float64,Variable}
Base.convert(::Type{QuadExpr}, v::Union{Real,Variable,AffExpr}) = QuadExpr(Variable[], Variable[], Float64[], AffExpr(v))
QuadExpr() = zero(QuadExpr)

function MOI.ScalarQuadraticFunction(q::QuadExpr)
    return MOI.ScalarQuadraticFunction(instancereference.(q.aff.vars), q.aff.coeffs, instancereference.(q.qvars1), instancereference.(q.qvars2), q.qcoeffs, q.aff.constant)
end

function QuadExpr(m::Model, f::MOI.ScalarQuadraticFunction)
    return QuadExpr(Variable.(m,f.quadratic_rowvariables), Variable.(m, f.quadratic_colvariables), f.quadratic_coefficients, AffExpr(Variable.(m,f.affine_variables), f.affine_coefficients, f.constant))
end

function setobjective(m::Model, sense::Symbol, a::QuadExpr)
    @assert !m.solverinstanceattached # TODO
    if sense == :Min
        moisense = MOI.MinSense
    else
        @assert sense == :Max
        moisense = MOI.MaxSense
    end
    MOI.set!(m.instance, MOI.ObjectiveSense(), moisense)
    MOI.set!(m.instance, MOI.ObjectiveFunction(), MOI.ScalarQuadraticFunction(a))
    nothing
end

"""
    objectivefunction(m::Model, ::Type{QuadExpr})

Return a `QuadExpr` object representing the objective function.
Error if the objective is not quadratic.
"""
function objectivefunction(m::Model, ::Type{QuadExpr})
    f = MOI.get(m.instance, MOI.ObjectiveFunction())::MOI.ScalarQuadraticFunction
    return QuadExpr(m, f)
end

# Copy a quadratic expression to a new model by converting all the
# variables to the new model's variables
function Base.copy(q::QuadExpr, new_model::Model)
    QuadExpr(copy(q.qvars1, new_model), copy(q.qvars2, new_model),
                copy(q.qcoeffs), copy(q.aff, new_model))
end

# """
#     getvalue(a::QuadExpr)
#
# Evaluate a `QuadExpr` given the current solution values.
# """
# function getvalue(a::QuadExpr)
#     ret = getvalue(a.aff)
#     for it in 1:length(a.qvars1)
#         ret += a.qcoeffs[it] * getvalue(a.qvars1[it]) * getvalue(a.qvars2[it])
#     end
#     return ret
# end
# getvalue(arr::Array{QuadExpr}) = map(getvalue, arr)



##########################################################################
# TODO: GenericQuadExprConstraint

struct QuadExprConstraint{S <: MOI.AbstractScalarSet} <: AbstractConstraint
    func::QuadExpr
    set::S
end

"""
    addconstraint(m::Model, c::QuadExprConstraint)

Add a `QuadExpr` constraint to `Model m`.
"""
function addconstraint(m::Model, c::QuadExprConstraint)
    @assert !m.solverinstanceattached # TODO
    cref = MOI.addconstraint!(m.instance, MOI.ScalarQuadraticFunction(c.func), c.set)
    return ConstraintRef(m, cref)
end

function constraintobject(cref::ConstraintRef{Model}, ::Type{QuadExpr}, ::Type{SetType}) where {SetType <: MOI.AbstractScalarSet}
    m = cref.m
    f = MOI.get(m.instance, MOI.ConstraintFunction(), cref.instanceref)::MOI.ScalarQuadraticFunction
    s = MOI.get(m.instance, MOI.ConstraintSet(), cref.instanceref)::SetType
    return QuadExprConstraint(QuadExpr(m, f), s)
end

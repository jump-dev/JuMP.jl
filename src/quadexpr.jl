#  Copyright 2015, Iain Dunning, Joey Huchette, Miles Lubin, and contributors
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
# - GenericQuadConstraint       ∑qᵢⱼ xᵢⱼ  +  ∑ aᵢ xᵢ  +  c  [≤,≥]  0
#   - QuadConstraint            Alias for (Float64, Variable)
# Operator overloads in src/operators.jl
#############################################################################


#############################################################################
# GenericQuadExpr
# ∑qᵢⱼ xᵢⱼ  +  ∑ aᵢ xᵢ  +  c
type GenericQuadExpr{CoefType,VarType}
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

# Alias for (Float64, Variable)
typealias QuadExpr GenericQuadExpr{Float64,Variable}
Base.convert(::Type{QuadExpr}, v::Union{Real,Variable,AffExpr}) = QuadExpr(Variable[], Variable[], Float64[], AffExpr(v))
QuadExpr() = zero(QuadExpr)

function setObjective(m::Model, sense::Symbol, q::QuadExpr)
    m.obj = q
    if m.internalModelLoaded
        if method_exists(MathProgBase.setquadobjterms!, (typeof(m.internalModel), Vector{Cint}, Vector{Cint}, Vector{Float64}))
            verify_ownership(m, m.obj.qvars1)
            verify_ownership(m, m.obj.qvars2)
            MathProgBase.setquadobjterms!(m.internalModel, Cint[v.col for v in m.obj.qvars1], Cint[v.col for v in m.obj.qvars2], m.obj.qcoeffs)
        else
            # we don't (yet) support hot-starting QCQP solutions
            Base.warn_once("JuMP does not yet support adding a quadratic objective to an existing model. Hot-start is disabled.")
            m.internalModelLoaded = false
        end
    end
    setObjectiveSense(m, sense)
end

# Copy a quadratic expression to a new model by converting all the
# variables to the new model's variables
function Base.copy(q::QuadExpr, new_model::Model)
    QuadExpr(copy(q.qvars1, new_model), copy(q.qvars2, new_model),
                copy(q.qcoeffs), copy(q.aff, new_model))
end

function getValue(a::QuadExpr)
    ret = getValue(a.aff)
    for it in 1:length(a.qvars1)
        ret += a.qcoeffs[it] * getValue(a.qvars1[it]) * getValue(a.qvars2[it])
    end
    return ret
end
getValue(arr::Array{QuadExpr}) = map(getValue, arr)



##########################################################################
# GenericQuadConstraint
# ∑qᵢⱼ xᵢⱼ  +  ∑ aᵢ xᵢ  +  c  [≤,≥]  0
# As RHS is implicitly taken to be zero, we store only LHS and sense
type GenericQuadConstraint{QuadType} <: JuMPConstraint
    terms::QuadType
    sense::Symbol
end
Base.copy{CON<:GenericQuadConstraint}(c::CON, new_model::Model) = CON(copy(c.terms, new_model), c.sense)


# Alias for (Float64, Variable)
typealias QuadConstraint GenericQuadConstraint{QuadExpr}

function Base.copy(c::QuadConstraint, new_model::Model)
    return QuadConstraint(copy(c.terms, new_model), c.sense)
end

function addConstraint(m::Model, c::QuadConstraint)
    push!(m.quadconstr,c)
    if m.internalModelLoaded
        if method_exists(MathProgBase.addquadconstr!, (typeof(m.internalModel),
                                                       Vector{Cint},
                                                       Vector{Float64},
                                                       Vector{Cint},
                                                       Vector{Cint},
                                                       Vector{Float64},
                                                       Char,
                                                       Float64))
            if !((s = string(c.sense)[1]) in ['<', '>', '='])
                error("Invalid sense for quadratic constraint")
            end
            terms = c.terms
            verify_ownership(m, terms.qvars1)
            verify_ownership(m, terms.qvars2)
            MathProgBase.addquadconstr!(m.internalModel,
                                        Cint[v.col for v in c.terms.aff.vars],
                                        c.terms.aff.coeffs,
                                        Cint[v.col for v in c.terms.qvars1],
                                        Cint[v.col for v in c.terms.qvars2],
                                        c.terms.qcoeffs,
                                        s,
                                        -c.terms.aff.constant)
        else
            Base.warn_once("Solver does not appear to support adding quadratic constraints to an existing model. Hot-start is disabled.")
            m.internalModelLoaded = false
        end
    end
    return ConstraintRef{QuadConstraint}(m,length(m.quadconstr))
end
addConstraint(m::Model, c::Array{QuadConstraint}) =
    error("Vectorized constraint added without elementwise comparisons. Try using one of (.<=,.>=,.==).")

function addVectorizedConstraint(m::Model, v::Array{QuadConstraint})
    ret = Array(ConstraintRef{QuadConstraint}, size(v))
    for I in eachindex(v)
        ret[I] = addConstraint(m, v[I])
    end
    ret
end

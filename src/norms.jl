#  Copyright 2015, Iain Dunning, Joey Huchette, Miles Lubin, and contributors
#  This Source Code Form is subject to the terms of the Mozilla Public
#  License, v. 2.0. If a copy of the MPL was not distributed with this
#  file, You can obtain one at http://mozilla.org/MPL/2.0/.
#############################################################################
# JuMP
# An algebraic modeling language for Julia
# See http://github.com/JuliaOpt/JuMP.jl
#############################################################################
# src/norms.jl
# Defines all types relating to norms
# - GenericNorm                 Lp norm of vector of GenericAffExprs
#   - Norm                      Lp norm of vector of AffExprs
# - GenericNormExpr             Expression of norms and affine expressions
#   - GenericSOCExpr            Alias for P=2
#       - SOCExpr               Alias for P=2 and AffExprs
# - GenericSOCConstraint        Constraint of form α||Ax-b||₂ + cᵀx + d ≤ 0
#   - SOCConstraint             Alias for AffExprs
# Operator overloads in src/operators.jl
#############################################################################


#############################################################################
# GenericNorm
# Container for Lp norms, ‖Ax‖p
type GenericNorm{P,C,V}
    terms::Vector{GenericAffExpr{C,V}}
end
# Preferred constructor, validates the norm type based on the expression
function GenericNorm{C,V}(P, terms::Vector{GenericAffExpr{C,V}})
    if C == Float64 && V == Variable
        # JuMP doesn't support anything else than the L2 Norm
        # Throw an error now before going any further
        P ==   1 && error("JuMP doesn't support L₁ norms.")
        P == Inf && error("JuMP doesn't support L∞ norms.")
        P !=   2 && error("JuMP only supports L₂ norms.")
    end
    GenericNorm{P,C,V}(terms)
end
Base.copy{P,C,V}(x::GenericNorm{P,C,V}) = GenericNorm{P,C,V}(copy(x.terms))

# Handle the norm() function by flattening arguments into a vector
Base.norm{V<:AbstractJuMPScalar}(x::V,           p=2) = vecnorm(x,p)
Base.norm{V<:AbstractJuMPScalar}(x::Array{V},    p=2) = vecnorm(x,p)
Base.norm{V<:AbstractJuMPScalar}(x::JuMPArray{V},p=2) = vecnorm(x,p)
Base.norm{V<:AbstractJuMPScalar}(x::JuMPDict{V}, p=2) = vecnorm(x,p)
Base.norm{C,V}(x::GenericAffExpr{C,V},           p=2) = vecnorm(x,p)
Base.norm{C,V}(x::Array{GenericAffExpr{C,V}},    p=2) = vecnorm(x,p)
Base.norm{C,V}(x::JuMPArray{GenericAffExpr{C,V}},p=2) = vecnorm(x,p)
Base.norm{C,V}(x::JuMPDict{GenericAffExpr{C,V}}, p=2) = vecnorm(x,p)

_vecaff(C,V,x) = map(GenericAffExpr{C,V},vec(x))
Base.vecnorm{V<:AbstractJuMPScalar}(x::V,           p=2) = GenericNorm(p, [GenericAffExpr{Float64,V}(x)] )
Base.vecnorm{V<:AbstractJuMPScalar}(x::Array{V},    p=2) = GenericNorm(p, _vecaff(Float64,V,x) )
Base.vecnorm{V<:AbstractJuMPScalar}(x::JuMPArray{V},p=2) = GenericNorm(p, _vecaff(Float64,V,x.innerArray) )
Base.vecnorm{V<:AbstractJuMPScalar}(x::JuMPDict{V}, p=2) = GenericNorm(p, _vecaff(Float64,V,collect(values(x))) )
Base.vecnorm{C,V}(x::GenericAffExpr{C,V},           p=2) = GenericNorm(p, [x])
Base.vecnorm{C,V}(x::Array{GenericAffExpr{C,V}},    p=2) = GenericNorm(p, vec(x))
Base.vecnorm{C,V}(x::JuMPArray{GenericAffExpr{C,V}},p=2) = GenericNorm(p, vec(x.innerArray))
Base.vecnorm{C,V}(x::JuMPDict{GenericAffExpr{C,V}}, p=2) = GenericNorm(p, collect(values(x)))

# Called by the parseNorm macro for e.g. norm2{...}
# If the arguments are tightly typed, just pass to the constructor
_build_norm{C,V}(P,terms::Vector{GenericAffExpr{C,V}}) = GenericNorm(P,terms)
# The terms vector produced by the macro may not be tightly typed,
# so we need to repackage in a tighter typed vector ourselves.
# This function is needed for performance reasons on only 0.3, as
# the following works just as well on 0.4:
# _build_norm(Lp, terms::Vector{GenericAffExpr}) = _build_norm(Lp, [terms...])
function _build_norm(Lp, terms::Vector{GenericAffExpr})
    if length(terms) == 0
        _build_norm(Lp,AffExpr[])  # Punt
    else
        new_terms = Array(typeof(terms[1]), length(terms))
        for i in 1:length(terms)
            new_terms[i] = terms[i]
        end
        _build_norm(Lp,new_terms)
    end
end

# Alias for AffExprs. Short-hand used in operator overloads, etc.
typealias Norm{P} GenericNorm{P,Float64,Variable}


##########################################################################
# GenericNormExpr
# Container for expressions containing GenericNorm and GenericAffExpr
type GenericNormExpr{P,C,V}
    norm::GenericNorm{P,C,V}
    coeff::C
    aff::GenericAffExpr{C,V}
end

GenericNormExpr{P,C,V}(norm::GenericNorm{P,C,V}) =
    GenericNormExpr{C,V}(norm, one(C), zero(GenericAffExpr{C,V}))
Base.copy{P,C,V}(x::GenericNormExpr{P,C,V}) =
    GenericNormExpr{P,C,V}(copy(x.norm), copy(x.coeff), copy(x.aff))
Base.convert{P,C,V}(::Type{GenericNormExpr{P,C,V}}, x::GenericNorm{P,C,V}) =
    GenericNormExpr{P,C,V}(x, one(C), zero(GenericAffExpr{C,V}))

# Alias for ‖Ax‖₂ case
typealias GenericSOCExpr{C,V} GenericNormExpr{2,C,V}

# Alias for ‖Ax‖₂ and AffExpr case
typealias SOCExpr GenericSOCExpr{Float64,Variable}


##########################################################################
# GenericSOCConstraint
# Second-order cone constraint of form 
# α||Ax-b||₂ + cᵀx + d ≤ 0
type GenericSOCConstraint{T<:GenericSOCExpr} <: JuMPConstraint
    normexpr::T
    function GenericSOCConstraint{T}(normexpr::T)
        if normexpr.coeff < 0
            # The coefficient in front of the norm is negative, which
            # means we have `norm >= c`, which is not convex.
            error("Invalid second-order cone constraint $(normexpr) ≤ 0")
        end
        new(normexpr)
    end
end

# Alias for the AffExpr case
typealias SOCConstraint GenericSOCConstraint{SOCExpr}

function addConstraint(m::Model, c::SOCConstraint)
    push!(m.socconstr,c)
    m.internalModelLoaded = false
    ConstraintRef{SOCConstraint}(m,length(m.socconstr))
end
#  Copyright 2017, Iain Dunning, Joey Huchette, Miles Lubin, and contributors
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
mutable struct GenericNorm{P,C,V}
    terms::Vector{GenericAffExpr{C,V}}
end
# Preferred constructor, validates the norm type based on the expression
function GenericNorm(P, terms::Vector{GenericAffExpr{C,V}}) where {C,V}
    if C == Float64 && V == Variable
        # JuMP doesn't support anything else than the L2 Norm
        # Throw an error now before going any further
        P ==   1 && error("JuMP doesn't support L₁ norms.")
        P == Inf && error("JuMP doesn't support L∞ norms.")
        P !=   2 && error("JuMP only supports L₂ norms.")
    end
    GenericNorm{P,C,V}(terms)
end
function GenericNorm(P, terms::AbstractVector{GenericAffExpr{C,V}}) where {C, V}
    GenericNorm(P, [terms[i] for i in eachindex(terms)])
end
Base.copy(x::GenericNorm{P,C,V}) where {P,C,V} = GenericNorm{P,C,V}(copy(x.terms))

# Handle the norm() function by flattening arguments into a vector
if VERSION < v"0.7-"
    LinearAlgebra.norm(x::V,                  p::Real=2) where {V<:AbstractJuMPScalar} = vecnorm(x,p)
    LinearAlgebra.norm(x::AbstractVector{V},  p::Real=2) where {V<:AbstractJuMPScalar} = vecnorm(x,p)
    LinearAlgebra.norm(x::AbstractMatrix{V},  p::Real=2) where {V<:AbstractJuMPScalar} = vecnorm(x,p)
    LinearAlgebra.norm(x::JuMPArray{V},       p::Real=2) where {V<:AbstractJuMPScalar} = vecnorm(x,p)
    LinearAlgebra.norm(x::JuMPDict{V},        p::Real=2) where {V<:AbstractJuMPScalar} = vecnorm(x,p)
    LinearAlgebra.norm(x::GenericAffExpr{C,V},                  p::Real=2) where {C,V} = vecnorm(x,p)
    LinearAlgebra.norm(x::AbstractVector{GenericAffExpr{C,V}},  p::Real=2) where {C,V} = vecnorm(x,p)
    LinearAlgebra.norm(x::AbstractMatrix{GenericAffExpr{C,V}},  p::Real=2) where {C,V} = vecnorm(x,p)
    LinearAlgebra.norm(x::JuMPArray{GenericAffExpr{C,V}},       p::Real=2) where {C,V} = vecnorm(x,p)
    LinearAlgebra.norm(x::JuMPDict{GenericAffExpr{C,V}},        p::Real=2) where {C,V} = vecnorm(x,p)
    
    LinearAlgebra.vecnorm(x::V,                   p::Real=2) where {V<:AbstractJuMPScalar} = GenericNorm(p, [GenericAffExpr{Float64,V}(x)] )
    LinearAlgebra.vecnorm(x::AbstractArray{V},    p::Real=2) where {V<:AbstractJuMPScalar} = GenericNorm(p, _vecaff(Float64,V,x) )
    LinearAlgebra.vecnorm(x::JuMPArray{V},        p::Real=2) where {V<:AbstractJuMPScalar} = GenericNorm(p, _vecaff(Float64,V,x.innerArray) )
    LinearAlgebra.vecnorm(x::JuMPDict{V},         p::Real=2) where {V<:AbstractJuMPScalar} = GenericNorm(p, _vecaff(Float64,V,collect(values(x))) )
    LinearAlgebra.vecnorm(x::GenericAffExpr{C,V},                   p::Real=2) where {C,V} = GenericNorm(p, [x])
    LinearAlgebra.vecnorm(x::AbstractArray{GenericAffExpr{C,V}},    p::Real=2) where {C,V} = GenericNorm(p, vec(x))
    LinearAlgebra.vecnorm(x::JuMPArray{GenericAffExpr{C,V}},        p::Real=2) where {C,V} = GenericNorm(p, vec(x.innerArray))
    LinearAlgebra.vecnorm(x::JuMPDict{GenericAffExpr{C,V}},         p::Real=2) where {C,V} = GenericNorm(p, collect(values(x)))
else
    LinearAlgebra.norm(x::V,                  p::Real=2) where {V<:AbstractJuMPScalar} = GenericNorm(p, [GenericAffExpr{Float64,V}(x)] )
    LinearAlgebra.norm(x::AbstractArray{V},   p::Real=2) where {V<:AbstractJuMPScalar} = GenericNorm(p, _vecaff(Float64,V,x) )
    LinearAlgebra.norm(x::JuMPArray{V},       p::Real=2) where {V<:AbstractJuMPScalar} = GenericNorm(p, _vecaff(Float64,V,x.innerArray) )
    LinearAlgebra.norm(x::JuMPDict{V},        p::Real=2) where {V<:AbstractJuMPScalar} = GenericNorm(p, _vecaff(Float64,V,collect(values(x))) )
    LinearAlgebra.norm(x::GenericAffExpr{C,V},                  p::Real=2) where {C,V} = GenericNorm(p, [x])
    LinearAlgebra.norm(x::AbstractArray{GenericAffExpr{C,V}},   p::Real=2) where {C,V} = GenericNorm(p, vec(x))
    LinearAlgebra.norm(x::JuMPArray{GenericAffExpr{C,V}},       p::Real=2) where {C,V} = GenericNorm(p, vec(x.innerArray))
    LinearAlgebra.norm(x::JuMPDict{GenericAffExpr{C,V}},        p::Real=2) where {C,V} = GenericNorm(p, collect(values(x)))
end
_vecaff(C,V,x) = map(GenericAffExpr{C,V},vec(x))

# Called by the parseNorm macro for e.g. norm2{...}
# If the arguments are tightly typed, just pass to the constructor
_build_norm(P,terms::Vector{GenericAffExpr{C,V}}) where {C,V} = GenericNorm(P,terms)
# The terms vector produced by the macro may not be tightly typed,
# so we need to repackage in a tighter typed vector ourselves.
# This function is needed for performance reasons on only 0.3, as
# the following works just as well on 0.4:
_build_norm(Lp, terms::Vector{GenericAffExpr}) = _build_norm(Lp, [terms...])

# Alias for AffExprs. Short-hand used in operator overloads, etc.
Norm{P} = GenericNorm{P,Float64,Variable}

getvalue(n::GenericNorm{P,C,V}) where {P,C,V} = norm(getvalue(n.terms),P)

##########################################################################
# GenericNormExpr
# Container for expressions containing GenericNorm and GenericAffExpr
mutable struct GenericNormExpr{P,C,V}
    norm::GenericNorm{P,C,V}
    coeff::C
    aff::GenericAffExpr{C,V}
end

GenericNormExpr(norm::GenericNorm{P,C,V}) where {P,C,V} =
    GenericNormExpr{C,V}(norm, one(C), zero(GenericAffExpr{C,V}))
Base.copy(x::GenericNormExpr{P,C,V}) where {P,C,V} =
    GenericNormExpr{P,C,V}(copy(x.norm), copy(x.coeff), copy(x.aff))
Base.convert(::Type{GenericNormExpr{P,C,V}}, x::GenericNorm{P,C,V}) where {P,C,V} =
    GenericNormExpr{P,C,V}(x, one(C), zero(GenericAffExpr{C,V}))

# Alias for ‖Ax‖₂ case
GenericSOCExpr{C,V} = GenericNormExpr{2,C,V}

# Alias for ‖Ax‖₂ and AffExpr case
const SOCExpr = GenericSOCExpr{Float64,Variable}

getvalue(n::GenericNormExpr) = n.coeff * getvalue(n.norm) + getvalue(n.aff)

##########################################################################
# GenericSOCConstraint
# Second-order cone constraint of form
# α||Ax-b||₂ + cᵀx + d ≤ 0
mutable struct GenericSOCConstraint{T<:GenericSOCExpr} <: AbstractConstraint
    normexpr::T
    function GenericSOCConstraint{T}(normexpr::T) where T
        if normexpr.coeff < 0
            # The coefficient in front of the norm is negative, which
            # means we have `norm >= c`, which is not convex.
            error("Invalid second-order cone constraint $(normexpr) ≤ 0")
        end
        new{T}(normexpr)
    end
end

# Alias for the AffExpr case
const SOCConstraint = GenericSOCConstraint{SOCExpr}

function addconstraint(m::Model, c::SOCConstraint)
    push!(m.socconstr,c)
    m.internalModelLoaded = false
    ConstraintRef{Model,SOCConstraint}(m,length(m.socconstr))
end

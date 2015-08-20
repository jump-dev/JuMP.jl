#  Copyright 2015, Iain Dunning, Joey Huchette, Miles Lubin, and contributors
#  This Source Code Form is subject to the terms of the Mozilla Public
#  License, v. 2.0. If a copy of the MPL was not distributed with this
#  file, You can obtain one at http://mozilla.org/MPL/2.0/.
#############################################################################
# JuMP
# An algebraic modelling langauge for Julia
# See http://github.com/JuliaOpt/JuMP.jl
#############################################################################
# src/sos.jl
# Defines SOSConstraint for "special ordered sets", and related methods
# SOS1: at most one of the variables in the set is positive, the rest are 0
# SOS2: at most two can be nonzero, if they a consecutive in the set
#############################################################################


type SOSConstraint <: JuMPConstraint
    terms::Vector{Variable}
    weights::Vector{Float64}
    sostype::Symbol
end
Base.copy(sos::SOSConstraint, new_model::Model) =
    SOSConstraint(Variable[copy(v, new_model) for v in sos.terms],
                    copy(sos.weights), sos.sostype)

# Given a vector of affine expressions, extract a vector of the single
# variable in each expression and a vector of their coefficients
function constructSOS(m::Model, coll::Vector{AffExpr})
    nvar = length(coll)
    vars = Array(Variable, nvar)
    weight = Array(Float64, nvar)
    for (i,aff) in enumerate(coll)
        if (length(aff.vars) != 1) || (aff.constant != 0)
            error("Must specify set in terms of single variables")
        end
        vars[i] = aff.vars[1]
        vars[i].m === m || error("Variable in constraint is not owned by the model")
        weight[i] = aff.coeffs[1]
    end
    return vars, weight
end


addSOS1(m::Model, coll) = addSOS1(m, convert(Vector{AffExpr}, coll))

function addSOS1(m::Model, coll::Vector{AffExpr})
    vars, weight = constructSOS(m,coll)
    push!(m.sosconstr, SOSConstraint(vars, weight, :SOS1))
    if m.internalModelLoaded
        idx = Int[v.col for v in vars]
        if applicable(MathProgBase.addsos1!, m.internalModel, idx, weight)
            MathProgBase.addsos1!(m.internalModel, idx, weight)
        else
            Base.warn_once("Solver does not appear to support adding constraints to an existing model. Hot-start is disabled.")
            m.internalModelLoaded = false
        end
    end
    return ConstraintRef{SOSConstraint}(m,length(m.sosconstr))
end

addSOS2(m::Model, coll) = addSOS2(m, convert(Vector{AffExpr}, coll))

function addSOS2(m::Model, coll::Vector{AffExpr})
    vars, weight = constructSOS(m,coll)
    push!(m.sosconstr, SOSConstraint(vars, weight, :SOS2))
    if m.internalModelLoaded
        idx = Int[v.col for v in vars]
        if applicable(MathProgBase.addsos2!, m.internalModel, idx, weight)
            MathProgBase.addsos2!(m.internalModel, idx, weight)
        else
            Base.warn_once("Solver does not appear to support adding constraints to an existing model. Hot-start is disabled.")
            m.internalModelLoaded = false
        end
    end
    return ConstraintRef{SOSConstraint}(m,length(m.sosconstr))
end
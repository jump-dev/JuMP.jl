#  Copyright 2017, Iain Dunning, Joey Huchette, Miles Lubin, and contributors
#  This Source Code Form is subject to the terms of the Mozilla Public
#  License, v. 2.0. If a copy of the MPL was not distributed with this
#  file, You can obtain one at http://mozilla.org/MPL/2.0/.
#############################################################################
# JuMP
# An algebraic modeling language for Julia
# See http://github.com/JuliaOpt/JuMP.jl
#############################################################################

# This file contains an example bridge used for tests.

using JuMP
const MOIB = MOI.Bridges
const MOIBC = MOI.Bridges.Constraint

"""
    Nonnegative <: MOI.AbstractScalarSet

Scalar set of nonnegative numbers.
"""
struct Nonnegative <: MOI.AbstractScalarSet end

"""
    NonnegativeBridge{T}

The `NonnegativeBridge` replaces a constraint `func`-in-`Nonnegative` into
`func`-in-`GreaterThan{T}`.
"""
struct NonnegativeBridge{T, F<:MOI.AbstractScalarFunction} <: MOIBC.AbstractBridge
    constraint_index::MOI.ConstraintIndex{F, MOI.GreaterThan{T}}
end

function MOIBC.bridge_constraint(::Type{NonnegativeBridge{T, F}},
                                 model,
                                 f::F,
                                 s::Nonnegative) where {T, F}
    ci = MOIU.normalize_and_add_constraint(model, f, MOI.GreaterThan(zero(T)))
    return NonnegativeBridge{T, F}(ci)
end

function MOI.supports_constraint(::Type{NonnegativeBridge{T}},
                                 ::Type{<:MOI.AbstractScalarFunction},
                                 ::Type{Nonnegative}) where T
    return true
end

function MOIB.added_constrained_variable_types(::Type{<:NonnegativeBridge})
    return Tuple{DataType}[]
end
function MOIB.added_constraint_types(::Type{NonnegativeBridge{T, F}}) where {T, F}
    return [(F, MOI.GreaterThan{T})]
end
function MOIBC.concrete_bridge_type(::Type{NonnegativeBridge{T}},
                                    F::Type{<:MOI.AbstractScalarFunction},
                                    ::Type{Nonnegative}) where T
    # In the constructor, the function `f` of type `F` is passed to
    # `MOIU.normalize_and_add_constraint` which removes the constrant from `f` but
    # does not change its type so the type of the function in `MOI.GreaterThan`
    # will also be `F`.
    return NonnegativeBridge{T, F}
end

# Attributes, Bridge acting as an model
function MOI.get(::NonnegativeBridge{T, F},
                 ::MOI.NumberOfConstraints{F, MOI.GreaterThan{T}}) where {T, F}
    return 1
end
function MOI.get(bridge::NonnegativeBridge{T, F},
                 ::MOI.ListOfConstraintIndices{F, MOI.GreaterThan{T}}) where {T, F}
    return [bridge.constraint_index]
end

# Indices
function MOI.delete(model::MOI.ModelLike, bridge::NonnegativeBridge)
    MOI.delete(model, bridge.constraint_index)
end

# Attributes, Bridge acting as a constraint
function MOI.get(model::MOI.ModelLike,
                 attr::Union{MOI.ConstraintPrimal, MOI.ConstraintDual},
                 bridge::NonnegativeBridge)
    return MOI.get(model, attr, bridge.constraint_index)
end

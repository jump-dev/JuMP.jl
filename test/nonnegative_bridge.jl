#  Copyright 2017, Iain Dunning, Joey Huchette, Miles Lubin, and contributors
#  This Source Code Form is subject to the terms of the Mozilla Public
#  License, v. 2.0. If a copy of the MPL was not distributed with this
#  file, You can obtain one at https://mozilla.org/MPL/2.0/.
#############################################################################
# JuMP
# An algebraic modeling language for Julia
# See https://github.com/jump-dev/JuMP.jl
#############################################################################

"""
    Nonnegative <: MOI.AbstractScalarSet

Scalar set of nonnegative numbers.
"""
struct Nonnegative <: MOI.AbstractScalarSet end

MOI.copy(set::Nonnegative) = set

struct BridgeMe{T,S}
    set::S
end

function JuMP.build_constraint(
    error::Function,
    f,
    set::BridgeMe{T,Nonnegative},
) where {T}
    return BridgeableConstraint(
        JuMP.build_constraint(error, f, set.set),
        NonnegativeBridge;
        coefficient_type = T,
    )
end

"""
    NonnegativeBridge{T}

The `NonnegativeBridge` replaces a constraint `func`-in-`Nonnegative` into
`func`-in-`GreaterThan{T}`.
"""
struct NonnegativeBridge{T,F<:MOI.AbstractScalarFunction} <:
       MOI.Bridges.Constraint.AbstractBridge
    constraint_index::MOI.ConstraintIndex{F,MOI.GreaterThan{T}}
end

function MOI.Bridges.Constraint.bridge_constraint(
    ::Type{NonnegativeBridge{T,F}},
    model::MOI.ModelLike,
    f::F,
    s::Nonnegative,
) where {T,F}
    ci = MOI.Utilities.normalize_and_add_constraint(
        model,
        f,
        MOI.GreaterThan(zero(T)),
    )
    return NonnegativeBridge{T,F}(ci)
end

function MOI.supports_constraint(
    ::Type{NonnegativeBridge{T}},
    ::Type{<:MOI.AbstractScalarFunction},
    ::Type{Nonnegative},
) where {T}
    return true
end

function MOI.Bridges.added_constrained_variable_types(
    ::Type{<:NonnegativeBridge},
)
    return Tuple{Type}[]
end

function MOI.Bridges.added_constraint_types(
    ::Type{NonnegativeBridge{T,F}},
) where {T,F}
    return [(F, MOI.GreaterThan{T})]
end

function MOI.Bridges.Constraint.concrete_bridge_type(
    ::Type{NonnegativeBridge{T}},
    F::Type{<:MOI.AbstractScalarFunction},
    ::Type{Nonnegative},
) where {T}
    # In the constructor, the function `f` of type `F` is passed to
    # `MOI.Utilities.normalize_and_add_constraint` which removes the constraint
    # from `f` but does not change its type so the type of the function in
    # `MOI.GreaterThan` will also be `F`.
    return NonnegativeBridge{T,F}
end

function MOI.get(
    ::NonnegativeBridge{T,F},
    ::MOI.NumberOfConstraints{F,MOI.GreaterThan{T}},
)::Int64 where {T,F}
    return 1
end

function MOI.get(
    bridge::NonnegativeBridge{T,F},
    ::MOI.ListOfConstraintIndices{F,MOI.GreaterThan{T}},
) where {T,F}
    return [bridge.constraint_index]
end

function MOI.delete(model::MOI.ModelLike, bridge::NonnegativeBridge)
    MOI.delete(model, bridge.constraint_index)
    return
end

function MOI.get(
    model::MOI.ModelLike,
    attr::Union{MOI.ConstraintPrimal,MOI.ConstraintDual},
    bridge::NonnegativeBridge,
)
    return MOI.get(model, attr, bridge.constraint_index)
end

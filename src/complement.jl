#  Copyright 2017, Iain Dunning, Joey Huchette, Miles Lubin, and contributors
#  This Source Code Form is subject to the terms of the Mozilla Public
#  License, v. 2.0. If a copy of the MPL was not distributed with this
#  file, You can obtain one at https://mozilla.org/MPL/2.0/.

# This file extends JuMP to complementarity constraints. It is a good example of
# how JuMP can be extended. In particular, it allows constraints of the form:
#
# @constraint(model, complements(2x - 1, x))
# @constraint(model, 2x - 1 ⟂ x)

function _build_complements_constraint(
    errorf::Function,
    F::AbstractArray{<:Union{Real,AbstractJuMPScalar}},
    x::AbstractArray{<:AbstractVariableRef},
)
    if size(F) != size(x)
        errorf(
            "size of mapping does not match size of variables: " *
            "$(size(F)) != $(size(x)).",
        )
    end
    return VectorConstraint([F; x], MOI.Complements(2 * length(F)))
end

function _build_complements_constraint(
    errorf::Function,
    F::Containers.SparseAxisArray{<:Union{Real,AbstractJuMPScalar}},
    x::Containers.SparseAxisArray{<:AbstractVariableRef},
)
    elements = [F[i] for i in eachindex(F)]
    for i in eachindex(F)
        if haskey(x, i)
            push!(elements, x[i])
        else
            errorf("keys of the SparseAxisArray's do not match.")
        end
    end
    return VectorConstraint(elements, MOI.Complements(length(elements)))
end

function _build_complements_constraint(
    errorf::Function,
    ::AbstractArray{<:Union{Real,AbstractJuMPScalar}},
    ::AbstractArray{<:AbstractJuMPScalar},
)
    return errorf("second term must be an array of variables.")
end

function _build_complements_constraint(
    ::Function,
    F::Union{Real,AbstractJuMPScalar},
    x::AbstractVariableRef,
)
    return VectorConstraint([F, x], MOI.Complements(2))
end

function _build_complements_constraint(
    errorf::Function,
    ::Union{Real,AbstractJuMPScalar},
    ::AbstractJuMPScalar,
)
    return errorf("second term must be a variable.")
end

function parse_constraint_call(
    errorf::Function,
    ::Bool,
    ::Union{Val{:complements},Val{:⟂}},
    F,
    x,
)
    f, parse_code = _rewrite_expression(F)
    return parse_code, :(_build_complements_constraint($errorf, $f, $(esc(x))))
end

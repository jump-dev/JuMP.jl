#  Copyright 2017, Iain Dunning, Joey Huchette, Miles Lubin, and contributors
#  This Source Code Form is subject to the terms of the Mozilla Public
#  License, v. 2.0. If a copy of the MPL was not distributed with this
#  file, You can obtain one at http://mozilla.org/MPL/2.0/.

# This file extends JuMP to complementarity constraints. It is a good example of
# how JuMP can be extended. In particular, it allows constraints of the form:
#
# @constraint(model, complements(2x - 1, x))
# @constraint(model, 2x - 1 ⟂ x)

function _build_complements_constraint(
    _error::Function,
    F::Vector{<:AbstractJuMPScalar},
    x::Vector{<:AbstractVariableRef}
)
    n = length(x)
    if length(F) != n
        _error(
            "Length of mapping ($(length(F))) is not equal to number of " *
            "matching variables ($(n))."
        )
    end
    return VectorConstraint([F; x], MOI.Complements(n))
end

function _build_complements_constraint(
    _error::Function,
    ::Vector{<:AbstractJuMPScalar},
    ::Vector{<:AbstractJuMPScalar},
)
    _error("second term must be a vector of variables.")
end

function _build_complements_constraint(
    ::Function, F::AbstractJuMPScalar, x::AbstractVariableRef
)
    return VectorConstraint([F, x], MOI.Complements(1))
end

function _build_complements_constraint(
    _error::Function, ::AbstractJuMPScalar, ::AbstractJuMPScalar
)
    _error("second term must be a variable.")
end

function parse_one_operator_constraint(
    _error::Function,
    ::Bool,
    ::Union{Val{:complements}, Val{:⟂}},
    F,
    x
)
    f, parse_code = _MA.rewrite(F)
    return parse_code, :(_build_complements_constraint($_error, $f, $(esc(x))))
end

#  Copyright 2017, Iain Dunning, Joey Huchette, Miles Lubin, and contributors
#  This Source Code Form is subject to the terms of the Mozilla Public
#  License, v. 2.0. If a copy of the MPL was not distributed with this
#  file, You can obtain one at http://mozilla.org/MPL/2.0/.

# This file extends JuMP to indicator constraints. It is a good example of how
# JuMP can be extended.

function _build_complements_constraint(
    _error::Function, variables::Vector{<:AbstractVariableRef},
    funcs::Vector{<:AbstractJuMPScalar})
    n = length(variables)
    if length(funcs) != n
        _error("Number of variables ($n) does not match number of constraints ($(length(funcs))).")
    end
    return VectorConstraint([variables; funcs], MOI.Complements(n))
end
function _build_complements_constraint(
    _error::Function, variable::AbstractVariableRef,
    func::AbstractJuMPScalar)
    return VectorConstraint([variable, func], MOI.Complements(1))
end

function _is_inequality(expr)
    return isexpr(expr, :call) && length(expr.args) == 3 && expr.args[1] in [:<=, :≤, :>=, :≥, :(.<=), :.≤, :(.>=), :.≥]
end
function parse_one_operator_constraint(_error::Function, _::Bool,
                                       s::Union{Val{:complements}, Val{:⟂}},
                                       lhs, rhs)
    if _is_inequality(lhs) && !_is_inequality(rhs)
        return parse_one_operator_constraint(_error, vectorized, s, rhs, lhs)
    end
    if !_is_inequality(rhs)
        _error("Expected one of the two sides of the complements to be an inequality.")
    end
    sense, vectorized = _check_vectorized(rhs.args[1])
    if sense in [:<=, :≤]
        l, r = rhs.args[3], rhs.args[2]
    else
        l, r = rhs.args[2], rhs.args[3]
    end
    sense = vectorized ? :.≥ : :≥
    _, rhs_parsecode, _build_call = parse_constraint(_error, sense, l, r)
    @assert length(_build_call.args) == 4
    @assert _build_call.args[1] == :build_constraint
    func = _build_call.args[3]
    return rhs_parsecode, :(_build_complements_constraint($_error, $(esc(lhs)), $func))
end

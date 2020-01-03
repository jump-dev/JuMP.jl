#  Copyright 2017, Iain Dunning, Joey Huchette, Miles Lubin, and contributors
#  This Source Code Form is subject to the terms of the Mozilla Public
#  License, v. 2.0. If a copy of the MPL was not distributed with this
#  file, You can obtain one at http://mozilla.org/MPL/2.0/.

# This file extends JuMP to indicator constraints. It is a good example of how
# JuMP can be extended.

function _build_complements_constraint(
    _error::Function, funcs::Vector{<:AbstractJuMPScalar},
    variables::Vector{<:AbstractVariableRef})
    n = length(variables)
    if length(funcs) != n
        _error("Number of variables ($n) does not match number of constraints ($(length(funcs))).")
    end
    return VectorConstraint([variables; funcs], MOI.Complements(n))
end
function _build_complements_constraint(
    _error::Function, func::AbstractJuMPScalar,
    variable::AbstractVariableRef)
    return VectorConstraint([variable, func], MOI.Complements(1))
end

function parse_one_operator_constraint(_error::Function, _::Bool,
                                       s::Union{Val{:complements}, Val{:⟂}},
                                       func_expr, variable)
    func, parse_code = _MA.rewrite(func_expr)
    return parse_code, :(_build_complements_constraint($_error, $func, $(esc(variable))))
end

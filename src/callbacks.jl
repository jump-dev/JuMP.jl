#  Copyright 2017, Iain Dunning, Joey Huchette, Miles Lubin, and contributors
#  This Source Code Form is subject to the terms of the Mozilla Public
#  License, v. 2.0. If a copy of the MPL was not distributed with this
#  file, You can obtain one at http://mozilla.org/MPL/2.0/.
#############################################################################
# JuMP
# An algebraic modeling language for Julia
# See http://github.com/JuliaOpt/JuMP.jl
#############################################################################

"""
    callback_value(cb_data, x::VariableRef)

Return the primal solution of a variable inside a callback.

`cb_data` is the argument to the callback function, and the type is dependent on
the solver.
"""
function callback_value(cb_data, x::VariableRef)
    # TODO(odow):
    # MOI defines `is_set_by_optimize(::CallbackVariablePrimal) = true`.
    # This causes problems for JuMP because it checks the termination_status to
    # see if optimize! has been called. Solutions are:
    # 1) defining is_set_by_optimize = false
    # 2) adding a flag to JuMP to store whether it is in a callback
    # 3) adding IN_OPTIMIZE to termination_status for callbacks
    # Once this is resolved, we can replace the current function with:
    #     MOI.get(owner_model(x), MOI.CallbackVariablePrimal(cb_data), x)
    return MOI.get(
        backend(owner_model(x)), MOI.CallbackVariablePrimal(cb_data), index(x)
    )
end

"""
    callback_value(cb_data, expr::GenericAffExpr)

Return the primal solution of an affine expression inside a callback by getting
the value for each variable appearing in the expression.

`cb_data` is the argument to the callback function, and the type is dependent on
the solver.
"""
function callback_value(cb_data, expr::GenericAffExpr)
    return expr.constant + sum(callback_value(cb_data, var) * coeff for (var, coeff) in expr.terms)
end

"""
    callback_value(cb_data, expr::GenericQuadExpr)

Return the primal solution of a quadratic expression inside a callback by getting
the value for each variable appearing in the expression.

`cb_data` is the argument to the callback function, and the type is dependent on
the solver.
"""
function callback_value(cb_data, expr::GenericQuadExpr)
    return callback_value(cb_data, expr.aff) + sum(callback_value(cb_data, var.first) * callback_value(cb_data, var.second) * coeff for (var, coeff) in expr.terms)
end

function MOI.submit(model::Model, cb::MOI.LazyConstraint, con::ScalarConstraint)
    return MOI.submit(backend(model), cb, moi_function(con.func), con.set)
end

function MOI.submit(model::Model, cb::MOI.UserCut, con::ScalarConstraint)
    return MOI.submit(backend(model), cb, moi_function(con.func), con.set)
end

function MOI.submit(
    model::Model,
    cb::MOI.HeuristicSolution,
    variables::Vector{VariableRef},
    values::Vector{Float64}
)
    return MOI.submit(backend(model), cb, index.(variables), values)
end

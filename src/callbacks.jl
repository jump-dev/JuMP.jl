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
    callback_node_status(cb_data, model::GenericModel)

Return an [`MOI.CallbackNodeStatusCode`](@ref) enum, indicating if the current
primal solution available from [`callback_value`](@ref) is integer feasible.
"""
function callback_node_status(cb_data, model::GenericModel)
    # TODO(odow):
    # MOI defines `is_set_by_optimize(::CallbackNodeStatus) = true`.
    # This causes problems for JuMP because it checks the termination_status to
    # see if optimize! has been called. Solutions are:
    # 1) defining is_set_by_optimize = false
    # 2) adding a flag to JuMP to store whether it is in a callback
    # 3) adding IN_OPTIMIZE to termination_status for callbacks
    # Once this is resolved, we can replace the current function with:
    #     MOI.get(model, MOI.CallbackNodeStatus(cb_data))
    return MOI.get(backend(model), MOI.CallbackNodeStatus(cb_data))
end

"""
    callback_value(cb_data, x::GenericVariableRef)

Return the primal solution of a variable inside a callback.

`cb_data` is the argument to the callback function, and the type is dependent on
the solver.
"""
function callback_value(cb_data, x::GenericVariableRef)
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
        backend(owner_model(x)),
        MOI.CallbackVariablePrimal(cb_data),
        index(x),
    )
end

"""
    callback_value(cb_data, expr::Union{GenericAffExpr, GenericQuadExpr})

Return the primal solution of an affine or quadratic expression inside a
callback by getting the value for each variable appearing in the expression.

`cb_data` is the argument to the callback function, and the type is dependent on
the solver.
"""
function callback_value(cb_data, expr::Union{GenericAffExpr,GenericQuadExpr})
    return value(expr) do x
        return callback_value(cb_data, x)
    end
end

function MOI.submit(
    model::GenericModel,
    cb::MOI.LazyConstraint,
    con::ScalarConstraint,
)
    return MOI.submit(backend(model), cb, moi_function(con.func), con.set)
end

function MOI.submit(model::GenericModel, cb::MOI.UserCut, con::ScalarConstraint)
    return MOI.submit(backend(model), cb, moi_function(con.func), con.set)
end

function MOI.submit(
    model::GenericModel{T},
    cb::MOI.HeuristicSolution,
    variables::Vector{GenericVariableRef{T}},
    values::Vector{<:Real},
) where {T}
    return MOI.submit(
        backend(model),
        cb,
        index.(variables),
        convert(Vector{T}, values),
    )
end

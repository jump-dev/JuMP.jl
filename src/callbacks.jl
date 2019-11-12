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
    return MOI.get(
        backend(owner_model(x)),
        MOI.CallbackVariablePrimal(cb_data),
        index(x)
    )
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

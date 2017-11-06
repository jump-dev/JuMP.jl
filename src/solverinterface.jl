#  Copyright 2017, Iain Dunning, Joey Huchette, Miles Lubin, and contributors
#  This Source Code Form is subject to the terms of the Mozilla Public
#  License, v. 2.0. If a copy of the MPL was not distributed with this
#  file, You can obtain one at http://mozilla.org/MPL/2.0/.

"""
    attach(m::JuMP.Model)

Attach the JuMP model `m` to a new solver instance. All instance data is transferred to the solver instance,
and while the solver instance remains attached, all modifications to the JuMP model are immediately mirrored
in the solver instance. See also `isattached()` and `detach()`.
"""
function attach(m::Model)
    @assert !m.solverinstanceattached
    m.solverinstance = MOI.SolverInstance(m.solver)
    solvervariables = MOI.addvariables!(m.solverinstance, numvar(m)) # TODO numvar shouldn't return unsigned int

    m.variabletosolvervariable = Dict{MOIVAR,MOIVAR}()
    # TODO: replace with ListOfVariableReferences()
    # Now we're assuming all instance variables are numbered sequentially
    for i in 1:numvar(m)
        m.variabletosolvervariable[MOIVAR(i)] = solvervariables[i]
    end

    MOI.set!(m.solverinstance, MOI.ObjectiveSense(), MOI.get(m.instance, MOI.ObjectiveSense()))
    MOI.set!(m.solverinstance, MOI.ObjectiveFunction(), MOI.get(m.instance, MOI.ObjectiveFunction()))

    # TODO: replace with ListOfConstraintReferences()
    # This would be a bit more transparent than the call below
    # TODO: keep track of constraint references!
    MOIU.broadcastcall( constrs -> for (cref, f, s) in constrs; MOI.addconstraint!(m.solverinstance, f, s) end, m.instance)

    m.solverinstanceattached = true
    nothing

end

"""
    isattached(m::JuMP.Model)::Bool

Return `true` if the JuMP model `m` is currently attached to a solver instance,
`false` otherwise.
"""
function isattached(m::Model)
    return m.solverinstanceattached
end

"""
    detach(m::JuMP.Model)

Detach the JuMP model `m` from the attached solver instance, after calling
`MathOptInterface.free!` on it to release any resources in use.
Throws an error if no solver instance is currently attached.
"""
function detach(m::Model)
    @assert m.solverinstanceattached
    m.solverinstanceattached = false
    MOI.free!(m.solverinstance)
    m.solverinstance = nothing
end

"""
    canget(m::JuMP.Model, attr::MathOptInterface.AbstractSolverInstanceAttribute)::Bool

Return `true` if one may query the attribute `attr` from the solver instance attached to the JuMP model,
false if not.
Throws an error if no solver instance is currently attached.
"""
function MOI.canget(m::Model, attr::MOI.AbstractSolverInstanceAttribute)
    @assert m.solverinstanceattached
    return MOI.canget(m.solverinstance, attr)
end

"""
    get(m::JuMP.Model, attr::MathOptInterface.AbstractSolverInstanceAttribute)

Return the value of the attribute `attr` from the solver instance attached to the JuMP model.
Throws an error if no solver instance is currently attached.
"""
function MOI.get(m::Model, attr::MOI.AbstractSolverInstanceAttribute)
    @assert m.solverinstanceattached
    return MOI.get(m.solverinstance, attr)
end

function solve(m::Model;
                ignore_solve_hook=(m.solvehook===nothing))

    # If the user or an extension has provided a solve hook, call
    # that instead of solving the model ourselves
    if !ignore_solve_hook
        return m.solvehook(m)
    end

    if !m.solverinstanceattached
        attach(m)
    end
    MOI.optimize!(m.solverinstance)

    empty!(m.variableresult)
    # If any variable has a result then all must have
    if MOI.canget(m.solverinstance, MOI.VariablePrimal(), first(m.variabletosolvervariable).second)
        for vref in keys(m.variabletosolvervariable)
            m.variableresult[Variable(m,vref)] = MOI.get(m.solverinstance, MOI.VariablePrimal(), m.variabletosolvervariable[vref])
        end
    end

    nothing

end

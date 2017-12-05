#  Copyright 2017, Iain Dunning, Joey Huchette, Miles Lubin, and contributors
#  This Source Code Form is subject to the terms of the Mozilla Public
#  License, v. 2.0. If a copy of the MPL was not distributed with this
#  file, You can obtain one at http://mozilla.org/MPL/2.0/.


"""
    copyconstraints!(m::JuMP.Model, ::Type{F}, ::Type{S}) where {F<:MOI.AbstractFunction, S<:MOI.AbstractSet}

Transfer the constraints of type `F`-in-`S` to the solver instance.
"""
function copyconstraints!(m::Model, ::Type{F}, ::Type{S}) where {F<:MOI.AbstractFunction, S<:MOI.AbstractSet}
    for cindex in MOI.get(m.instance, MOI.ListOfConstraintIndices{F, S}())
        f = MOI.get(m.instance, MOI.ConstraintFunction(), cindex)
        s = MOI.get(m.instance, MOI.ConstraintSet(), cindex)
        solvercindex = MOI.addconstraint!(m.solverinstance, f, s)
        @assert !haskey(m.constrainttosolverconstraint, cindex)
        m.constrainttosolverconstraint[cindex] = solvercindex
    end
end

"""
    attach(m::JuMP.Model)

Attach the JuMP model `m` to a new solver instance. All instance data is transferred to the solver instance,
and while the solver instance remains attached, all modifications to the JuMP model are immediately mirrored
in the solver instance. See also `isattached()` and `detach()`.
"""
function attach(m::Model, solverinstance::MOI.AbstractSolverInstance)
    @assert !m.solverinstanceattached
    m.solverinstance = solverinstance
    solvervariables = MOI.addvariables!(m.solverinstance, numvar(m)) # TODO numvar shouldn't return unsigned int

    m.variabletosolvervariable = Dict{MOIVAR,MOIVAR}()
    m.constrainttosolverconstraint = Dict{UInt64,UInt64}()
    # TODO: replace with ListOfVariableIndices()
    # Now we're assuming all instance variables are numbered sequentially
    for i in 1:numvar(m)
        m.variabletosolvervariable[MOIVAR(i)] = solvervariables[i]
    end

    MOI.set!(m.solverinstance, MOI.ObjectiveSense(), MOI.get(m.instance, MOI.ObjectiveSense()))
    MOI.set!(m.solverinstance, MOI.ObjectiveFunction(), MOI.get(m.instance, MOI.ObjectiveFunction()))

    for (F, S) in MOI.get(m.instance, MOI.ListOfConstraints())
        # do the rest in copyconstraints! which is type stable
        copyconstraints!(m, F, S)
    end

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
    canget(m::JuMP.Model, attr::MathOptInterface.AbstractInstanceAttribute)::Bool

Return `true` if one may query the attribute `attr` from the solver instance attached to the JuMP model,
false if not.
Throws an error if no solver instance is currently attached.
"""
function MOI.canget(m::Model, attr::MOI.AbstractInstanceAttribute)
    @assert m.solverinstanceattached
    return MOI.canget(m.solverinstance, attr)
end

"""
    get(m::JuMP.Model, attr::MathOptInterface.AbstractInstanceAttribute)

Return the value of the attribute `attr` from the solver instance attached to the JuMP model.
Throws an error if no solver instance is currently attached.
"""
function MOI.get(m::Model, attr::MOI.AbstractInstanceAttribute)
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

    @assert m.solverinstanceattached
    MOI.optimize!(m.solverinstance)

    empty!(m.variableresult)
    # If any variable has a result then all must have
    if MOI.canget(m.solverinstance, MOI.VariablePrimal(), first(m.variabletosolvervariable).second)
        for vindex in keys(m.variabletosolvervariable)
            m.variableresult[Variable(m,vindex)] = MOI.get(m.solverinstance, MOI.VariablePrimal(), m.variabletosolvervariable[vindex])
        end
    end

    nothing

end

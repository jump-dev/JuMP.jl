#  Copyright 2017, Iain Dunning, Joey Huchette, Miles Lubin, and contributors
#  This Source Code Form is subject to the terms of the Mozilla Public
#  License, v. 2.0. If a copy of the MPL was not distributed with this
#  file, You can obtain one at http://mozilla.org/MPL/2.0/.

# These methods directly map to InstanceManager methods.
# They cannot be called in Direct mode.
function MOIU.resetsolver!(m::Model, solver::MOI.AbstractSolverInstance)
    @assert mode(m) != Direct
    MOIU.resetsolver!(m.moibackend, solver)
end

function MOIU.resetsolver!(m::Model)
    @assert mode(m) != Direct
    MOIU.resetsolver!(m.moibackend)
end

function MOIU.dropsolver!(m::Model)
    @assert mode(m) != Direct
    MOIU.dropsolver!(m.moibackend)
end

function MOIU.attachsolver!(m::Model)
    @assert mode(m) != Direct
    copyresult = MOIU.attachsolver!(m.moibackend)
    # TODO: more reliable error reporting
    @assert copyresult.status == MOI.CopySuccess
    @assert m.moibackend.state == MOIU.AttachedSolver
    return copyresult
end


function solve(m::Model;
                ignore_solve_hook=(m.solvehook===nothing))
    # If the user or an extension has provided a solve hook, call
    # that instead of solving the model ourselves
    if !ignore_solve_hook
        return m.solvehook(m)
    end

    MOI.optimize!(m.moibackend)

    return
end

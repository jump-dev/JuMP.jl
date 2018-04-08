#  Copyright 2017, Iain Dunning, Joey Huchette, Miles Lubin, and contributors
#  This Source Code Form is subject to the terms of the Mozilla Public
#  License, v. 2.0. If a copy of the MPL was not distributed with this
#  file, You can obtain one at http://mozilla.org/MPL/2.0/.

# These methods directly map to CachingOptimizer methods.
# They cannot be called in Direct mode.
function MOIU.resetoptimizer!(m::Model, optimizer::MOI.AbstractOptimizer)
    @assert mode(m) != Direct
    MOIU.resetoptimizer!(m.moibackend, optimizer)
end

function MOIU.resetoptimizer!(m::Model)
    @assert mode(m) != Direct
    MOIU.resetoptimizer!(m.moibackend)
end

function MOIU.dropoptimizer!(m::Model)
    @assert mode(m) != Direct
    MOIU.dropoptimizer!(m.moibackend)
end

function MOIU.attachoptimizer!(m::Model)
    @assert mode(m) != Direct
    copyresult = MOIU.attachoptimizer!(m.moibackend)
    # TODO: more reliable error reporting
    @assert copyresult.status == MOI.CopySuccess
    @assert m.moibackend.state == MOIU.AttachedOptimizer
    return copyresult
end


function optimize(m::Model;
                ignore_optimize_hook=(m.optimizehook===nothing))
    # The NLPData is not kept in sync, so re-set it here.
    # TODO: Consider how to handle incremental solves.
    if m.nlpdata !== nothing
        MOI.set!(m, MOI.NLPBlock(), create_nlp_block_data(m))
    end

    # If the user or an extension has provided an optimize hook, call
    # that instead of solving the model ourselves
    if !ignore_optimize_hook
        return m.optimizehook(m)
    end

    MOI.optimize!(m.moibackend)

    return
end

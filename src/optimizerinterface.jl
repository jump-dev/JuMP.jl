#  Copyright 2017, Iain Dunning, Joey Huchette, Miles Lubin, and contributors
#  This Source Code Form is subject to the terms of the Mozilla Public
#  License, v. 2.0. If a copy of the MPL was not distributed with this
#  file, You can obtain one at http://mozilla.org/MPL/2.0/.

# These methods directly map to CachingOptimizer methods.
# They cannot be called in Direct mode.
function MOIU.resetoptimizer!(model::Model, optimizer::MOI.AbstractOptimizer)
    @assert mode(model) != Direct
    MOIU.resetoptimizer!(caching_optimizer(model), optimizer)
end

function MOIU.resetoptimizer!(model::Model)
    @assert mode(model) != Direct
    MOIU.resetoptimizer!(caching_optimizer(model))
end

function MOIU.dropoptimizer!(model::Model)
    @assert mode(model) != Direct
    MOIU.dropoptimizer!(caching_optimizer(model))
end

function MOIU.attachoptimizer!(model::Model)
    @assert mode(model) != Direct
    copyresult = MOIU.attachoptimizer!(caching_optimizer(model))
    # TODO: more reliable error reporting
    @assert copyresult.status == MOI.CopySuccess
    @assert caching_optimizer(model).state == MOIU.AttachedOptimizer
    return copyresult
end


function optimize(model::Model;
                ignore_optimize_hook=(model.optimizehook===nothing))
    # The NLPData is not kept in sync, so re-set it here.
    # TODO: Consider how to handle incremental solves.
    if model.nlpdata !== nothing
        MOI.set!(model, MOI.NLPBlock(), create_nlp_block_data(model))
        empty!(model.nlpdata.nlconstr_duals)
    end

    # If the user or an extension has provided an optimize hook, call
    # that instead of solving the model ourselves
    if !ignore_optimize_hook
        return model.optimizehook(model)
    end

    MOI.optimize!(model.moibackend)

    return
end

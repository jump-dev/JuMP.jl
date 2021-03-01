#  Copyright 2017, Iain Dunning, Joey Huchette, Miles Lubin, and contributors
#  This Source Code Form is subject to the terms of the Mozilla Public
#  License, v. 2.0. If a copy of the MPL was not distributed with this
#  file, You can obtain one at https://mozilla.org/MPL/2.0/.

"""
    error_if_direct_mode(model::Model, func::Symbol)

Errors if `model` is in direct mode during a call from the function named
`func`.

Used internally within JuMP, or by JuMP extensions who do not want to support
models in direct mode.
"""
function error_if_direct_mode(model::Model, func::Symbol)
    if mode(model) == DIRECT
        error("The `$func` function is not supported in DIRECT mode.")
    end
end

# These methods directly map to CachingOptimizer methods.

"""
    MOIU.reset_optimizer(model::Model, optimizer::MOI.AbstractOptimizer)

Call `MOIU.reset_optimizer` on the backend of `model`.

Cannot be called in direct mode.
"""
function MOIU.reset_optimizer(
    model::Model,
    optimizer::MOI.AbstractOptimizer,
    ::Bool = true,
)
    error_if_direct_mode(model, :reset_optimizer)
    return MOIU.reset_optimizer(backend(model), optimizer)
end

"""
    MOIU.reset_optimizer(model::Model)

Call `MOIU.reset_optimizer` on the backend of `model`.

Cannot be called in direct mode.
"""
function MOIU.reset_optimizer(model::Model)
    error_if_direct_mode(model, :reset_optimizer)
    return MOIU.reset_optimizer(backend(model))
end

"""
    MOIU.drop_optimizer(model::Model)

Call `MOIU.drop_optimizer` on the backend of `model`.

Cannot be called in direct mode.
"""
function MOIU.drop_optimizer(model::Model)
    error_if_direct_mode(model, :drop_optimizer)
    return MOIU.drop_optimizer(backend(model))
end

"""
    MOIU.attach_optimizer(model::Model)

Call `MOIU.attach_optimizer` on the backend of `model`.

Cannot be called in direct mode.
"""
function MOIU.attach_optimizer(model::Model)
    error_if_direct_mode(model, :attach_optimizer)
    return MOIU.attach_optimizer(backend(model))
end

"""
    set_optimizer(model::Model, optimizer_factory;
                  bridge_constraints::Bool=true)


Creates an empty `MathOptInterface.AbstractOptimizer` instance by calling
`optimizer_factory()` and sets it as the optimizer of `model`. Specifically,
`optimizer_factory` must be callable with zero arguments and return an empty
`MathOptInterface.AbstractOptimizer`.

If `bridge_constraints` is true, constraints that are not supported by the
optimizer are automatically bridged to equivalent supported constraints when
an appropriate transformation is defined in the `MathOptInterface.Bridges`
module or is defined in another module and is explicitly added.

See [`set_optimizer_attributes`](@ref) and [`set_optimizer_attribute`](@ref) for setting
solver-specific parameters of the optimizer.

## Examples
```julia
model = Model()
set_optimizer(model, GLPK.Optimizer)
```
"""
function set_optimizer(
    model::Model,
    optimizer_constructor;
    bridge_constraints::Bool = true,
)
    error_if_direct_mode(model, :set_optimizer)
    if bridge_constraints
        # We set `with_names=false` because the names are handled by the first
        # caching optimizer. If `default_copy_to` without names is supported,
        # no need for a second cache.
        optimizer = MOI.instantiate(
            optimizer_constructor,
            with_bridge_type = Float64,
            with_names = false,
        )
        for bridge_type in model.bridge_types
            _moi_add_bridge(optimizer, bridge_type)
        end
    else
        optimizer = MOI.instantiate(optimizer_constructor)
    end
    return MOIU.reset_optimizer(model, optimizer)
end

# Deprecation for JuMP v0.18 -> JuMP v0.19 transition
export solve
function solve(::Model)
    return error(
        "`solve` has been replaced by `optimize!`. Note that `solve` " *
        "used to return a `Symbol` summarizing the solution while " *
        "`optimize!` returns nothing and the status of the solution " *
        "is queried using `termination_status`, `primal_status` " *
        "and `dual_status`.",
    )
end

"""
    optimize!(model::Model;
              ignore_optimize_hook=(model.optimize_hook === nothing),
              kwargs...)

Optimize the model. If an optimizer has not been set yet (see
[`set_optimizer`](@ref)), a [`NoOptimizer`](@ref) error is thrown.

Keyword arguments `kwargs` are passed to the `optimize_hook`. An error is
thrown if `optimize_hook` is `nothing` and keyword arguments are provided.
"""
function optimize!(
    model::Model,
    # TODO: Remove the optimizer_factory and bridge_constraints
    # arguments when the deprecation error below is removed.
    optimizer_factory = nothing;
    bridge_constraints::Bool = true,
    ignore_optimize_hook = (model.optimize_hook === nothing),
    kwargs...,
)
    # The nlp_data is not kept in sync, so re-set it here.
    # TODO: Consider how to handle incremental solves.
    if model.nlp_data !== nothing
        MOI.set(model, MOI.NLPBlock(), _create_nlp_block_data(model))
        empty!(model.nlp_data.nlconstr_duals)
    end

    if optimizer_factory !== nothing
        # This argument was deprecated in JuMP 0.21.
        error(
            "The optimizer factory argument is no longer accepted by " *
            "`optimize!`. Call `set_optimizer` before `optimize!`.",
        )
    end

    # If the user or an extension has provided an optimize hook, call
    # that instead of solving the model ourselves
    if !ignore_optimize_hook
        return model.optimize_hook(model; kwargs...)
    end

    isempty(kwargs) || error(
        "Unrecognized keyword arguments: $(join([k[1] for k in kwargs], ", "))",
    )

    if mode(model) != DIRECT && MOIU.state(backend(model)) == MOIU.NO_OPTIMIZER
        throw(NoOptimizer())
    end

    try
        MOI.optimize!(backend(model))
    catch err
        # TODO: This error also be thrown also in MOI.set() if the solver is
        # attached. Currently we catch only the more common case. More generally
        # JuMP is missing a translation layer from MOI errors to JuMP errors.
        if err isa MOI.UnsupportedAttribute{MOI.NLPBlock}
            error(
                "The solver does not support nonlinear problems " *
                "(i.e., NLobjective and NLconstraint).",
            )
        else
            rethrow(err)
        end
    end

    return
end

"""
    compute_conflict!(model::Model)

Compute a conflict if the model is infeasible. If an optimizer has not
been set yet (see [`set_optimizer`](@ref)), a [`NoOptimizer`](@ref)
error is thrown.

The status of the conflict can be checked with the `MOI.ConflictStatus`
model attribute. Then, the status for each constraint can be queried with
the `MOI.ConstraintConflictStatus` attribute.
"""
function compute_conflict!(model::Model)
    if mode(model) != DIRECT && MOIU.state(backend(model)) == MOIU.NO_OPTIMIZER
        throw(NoOptimizer())
    end
    MOI.compute_conflict!(backend(model))
    return
end

"""
    result_count(model::Model)

Return the number of results available to query after a call to
[`optimize!`](@ref).
"""
function result_count(model::Model)::Int
    return MOI.get(model, MOI.ResultCount())
end

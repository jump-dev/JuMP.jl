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
    set_optimizer(
        model::Model,
        optimizer_factory;
        add_bridges::Bool = true,
    )

Creates an empty `MathOptInterface.AbstractOptimizer` instance by calling
`optimizer_factory()` and sets it as the optimizer of `model`. Specifically,
`optimizer_factory` must be callable with zero arguments and return an empty
`MathOptInterface.AbstractOptimizer`.

If `add_bridges` is true, constraints and objectives that are not supported by
the optimizer are automatically bridged to equivalent supported formulation.
Passing `add_bridges = false` can improve performance if the solver natively
supports all of the elements in `model`.

See [`set_optimizer_attributes`](@ref) and [`set_optimizer_attribute`](@ref) for
setting solver-specific parameters of the optimizer.

## Examples

```julia
model = Model()
set_optimizer(model, HiGHS.Optimizer)
set_optimizer(model, HiGHS.Optimizer; add_bridges = false)
```
"""
function set_optimizer(
    model::Model,
    optimizer_constructor;
    add_bridges::Bool = true,
)
    error_if_direct_mode(model, :set_optimizer)
    if add_bridges
        optimizer =
            MOI.instantiate(optimizer_constructor, with_bridge_type = Float64)
        for bridge_type in model.bridge_types
            _moi_add_bridge(optimizer, bridge_type)
        end
    else
        optimizer = MOI.instantiate(optimizer_constructor)
    end
    # Update the backend to create a new, concretely typed CachingOptimizer
    # using the existing `model_cache`.
    model.moi_backend =
        MOI.Utilities.CachingOptimizer(backend(model).model_cache, optimizer)
    return
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
    model::Model;
    ignore_optimize_hook = (model.optimize_hook === nothing),
    kwargs...,
)
    # The nlp_data is not kept in sync, so re-set it here.
    # TODO: Consider how to handle incremental solves.
    if model.nlp_data !== nothing
        block = MOI.NLPBlockData(model.nlp_data, index.(all_variables(model)))
        MOI.set(model, MOI.NLPBlock(), block)
    end
    # If the user or an extension has provided an optimize hook, call
    # that instead of solving the model ourselves
    if !ignore_optimize_hook
        return model.optimize_hook(model; kwargs...)
    end
    if !isempty(kwargs)
        error(
            "Unrecognized keyword arguments: $(join([k[1] for k in kwargs], ", "))",
        )
    end
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
    model.is_model_dirty = false
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
    if termination_status(model) == MOI.OPTIMIZE_NOT_CALLED
        return 0
    end
    return MOI.get(model, MOI.ResultCount())
end

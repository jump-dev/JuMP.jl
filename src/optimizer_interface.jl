#  Copyright 2017, Iain Dunning, Joey Huchette, Miles Lubin, and contributors
#  This Source Code Form is subject to the terms of the Mozilla Public
#  License, v. 2.0. If a copy of the MPL was not distributed with this
#  file, You can obtain one at https://mozilla.org/MPL/2.0/.

"""
    optimizer_with_attributes(optimizer_constructor, attrs::Pair...)

Groups an optimizer constructor with the list of attributes `attrs`. Note that
it is equivalent to `MOI.OptimizerWithAttributes`.

When provided to the `Model` constructor or to [`set_optimizer`](@ref), it
creates an optimizer by calling `optimizer_constructor()`, and then sets the
attributes using [`set_optimizer_attribute`](@ref).

## Example

```julia
model = Model(
    optimizer_with_attributes(
        Gurobi.Optimizer, "Presolve" => 0, "OutputFlag" => 1
    )
)
```
is equivalent to:
```julia
model = Model(Gurobi.Optimizer)
set_optimizer_attribute(model, "Presolve", 0)
set_optimizer_attribute(model, "OutputFlag", 1)
```

## Note

The string names of the attributes are specific to each solver. One should
consult the solver's documentation to find the attributes of interest.

See also: [`set_optimizer_attribute`](@ref), [`set_optimizer_attributes`](@ref),
[`get_optimizer_attribute`](@ref).
"""
function optimizer_with_attributes(optimizer_constructor, args::Pair...)
    return MOI.OptimizerWithAttributes(optimizer_constructor, args...)
end

"""
    set_optimizer_attribute(
        model::Union{Model,MOI.OptimizerWithAttributes},
        name::String,
        value,
    )

Sets solver-specific attribute identified by `name` to `value`.

Note that this is equivalent to
`set_optimizer_attribute(model, MOI.RawOptimizerAttribute(name), value)`.

## Example

```julia
set_optimizer_attribute(model, "SolverSpecificAttributeName", true)
```

See also: [`set_optimizer_attributes`](@ref), [`get_optimizer_attribute`](@ref).
"""
function set_optimizer_attribute(
    model::Union{Model,MOI.OptimizerWithAttributes},
    name::String,
    value,
)
    set_optimizer_attribute(model, MOI.RawOptimizerAttribute(name), value)
    return
end

# This method is needed for string types like String15 coming from a DataFrame.
function set_optimizer_attribute(
    model::Union{Model,MOI.OptimizerWithAttributes},
    name::AbstractString,
    value,
)
    set_optimizer_attribute(model, String(name), value)
    return
end

"""
    set_optimizer_attribute(
        model::Union{Model,MOI.OptimizerWithAttributes},
        attr::MOI.AbstractOptimizerAttribute,
        value,
    )

Set the solver-specific attribute `attr` in `model` to `value`.

## Example

```julia
set_optimizer_attribute(model, MOI.Silent(), true)
```

See also: [`set_optimizer_attributes`](@ref), [`get_optimizer_attribute`](@ref).
"""
function set_optimizer_attribute(
    model::Union{Model,MOI.OptimizerWithAttributes},
    attr::MOI.AbstractOptimizerAttribute,
    value,
)
    MOI.set(model, attr, value)
    return
end

"""
    set_optimizer_attributes(
        model::Union{Model,MOI.OptimizerWithAttributes},
        pairs::Pair...,
    )

Given a list of `attribute => value` pairs, calls
`set_optimizer_attribute(model, attribute, value)` for each pair.

## Example

```julia
model = Model(Ipopt.Optimizer)
set_optimizer_attributes(model, "tol" => 1e-4, "max_iter" => 100)
```
is equivalent to:
```julia
model = Model(Ipopt.Optimizer)
set_optimizer_attribute(model, "tol", 1e-4)
set_optimizer_attribute(model, "max_iter", 100)
```

See also: [`set_optimizer_attribute`](@ref), [`get_optimizer_attribute`](@ref).
"""
function set_optimizer_attributes(
    model::Union{Model,MOI.OptimizerWithAttributes},
    pairs::Pair...,
)
    for (name, value) in pairs
        set_optimizer_attribute(model, name, value)
    end
    return
end

"""
    get_optimizer_attribute(
        model::Union{Model,MOI.OptimizerWithAttributes},
        name::String,
    )

Return the value associated with the solver-specific attribute named `name`.

Note that this is equivalent to
`get_optimizer_attribute(model, MOI.RawOptimizerAttribute(name))`.

## Example

```julia
get_optimizer_attribute(model, "SolverSpecificAttributeName")
```

See also: [`set_optimizer_attribute`](@ref), [`set_optimizer_attributes`](@ref).
"""
function get_optimizer_attribute(
    model::Union{Model,MOI.OptimizerWithAttributes},
    name::String,
)
    return get_optimizer_attribute(model, MOI.RawOptimizerAttribute(name))
end

# This method is needed for string types like String15 coming from a DataFrame.
function get_optimizer_attribute(
    model::Union{Model,MOI.OptimizerWithAttributes},
    name::AbstractString,
)
    return get_optimizer_attribute(model, String(name))
end

"""
    get_optimizer_attribute(
        model::Union{Model,MOI.OptimizerWithAttributes},
        attr::MOI.AbstractOptimizerAttribute,
    )

Return the value of the solver-specific attribute `attr` in `model`.

## Example

```julia
get_optimizer_attribute(model, MOI.Silent())
```

See also: [`set_optimizer_attribute`](@ref), [`set_optimizer_attributes`](@ref).
"""
function get_optimizer_attribute(
    model::Union{Model,MOI.OptimizerWithAttributes},
    attr::MOI.AbstractOptimizerAttribute,
)
    return MOI.get(model, attr)
end

"""
    set_silent(model::Model)

Takes precedence over any other attribute controlling verbosity and requires the
solver to produce no output.

See also: [`unset_silent`](@ref).
"""
function set_silent(model::Model)
    return MOI.set(model, MOI.Silent(), true)
end

"""
    unset_silent(model::Model)

Neutralize the effect of the `set_silent` function and let the solver attributes
control the verbosity.

See also: [`set_silent`](@ref).
"""
function unset_silent(model::Model)
    return MOI.set(model, MOI.Silent(), false)
end

"""
    set_time_limit_sec(model::Model, limit::Float64)

Set the time limit (in seconds) of the solver.

Can be unset using [`unset_time_limit_sec`](@ref) or with `limit` set to
`nothing`.

See also: [`unset_time_limit_sec`](@ref), [`time_limit_sec`](@ref).
"""
function set_time_limit_sec(model::Model, limit::Real)
    return MOI.set(model, MOI.TimeLimitSec(), convert(Float64, limit))
end

function set_time_limit_sec(model::Model, ::Nothing)
    return unset_time_limit_sec(model)
end

"""
    unset_time_limit_sec(model::Model)

Unset the time limit of the solver.

See also: [`set_time_limit_sec`](@ref), [`time_limit_sec`](@ref).
"""
function unset_time_limit_sec(model::Model)
    return MOI.set(model, MOI.TimeLimitSec(), nothing)
end

"""
    time_limit_sec(model::Model)

Return the time limit (in seconds) of the `model`.

Returns `nothing` if unset.

See also: [`set_time_limit_sec`](@ref), [`unset_time_limit_sec`](@ref).
"""
function time_limit_sec(model::Model)
    return MOI.get(model, MOI.TimeLimitSec())
end

function _try_get_solver_name(model_like)
    try
        return MOI.get(model_like, MOI.SolverName())::String
    catch ex
        if isa(ex, ArgumentError)
            return "SolverName() attribute not implemented by the optimizer."
        else
            rethrow(ex)
        end
    end
end

"""
    solver_name(model::Model)

If available, returns the `SolverName` property of the underlying optimizer.

Returns `"No optimizer attached"` in `AUTOMATIC` or `MANUAL` modes when no
optimizer is attached.

Returns `"SolverName() attribute not implemented by the optimizer."` if the
attribute is not implemented.
"""
function solver_name(model::Model)
    if mode(model) != DIRECT && MOIU.state(backend(model)) == MOIU.NO_OPTIMIZER
        return "No optimizer attached."
    end
    return _try_get_solver_name(backend(model))
end

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
    return
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
    MOIU.reset_optimizer(backend(model), optimizer)
    return
end

"""
    MOIU.reset_optimizer(model::Model)

Call `MOIU.reset_optimizer` on the backend of `model`.

Cannot be called in direct mode.
"""
function MOIU.reset_optimizer(model::Model)
    error_if_direct_mode(model, :reset_optimizer)
    MOIU.reset_optimizer(backend(model))
    return
end

"""
    MOIU.drop_optimizer(model::Model)

Call `MOIU.drop_optimizer` on the backend of `model`.

Cannot be called in direct mode.
"""
function MOIU.drop_optimizer(model::Model)
    error_if_direct_mode(model, :drop_optimizer)
    MOIU.drop_optimizer(backend(model))
    return
end

"""
    MOIU.attach_optimizer(model::Model)

Call `MOIU.attach_optimizer` on the backend of `model`.

Cannot be called in direct mode.
"""
function MOIU.attach_optimizer(model::Model)
    error_if_direct_mode(model, :attach_optimizer)
    MOIU.attach_optimizer(backend(model))
    return
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
    (@nospecialize optimizer_constructor);
    add_bridges::Bool = true,
)
    error_if_direct_mode(model, :set_optimizer)
    if add_bridges
        optimizer =
            MOI.instantiate(optimizer_constructor; with_bridge_type = Float64)
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
    optimize!(
        model::Model;
        ignore_optimize_hook = (model.optimize_hook === nothing),
        _differentiation_backend::MOI.Nonlinear.AbstractAutomaticDifferentiation =
            MOI.Nonlinear.SparseReverseMode(),
        kwargs...,
    )

Optimize the model.

If an optimizer has not been set yet (see [`set_optimizer`](@ref)), a
[`NoOptimizer`](@ref) error is thrown.

If `ignore_optimize_hook == true`, the optimize hook is ignored and the model is
solved as if the hook was not set. Keyword arguments `kwargs` are passed to the
`optimize_hook`. An error is thrown if `optimize_hook` is `nothing` and keyword
arguments are provided.

## Experimental features

These features may change or be removed in any future version of JuMP.

Pass `_differentiation_backend` to set the
[`MOI.Nonlinear.AbstractAutomaticDifferentiation`](@ref) backend used to compute
derivatives of nonlinear programs.

If you require only `:ExprGraph`, it is more efficient to pass
`_differentiation_backend = MOI.Nonlinear.ExprGraphOnly()`.
"""
function optimize!(
    model::Model;
    ignore_optimize_hook = (model.optimize_hook === nothing),
    _differentiation_backend::MOI.Nonlinear.AbstractAutomaticDifferentiation = MOI.Nonlinear.SparseReverseMode(),
    kwargs...,
)
    # The nlp_model is not kept in sync, so re-set it here.
    # TODO: Consider how to handle incremental solves.
    if nonlinear_model(model) !== nothing
        evaluator = MOI.Nonlinear.Evaluator(
            nonlinear_model(model),
            _differentiation_backend,
            index.(all_variables(model)),
        )
        MOI.set(model, MOI.NLPBlock(), MOI.NLPBlockData(evaluator))
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
    termination_status(model::Model)

Return a [`MOI.TerminationStatusCode`](@ref) describing why the solver stopped
(i.e., the [`MOI.TerminationStatus`](@ref) attribute).
"""
function termination_status(model::Model)
    return MOI.get(model, MOI.TerminationStatus())::MOI.TerminationStatusCode
end

function MOI.get(model::Model, attr::MOI.TerminationStatus)
    if model.is_model_dirty && mode(model) != DIRECT
        return MOI.OPTIMIZE_NOT_CALLED
    end
    return MOI.get(backend(model), attr)
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

"""
    raw_status(model::Model)

Return the reason why the solver stopped in its own words (i.e., the
MathOptInterface model attribute `RawStatusString`).
"""
function raw_status(model::Model)
    if MOI.get(model, MOI.TerminationStatus()) == MOI.OPTIMIZE_NOT_CALLED
        return "optimize not called"
    end
    return MOI.get(model, MOI.RawStatusString())
end

function MOI.get(model::Model, attr::Union{MOI.PrimalStatus,MOI.DualStatus})
    if model.is_model_dirty && mode(model) != DIRECT
        return MOI.NO_SOLUTION
    end
    return MOI.get(backend(model), attr)
end

"""
    primal_status(model::Model; result::Int = 1)

Return a [`MOI.ResultStatusCode`](@ref) describing the status of the most recent
primal solution of the solver (i.e., the [`MOI.PrimalStatus`](@ref) attribute)
associated with the result index `result`.

See also: [`result_count`](@ref).
"""
function primal_status(model::Model; result::Int = 1)
    return MOI.get(model, MOI.PrimalStatus(result))::MOI.ResultStatusCode
end

"""
    dual_status(model::Model; result::Int = 1)

Return a [`MOI.ResultStatusCode`](@ref) describing the status of the most recent
dual solution of the solver (i.e., the [`MOI.DualStatus`](@ref) attribute)
associated with the result index `result`.

See also: [`result_count`](@ref).
"""
function dual_status(model::Model; result::Int = 1)
    return MOI.get(model, MOI.DualStatus(result))::MOI.ResultStatusCode
end

"""
    solve_time(model::Model)

If available, returns the solve time reported by the solver.
Returns "ArgumentError: ModelLike of type `Solver.Optimizer` does not support
accessing the attribute MathOptInterface.SolveTimeSec()" if the attribute is
not implemented.
"""
function solve_time(model::Model)
    return MOI.get(model, MOI.SolveTimeSec())
end

"""
    simplex_iterations(model::Model)

Gets the cumulative number of simplex iterations during the most-recent
optimization.

Solvers must implement `MOI.SimplexIterations()` to use this function.
"""
function simplex_iterations(model::Model)
    return MOI.get(model, MOI.SimplexIterations())
end

"""
    barrier_iterations(model::Model)

Gets the cumulative number of barrier iterations during the most recent
optimization.

Solvers must implement `MOI.BarrierIterations()` to use this function.
"""
function barrier_iterations(model::Model)
    return MOI.get(model, MOI.BarrierIterations())
end

"""
    node_count(model::Model)

Gets the total number of branch-and-bound nodes explored during the most recent
optimization in a Mixed Integer Program.

Solvers must implement `MOI.NodeCount()` to use this function.
"""
function node_count(model::Model)
    return MOI.get(model, MOI.NodeCount())
end

"""
    get(model::Model, attr::MathOptInterface.AbstractOptimizerAttribute)

Return the value of the attribute `attr` from the model's MOI backend.
"""
function MOI.get(model::Model, attr::MOI.AbstractOptimizerAttribute)
    return MOI.get(backend(model), attr)
end

"""
    struct OptimizeNotCalled <: Exception end

A result attribute cannot be queried before [`optimize!`](@ref) is called.
"""
struct OptimizeNotCalled <: Exception end

"""
    struct NoOptimizer <: Exception end

No optimizer is set. The optimizer can be provided to the [`Model`](@ref)
constructor or by calling [`set_optimizer`](@ref).
"""
struct NoOptimizer <: Exception end

# Throws an error if `optimize!` has not been called, i.e., if there is no
# optimizer attached or if the termination status is `MOI.OPTIMIZE_NOT_CALLED`.
function _moi_get_result(model::MOI.ModelLike, args...)
    if MOI.get(model, MOI.TerminationStatus()) == MOI.OPTIMIZE_NOT_CALLED
        throw(OptimizeNotCalled())
    end
    return MOI.get(model, args...)
end

function _moi_get_result(model::MOIU.CachingOptimizer, args...)
    if MOIU.state(model) == MOIU.NO_OPTIMIZER
        throw(NoOptimizer())
    elseif MOI.get(model, MOI.TerminationStatus()) == MOI.OPTIMIZE_NOT_CALLED
        throw(OptimizeNotCalled())
    end
    return MOI.get(model, args...)
end

"""
    get(model::Model, attr::MathOptInterface.AbstractModelAttribute)

Return the value of the attribute `attr` from the model's MOI backend.
"""
function MOI.get(model::Model, attr::MOI.AbstractModelAttribute)
    if !MOI.is_set_by_optimize(attr)
        return MOI.get(backend(model), attr)
    elseif model.is_model_dirty && mode(model) != DIRECT
        @warn(
            "The model has been modified since the last call to `optimize!` (" *
            "or `optimize!` has not been called yet). If you are iteratively " *
            "querying solution information and modifying a model, query all " *
            "the results first, then modify the model.",
        )
        throw(OptimizeNotCalled())
    end
    return _moi_get_result(backend(model), attr)
end

function MOI.get(
    model::Model,
    attr::MOI.AbstractVariableAttribute,
    v::VariableRef,
)
    check_belongs_to_model(v, model)
    if !MOI.is_set_by_optimize(attr)
        return MOI.get(backend(model), attr, index(v))
    elseif model.is_model_dirty && mode(model) != DIRECT
        @warn(
            "The model has been modified since the last call to `optimize!` (" *
            "or `optimize!` has not been called yet). If you are iteratively " *
            "querying solution information and modifying a model, query all " *
            "the results first, then modify the model.",
        )
        throw(OptimizeNotCalled())
    end
    return _moi_get_result(backend(model), attr, index(v))
end

function MOI.get(
    model::Model,
    attr::MOI.AbstractConstraintAttribute,
    cr::ConstraintRef,
)
    check_belongs_to_model(cr, model)
    if !MOI.is_set_by_optimize(attr)
        return MOI.get(backend(model), attr, index(cr))
    elseif model.is_model_dirty && mode(model) != DIRECT
        @warn(
            "The model has been modified since the last call to `optimize!` (" *
            "or `optimize!` has not been called yet). If you are iteratively " *
            "querying solution information and modifying a model, query all " *
            "the results first, then modify the model.",
        )
        throw(OptimizeNotCalled())
    end
    return _moi_get_result(backend(model), attr, index(cr))
end

function MOI.set(m::Model, attr::MOI.AbstractOptimizerAttribute, value)
    m.is_model_dirty = true
    return MOI.set(backend(m), attr, value)
end

function MOI.set(m::Model, attr::MOI.AbstractModelAttribute, value)
    m.is_model_dirty = true
    return MOI.set(backend(m), attr, value)
end

function MOI.set(
    model::Model,
    attr::MOI.AbstractVariableAttribute,
    v::VariableRef,
    value,
)
    check_belongs_to_model(v, model)
    model.is_model_dirty = true
    return MOI.set(backend(model), attr, index(v), value)
end

function MOI.set(
    model::Model,
    attr::MOI.AbstractConstraintAttribute,
    cr::ConstraintRef,
    value,
)
    check_belongs_to_model(cr, model)
    model.is_model_dirty = true
    return MOI.set(backend(model), attr, index(cr), value)
end

_moi_optimizer_index(model::MOI.AbstractOptimizer, index::MOI.Index) = index

function _moi_optimizer_index(model::MOIU.CachingOptimizer, index::MOI.Index)
    if MOIU.state(model) == MOIU.NO_OPTIMIZER
        throw(NoOptimizer())
    elseif MOIU.state(model) == MOIU.EMPTY_OPTIMIZER
        error(
            "There is no `optimizer_index` as the optimizer is not ",
            "synchronized with the cached model. Call ",
            "`MOIU.attach_optimizer(model)` to synchronize it.",
        )
    end
    @assert MOIU.state(model) == MOIU.ATTACHED_OPTIMIZER
    return _moi_optimizer_index(
        model.optimizer,
        model.model_to_optimizer_map[index],
    )
end

function _moi_optimizer_index(
    model::MOI.Bridges.LazyBridgeOptimizer,
    index::MOI.Index,
)
    if index isa MOI.ConstraintIndex && MOI.Bridges.is_bridged(model, index)
        error(
            "There is no `optimizer_index` for $(typeof(index)) constraints",
            " because they are bridged.",
        )
    end
    return _moi_optimizer_index(model.model, index)
end

"""
    optimizer_index(x::VariableRef)::MOI.VariableIndex
    optimizer_index(x::ConstraintRef{Model})::MOI.ConstraintIndex

Return the index that corresponds to `x` in the optimizer model.

Throws [`NoOptimizer`](@ref) if no optimizer is set, and throws an
`ErrorException` if the optimizer is set but is not attached.
"""
function optimizer_index(x::Union{VariableRef,ConstraintRef{Model}})
    model = owner_model(x)
    if mode(model) == DIRECT
        return index(x)
    end
    return _moi_optimizer_index(backend(model), index(x))
end

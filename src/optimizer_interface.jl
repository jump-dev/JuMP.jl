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
attributes using [`set_attribute`](@ref).

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
set_attribute(model, "Presolve", 0)
set_attribute(model, "OutputFlag", 1)
```

## Note

The string names of the attributes are specific to each solver. One should
consult the solver's documentation to find the attributes of interest.

See also: [`set_attribute`](@ref), [`get_attribute`](@ref).
"""
function optimizer_with_attributes(optimizer_constructor, args::Pair...)
    return MOI.OptimizerWithAttributes(optimizer_constructor, args...)
end

"""
    set_optimizer_attribute(
        model::Union{Model,MOI.OptimizerWithAttributes},
        attr::Union{AbstractString,MOI.AbstractOptimizerAttribute},
        value,
    )

Set the solver-specific attribute `attr` in `model` to `value`.

If `attr` is an `AbstractString`, this is equivalent to
`set_optimizer_attribute(model, MOI.RawOptimizerAttribute(name), value)`.

!!! compat
    This method will remain in all v1.X releases of JuMP, but it may be removed
    in a future v2.0 release. We recommend using [`set_attribute`](@ref) instead.

## Example

```julia
set_optimizer_attribute(model, MOI.Silent(), true)
```

See also: [`set_optimizer_attributes`](@ref), [`get_optimizer_attribute`](@ref).
"""
set_optimizer_attribute(model, attr, value) = set_attribute(model, attr, value)

"""
    set_optimizer_attributes(
        model::Union{Model,MOI.OptimizerWithAttributes},
        pairs::Pair...,
    )

Given a list of `attribute => value` pairs, calls
`set_optimizer_attribute(model, attribute, value)` for each pair.

!!! compat
    This method will remain in all v1.X releases of JuMP, but it may be removed
    in a future v2.0 release. We recommend using [`set_attributes`](@ref) instead.

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
        set_attribute(model, name, value)
    end
    return
end

"""
    get_optimizer_attribute(
        model::Union{Model,MOI.OptimizerWithAttributes},
        attr::Union{AbstractString,MOI.AbstractOptimizerAttribute},
    )

Return the value associated with the solver-specific attribute `attr`.

If `attr` is an `AbstractString`, this is equivalent to
`get_optimizer_attribute(model, MOI.RawOptimizerAttribute(name))`.

!!! compat
    This method will remain in all v1.X releases of JuMP, but it may be removed
    in a future v2.0 release. We recommend using [`get_attribute`](@ref) instead.

## Example

```julia
get_optimizer_attribute(model, "SolverSpecificAttributeName")
```

See also: [`set_optimizer_attribute`](@ref), [`set_optimizer_attributes`](@ref).
"""
get_optimizer_attribute(model, attr) = get_attribute(model, attr)

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

See [`set_attribute`](@ref) for setting solver-specific parameters of the
optimizer.

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
            _moi_call_bridge_function(
                MOI.Bridges.add_bridge,
                optimizer,
                bridge_type{Float64},
            )
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
    MOI.set(backend(m), attr, value)
    return
end

function MOI.set(m::Model, attr::MOI.AbstractModelAttribute, value)
    m.is_model_dirty = true
    MOI.set(backend(m), attr, value)
    return
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

"""
    get_attribute(model::Model, attr::MOI.AbstractModelAttribute)
    get_attribute(x::VariableRef, attr::MOI.AbstractVariableAttribute)
    get_attribute(cr::ConstraintRef, attr::MOI.AbstractConstraintAttribute)

Get the value of a solver-specifc attribute `attr`.

This is equivalent to calling [`MOI.get`](@ref) with the associated MOI model
and, for variables and constraints, with the associated [`MOI.VariableIndex`](@ref)
or [`MOI.ConstraintIndex`](@ref).

## Example

```julia
using JuMP
model = Model()
@variable(model, x)
@constraint(model, c, 2 * x <= 1)
get_attribute(model, MOI.Name())
get_attribute(x, MOI.VariableName())
get_attribute(c, MOI.ConstraintName())
```
"""
function get_attribute(model::Model, attr::MOI.AbstractModelAttribute)
    return MOI.get(model, attr)
end

function get_attribute(x::VariableRef, attr::MOI.AbstractVariableAttribute)
    return MOI.get(owner_model(x), attr, x)
end

function get_attribute(cr::ConstraintRef, attr::MOI.AbstractConstraintAttribute)
    return MOI.get(owner_model(cr), attr, cr)
end

"""
    get_attribute(
        model::Union{Model,MOI.OptimizerWithAttributes},
        attr::Union{AbstractString,MOI.AbstractOptimizerAttribute},
    )

Get the value of a solver-specifc attribute `attr`.

This is equivalent to calling [`MOI.get`](@ref) with the associated MOI model.

If `attr` is an `AbstractString`, it is converted to
[`MOI.RawOptimizerAttribute`](@ref).

## Example

```julia
using JuMP, HiGHS
opt = optimizer_with_attributes(HiGHS.Optimizer, "output_flag" => 0)
model = Model(HiGHS.Optimizer)
get_attribute(model, "output_flag")
get_attribute(model, MOI.RawOptimizerAttribute("output_flag"))
get_attribute(opt, "output_flag")
get_attribute(opt, MOI.RawOptimizerAttribute("output_flag"))
```
"""
function get_attribute(
    model::Union{Model,MOI.OptimizerWithAttributes},
    attr::MOI.AbstractOptimizerAttribute,
)
    return MOI.get(model, attr)
end

function get_attribute(
    model::Union{Model,MOI.OptimizerWithAttributes},
    name::String,
)
    return get_attribute(model, MOI.RawOptimizerAttribute(name))
end

# This method is needed for string types like String15 coming from a DataFrame.
function get_attribute(
    model::Union{Model,MOI.OptimizerWithAttributes},
    name::AbstractString,
)
    return get_attribute(model, String(name))
end

"""
    set_attribute(model::Model, attr::MOI.AbstractModelAttribute, value)
    set_attribute(x::VariableRef, attr::MOI.AbstractVariableAttribute, value)
    set_attribute(cr::ConstraintRef, attr::MOI.AbstractConstraintAttribute, value)

Set the value of a solver-specifc attribute `attr` to `value`.

This is equivalent to calling [`MOI.set`](@ref) with the associated MOI model
and, for variables and constraints, with the associated [`MOI.VariableIndex`](@ref)
or [`MOI.ConstraintIndex`](@ref).

## Example

```julia
using JuMP
model = Model()
@variable(model, x)
@constraint(model, c, 2 * x <= 1)
set_attribute(model, MOI.Name(), "model_new")
set_attribute(x, MOI.VariableName(), "x_new")
set_attribute(c, MOI.ConstraintName(), "c_new")
```
"""
function set_attribute(model::Model, attr::MOI.AbstractModelAttribute, value)
    MOI.set(model, attr, value)
    return
end

function set_attribute(
    x::VariableRef,
    attr::MOI.AbstractVariableAttribute,
    value,
)
    MOI.set(owner_model(x), attr, x, value)
    return
end

function set_attribute(
    cr::ConstraintRef,
    attr::MOI.AbstractConstraintAttribute,
    value,
)
    MOI.set(owner_model(cr), attr, cr, value)
    return
end

"""
    set_attribute(
        model::Union{Model,MOI.OptimizerWithAttributes},
        attr::Union{AbstractString,MOI.AbstractOptimizerAttribute},
        value,
    )

Set the value of a solver-specifc attribute `attr` to `value`.

This is equivalent to calling [`MOI.set`](@ref) with the associated MOI model.

If `attr` is an `AbstractString`, it is converted to
[`MOI.RawOptimizerAttribute`](@ref).

## Example

```julia
using JuMP, HiGHS
opt = optimizer_with_attributes(HiGHS.Optimizer, "output_flag" => 0)
model = Model(HiGHS.Optimizer)
set_attribute(model, "output_flag", 1)
set_attribute(model, MOI.RawOptimizerAttribute("output_flag"), 1)
set_attribute(opt, "output_flag", 1)
set_attribute(opt, MOI.RawOptimizerAttribute("output_flag"), 1)
```
"""
function set_attribute(
    model::Union{Model,MOI.OptimizerWithAttributes},
    attr::MOI.AbstractOptimizerAttribute,
    value,
)
    MOI.set(model, attr, value)
    return
end

function set_attribute(
    model::Union{Model,MOI.OptimizerWithAttributes},
    name::String,
    value,
)
    set_attribute(model, MOI.RawOptimizerAttribute(name), value)
    return
end

# This method is needed for string types like String15 coming from a DataFrame.
function set_attribute(
    model::Union{Model,MOI.OptimizerWithAttributes},
    name::AbstractString,
    value,
)
    set_attribute(model, String(name), value)
    return
end

"""
    set_attributes(
        destination::Union{
            Model,
            MOI.OptimizerWithAttributes,
            VariableRef,
            ConstraintRef,
        },
        pairs::Pair...,
    )

Given a list of `attribute => value` pairs, calls
`set_attribute(destination, attribute, value)` for each pair.

## Example

```julia
model = Model(Ipopt.Optimizer)
set_attributes(model, "tol" => 1e-4, "max_iter" => 100)
```
is equivalent to:
```julia
model = Model(Ipopt.Optimizer)
set_attribute(model, "tol", 1e-4)
set_attribute(model, "max_iter", 100)
```

See also: [`set_attribute`](@ref), [`get_attribute`](@ref).
"""
function set_attributes(
    destination::Union{
        Model,
        MOI.OptimizerWithAttributes,
        VariableRef,
        ConstraintRef,
    },
    pairs::Pair...,
)
    for (name, value) in pairs
        set_attribute(destination, name, value)
    end
    return
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

struct ObjectiveFunctionAttribute{A,F}
    attr::A
end

"""
    struct ObjectiveDualStart <: MOI.AbstractModelAttribute end

If the objective function had a dual, it would be `-1` for the Lagrangian
function to be the same.
When the `MOI.Bridges.Objective.SlackBridge` is used, it creates a constraint.
The dual of this constraint is therefore `-1` as well.
When setting this attribute, it allows to set the constraint dual of this
constraint.
"""
struct ObjectiveDualStart <: MOI.AbstractModelAttribute end
# Defining it for `MOI.set` leads to ambiguity
function MOI.throw_set_error_fallback(
    ::MOI.ModelLike,
    ::ObjectiveDualStart,
    value,
)
    return nothing
end

"""
    struct ObjectiveSlackGapPrimalStart <: MOI.AbstractModelAttribute end

If the objective function had a dual, it would be `-1` for the Lagrangian
function to be the same.
When the `MOI.Bridges.Objective.SlackBridge` is used, it creates a constraint.
The dual of this constraint is therefore `-1` as well.
When setting this attribute, it allows to set the constraint dual of this
constraint.
"""
struct ObjectiveSlackGapPrimalStart <: MOI.AbstractModelAttribute end
function MOI.throw_set_error_fallback(
    ::MOI.ModelLike,
    ::ObjectiveSlackGapPrimalStart,
    value,
)
    return nothing
end

function MOI.set(
    b::MOI.Bridges.AbstractBridgeOptimizer,
    attr::ObjectiveFunctionAttribute{A,F},
    value,
) where {A,F}
    obj_attr = MOI.ObjectiveFunction{F}()
    if MOI.Bridges.is_bridged(b, obj_attr)
        return MOI.set(
            MOI.Bridges.recursive_model(b),
            attr,
            MOI.Bridges.bridge(b, obj_attr),
            value,
        )
    else
        return MOI.set(b.model, attr.attr, value)
    end
end

function MOI.set(
    b::MOI.Bridges.AbstractBridgeOptimizer,
    attr::Union{ObjectiveDualStart,ObjectiveSlackGapPrimalStart},
    value,
)
    if MOI.Bridges.is_objective_bridged(b)
        F = MOI.Bridges.Objective.function_type(
            MOI.Bridges.Objective.bridges(b),
        )
        return MOI.set(
            b,
            ObjectiveFunctionAttribute{typeof(attr),F}(attr),
            value,
        )
    else
        return MOI.set(b.model, attr, value)
    end
end

function MOI.set(
    model::MOI.ModelLike,
    ::ObjectiveFunctionAttribute{ObjectiveDualStart},
    b::MOI.Bridges.Objective.SlackBridge,
    value,
)
    return MOI.set(model, MOI.ConstraintDualStart(), b.constraint, value)
end

function MOI.set(
    model::MOI.ModelLike,
    ::ObjectiveFunctionAttribute{ObjectiveSlackGapPrimalStart},
    b::MOI.Bridges.Objective.SlackBridge{T},
    value,
) where {T}
    # `f(x) - slack = value` so `slack = f(x) - value`
    fun = MOI.get(model, MOI.ConstraintFunction(), b.constraint)
    set = MOI.get(model, MOI.ConstraintSet(), b.constraint)
    MOI.Utilities.operate!(-, T, fun, MOI.constant(set))
    # `fun = f - slack` so we remove the term `-slack` to get `f`
    f = MOI.Utilities.remove_variable(fun, b.slack)
    f_val = MOI.Utilities.eval_variables(f) do v
        return MOI.get(model, MOI.VariablePrimalStart(), v)
    end
    MOI.set(model, MOI.VariablePrimalStart(), b.slack, f_val - value)
    return MOI.set(model, MOI.ConstraintPrimalStart(), b.constraint, value)
end

"""
    set_start_values(
        model::Model;
        variable_primal_start::Union{Nothing,Function} = value,
        constraint_primal_start::Union{Nothing,Function} = value,
        constraint_dual_start::Union{Nothing,Function} = dual,
        nonlinear_dual_start::Union{Nothing,Function} = nonlinear_dual_start_value,
    )

Set the primal and dual starting values in `model` using the functions provided.

If any keyword argument is `nothing`, the corresponding start value is skipped.

If the optimizer does not support setting the starting value, the value will be
skipped.

## `variable_primal_start`

This function controls the primal starting solution for the variables. It is
equivalent to calling [`set_start_value`](@ref) for each variable, or setting
the [`MOI.VariablePrimalStart`](@ref) attribute.

If it is a function, it must have the form `variable_primal_start(x::VariableRef)`
that maps each variable `x` to the starting primal value.

The default is [`value`](@ref).

## `constraint_primal_start`

This function controls the primal starting solution for the constraints. It is
equivalent to calling [`set_start_value`](@ref) for each constraint, or setting
the [`MOI.ConstraintPrimalStart`](@ref) attribute.

If it is a function, it must have the form `constraint_primal_start(ci::ConstraintRef)`
that maps each constraint `ci` to the starting primal value.

The default is [`value`](@ref).

## `constraint_dual_start`

This function controls the dual starting solution for the constraints. It is
equivalent to calling [`set_dual_start_value`](@ref) for each constraint, or
setting the [`MOI.ConstraintDualStart`](@ref) attribute.

If it is a function, it must have the form `constraint_dual_start(ci::ConstraintRef)`
that maps each constraint `ci` to the starting dual value.

The default is [`dual`](@ref).

## `nonlinear_dual_start`

This function controls the dual starting solution for the nonlinear constraints
It is equivalent to calling [`set_nonlinear_dual_start_value`](@ref).

If it is a function, it must have the form `nonlinear_dual_start(model::Model)`
that returns a vector corresponding to the dual start of the constraints.

The default is [`nonlinear_dual_start_value`](@ref).
"""
function set_start_values(
    model::Model;
    variable_primal_start::Union{Nothing,Function} = value,
    constraint_primal_start::Union{Nothing,Function} = value,
    constraint_dual_start::Union{Nothing,Function} = dual,
    nonlinear_dual_start::Union{Nothing,Function} = nonlinear_dual_start_value,
)
    variable_primal = Dict{VariableRef,Float64}()
    support_variable_primal = MOI.supports(
        backend(model),
        MOI.VariablePrimalStart(),
        MOI.VariableIndex,
    )
    if support_variable_primal && variable_primal_start !== nothing
        for x in all_variables(model)
            variable_primal[x] = variable_primal_start(x)
        end
    end
    constraint_primal = Dict{ConstraintRef,Float64}()
    constraint_dual = Dict{ConstraintRef,Float64}()
    for (F, S) in list_of_constraint_types(model)
        _get_start_values(
            model,
            F,
            S,
            constraint_primal,
            constraint_primal_start,
            constraint_dual,
            constraint_dual_start,
        )
    end
    if nonlinear_dual_start !== nothing && num_nonlinear_constraints(model) > 0
        if MOI.supports(backend(model), MOI.NLPBlockDualStart())
            nlp_dual_start = nonlinear_dual_start(model)
            set_nonlinear_dual_start_value(model, nlp_dual_start)
        end
    end
    for (x, primal_start) in variable_primal
        set_start_value(x, primal_start)
    end
    for (ci, primal_start) in constraint_primal
        set_start_value(ci, primal_start)
    end
    for (ci, dual_start) in constraint_dual
        set_dual_start_value(ci, dual_start)
    end
    MOI.set(model, ObjectiveDualStart(), -1.0)
    MOI.set(model, ObjectiveSlackGapPrimalStart(), 0.0)
    return
end

function _get_start_values(
    model,
    ::Type{F},
    ::Type{S},
    constraint_primal,
    constraint_primal_start::Union{Nothing,Function},
    constraint_dual,
    constraint_dual_start::Union{Nothing,Function},
) where {F,S}
    moi_model = backend(model)
    CI = MOI.ConstraintIndex{moi_function_type(F),S}
    support_constraint_primal =
        MOI.supports(moi_model, MOI.ConstraintPrimalStart(), CI)
    support_constraint_dual =
        MOI.supports(moi_model, MOI.ConstraintDualStart(), CI)
    for ci in all_constraints(model, F, S)
        if support_constraint_primal && constraint_primal_start !== nothing
            constraint_primal[ci] = constraint_primal_start(ci)
        end
        if support_constraint_dual && constraint_dual_start !== nothing
            constraint_dual[ci] = constraint_dual_start(ci)
        end
    end
    return
end

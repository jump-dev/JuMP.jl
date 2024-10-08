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

See also: [`set_attribute`](@ref), [`get_attribute`](@ref).

## Note

The string names of the attributes are specific to each solver. One should
consult the solver's documentation to find the attributes of interest.

## Example

```jldoctest
julia> import HiGHS

julia> optimizer = optimizer_with_attributes(
           HiGHS.Optimizer, "presolve" => "off", MOI.Silent() => true,
       );

julia> model = Model(optimizer);
```

is equivalent to:

```jldoctest
julia> import HiGHS

julia> model = Model(HiGHS.Optimizer);

julia> set_attribute(model, "presolve", "off")

julia> set_attribute(model, MOI.Silent(), true)
```
"""
function optimizer_with_attributes(optimizer_constructor, args::Pair...)
    return MOI.OptimizerWithAttributes(optimizer_constructor, args...)
end

"""
    set_optimizer_attribute(
        model::Union{GenericModel,MOI.OptimizerWithAttributes},
        attr::Union{AbstractString,MOI.AbstractOptimizerAttribute},
        value,
    )

Set the solver-specific attribute `attr` in `model` to `value`.

If `attr` is an `AbstractString`, this is equivalent to
`set_optimizer_attribute(model, MOI.RawOptimizerAttribute(name), value)`.

!!! compat
    This method will remain in all v1.X releases of JuMP, but it may be removed
    in a future v2.0 release. We recommend using [`set_attribute`](@ref) instead.

See also: [`set_optimizer_attributes`](@ref), [`get_optimizer_attribute`](@ref).

## Example

```jldoctest
julia> model = Model();

julia> set_optimizer_attribute(model, MOI.Silent(), true)
```
"""
set_optimizer_attribute(model, attr, value) = set_attribute(model, attr, value)

"""
    set_optimizer_attributes(
        model::Union{GenericModel,MOI.OptimizerWithAttributes},
        pairs::Pair...,
    )

Given a list of `attribute => value` pairs, calls
`set_optimizer_attribute(model, attribute, value)` for each pair.

!!! compat
    This method will remain in all v1.X releases of JuMP, but it may be removed
    in a future v2.0 release. We recommend using [`set_attributes`](@ref) instead.

See also: [`set_optimizer_attribute`](@ref), [`get_optimizer_attribute`](@ref).

## Example

```jldoctest
julia> import Ipopt

julia> model = Model(Ipopt.Optimizer);

julia> set_optimizer_attributes(model, "tol" => 1e-4, "max_iter" => 100)
```
is equivalent to:
```jldoctest
julia> import Ipopt

julia> model = Model(Ipopt.Optimizer);

julia> set_optimizer_attribute(model, "tol", 1e-4)

julia> set_optimizer_attribute(model, "max_iter", 100)
```
"""
function set_optimizer_attributes(
    model::Union{GenericModel,MOI.OptimizerWithAttributes},
    pairs::Pair...,
)
    for (name, value) in pairs
        set_attribute(model, name, value)
    end
    return
end

"""
    get_optimizer_attribute(
        model::Union{GenericModel,MOI.OptimizerWithAttributes},
        attr::Union{AbstractString,MOI.AbstractOptimizerAttribute},
    )

Return the value associated with the solver-specific attribute `attr`.

If `attr` is an `AbstractString`, this is equivalent to
`get_optimizer_attribute(model, MOI.RawOptimizerAttribute(name))`.

!!! compat
    This method will remain in all v1.X releases of JuMP, but it may be removed
    in a future v2.0 release. We recommend using [`get_attribute`](@ref) instead.

See also: [`set_optimizer_attribute`](@ref), [`set_optimizer_attributes`](@ref).

## Example

```jldoctest
julia> import Ipopt

julia> model = Model(Ipopt.Optimizer);

julia> get_optimizer_attribute(model, MOI.Silent())
false
```
"""
get_optimizer_attribute(model, attr) = get_attribute(model, attr)

"""
    set_silent(model::GenericModel)

Takes precedence over any other attribute controlling verbosity and requires the
solver to produce no output.

See also: [`unset_silent`](@ref).

## Example

```jldoctest
julia> import Ipopt

julia> model = Model(Ipopt.Optimizer);

julia> set_silent(model)

julia> get_attribute(model, MOI.Silent())
true

julia> unset_silent(model)

julia> get_attribute(model, MOI.Silent())
false
```
"""
function set_silent(model::GenericModel)
    return MOI.set(model, MOI.Silent(), true)
end

"""
    unset_silent(model::GenericModel)

Neutralize the effect of the `set_silent` function and let the solver attributes
control the verbosity.

See also: [`set_silent`](@ref).

## Example

```jldoctest
julia> import Ipopt

julia> model = Model(Ipopt.Optimizer);

julia> set_silent(model)

julia> get_attribute(model, MOI.Silent())
true

julia> unset_silent(model)

julia> get_attribute(model, MOI.Silent())
false
```
"""
function unset_silent(model::GenericModel)
    return MOI.set(model, MOI.Silent(), false)
end

"""
    set_time_limit_sec(model::GenericModel, limit::Float64)

Set the time limit (in seconds) of the solver.

Can be unset using [`unset_time_limit_sec`](@ref) or with `limit` set to
`nothing`.

See also: [`unset_time_limit_sec`](@ref), [`time_limit_sec`](@ref).

## Example

```jldoctest
julia> import Ipopt

julia> model = Model(Ipopt.Optimizer);

julia> time_limit_sec(model)

julia> set_time_limit_sec(model, 60.0)

julia> time_limit_sec(model)
60.0

julia> unset_time_limit_sec(model)

julia> time_limit_sec(model)
```
"""
function set_time_limit_sec(model::GenericModel, limit::Real)
    return MOI.set(model, MOI.TimeLimitSec(), convert(Float64, limit))
end

function set_time_limit_sec(model::GenericModel, ::Nothing)
    return unset_time_limit_sec(model)
end

"""
    unset_time_limit_sec(model::GenericModel)

Unset the time limit of the solver.

See also: [`set_time_limit_sec`](@ref), [`time_limit_sec`](@ref).

## Example

```jldoctest
julia> import Ipopt

julia> model = Model(Ipopt.Optimizer);

julia> time_limit_sec(model)

julia> set_time_limit_sec(model, 60.0)

julia> time_limit_sec(model)
60.0

julia> unset_time_limit_sec(model)

julia> time_limit_sec(model)
```
"""
function unset_time_limit_sec(model::GenericModel)
    return MOI.set(model, MOI.TimeLimitSec(), nothing)
end

"""
    time_limit_sec(model::GenericModel)

Return the time limit (in seconds) of the `model`.

Returns `nothing` if unset.

See also: [`set_time_limit_sec`](@ref), [`unset_time_limit_sec`](@ref).

## Example

```jldoctest
julia> import Ipopt

julia> model = Model(Ipopt.Optimizer);

julia> time_limit_sec(model)

julia> set_time_limit_sec(model, 60.0)

julia> time_limit_sec(model)
60.0

julia> unset_time_limit_sec(model)

julia> time_limit_sec(model)
```
"""
function time_limit_sec(model::GenericModel)
    return MOI.get(model, MOI.TimeLimitSec())
end

function _try_get_solver_name(model_like)
    try
        return MOI.get(model_like, MOI.SolverName())::String
    catch ex
        if isa(ex, ArgumentError) || isa(ex, MOI.GetAttributeNotAllowed)
            return "SolverName() attribute not implemented by the optimizer."
        else
            rethrow(ex)
        end
    end
end

"""
    solver_name(model::GenericModel)

If available, returns the [`MOI.SolverName`](@ref) property of the underlying
optimizer.

Returns `"No optimizer attached."` in `AUTOMATIC` or `MANUAL` modes when no
optimizer is attached.

Returns `"SolverName() attribute not implemented by the optimizer."` if the
attribute is not implemented.

## Example

```jldoctest
julia> import Ipopt

julia> model = Model(Ipopt.Optimizer);

julia> solver_name(model)
"Ipopt"

julia> model = Model();

julia> solver_name(model)
"No optimizer attached."

julia> model = Model(MOI.FileFormats.MPS.Model);

julia> solver_name(model)
"SolverName() attribute not implemented by the optimizer."
```
"""
function solver_name(model::GenericModel)
    if mode(model) != DIRECT && MOIU.state(backend(model)) == MOIU.NO_OPTIMIZER
        return "No optimizer attached."
    end
    return _try_get_solver_name(backend(model))
end

"""
    error_if_direct_mode(model::GenericModel, func::Symbol)

Errors if `model` is in direct mode during a call from the function named
`func`.

Used internally within JuMP, or by JuMP extensions who do not want to support
models in direct mode.

## Example

```jldoctest
julia> import HiGHS

julia> model = direct_model(HiGHS.Optimizer());

julia> error_if_direct_mode(model, :foo)
ERROR: The `foo` function is not supported in DIRECT mode.
Stacktrace:
[...]
```
"""
function error_if_direct_mode(model::GenericModel, func::Symbol)
    if mode(model) == DIRECT
        error("The `$func` function is not supported in DIRECT mode.")
    end
    return
end

# These methods directly map to CachingOptimizer methods.

"""
    MOIU.reset_optimizer(model::GenericModel, optimizer::MOI.AbstractOptimizer)

Call `MOIU.reset_optimizer` on the backend of `model`.

Cannot be called in direct mode.
"""
function MOIU.reset_optimizer(
    model::GenericModel,
    optimizer::MOI.AbstractOptimizer,
    ::Bool = true,
)
    error_if_direct_mode(model, :reset_optimizer)
    MOIU.reset_optimizer(backend(model), optimizer)
    return
end

"""
    MOIU.reset_optimizer(model::GenericModel)

Call `MOIU.reset_optimizer` on the backend of `model`.

Cannot be called in direct mode.
"""
function MOIU.reset_optimizer(model::GenericModel)
    error_if_direct_mode(model, :reset_optimizer)
    if MOI.Utilities.state(backend(model)) == MOI.Utilities.ATTACHED_OPTIMIZER
        MOIU.reset_optimizer(backend(model))
    end
    return
end

"""
    MOIU.drop_optimizer(model::GenericModel)

Call `MOIU.drop_optimizer` on the backend of `model`.

Cannot be called in direct mode.
"""
function MOIU.drop_optimizer(model::GenericModel)
    error_if_direct_mode(model, :drop_optimizer)
    MOIU.drop_optimizer(backend(model))
    return
end

"""
    MOIU.attach_optimizer(model::GenericModel)

Call `MOIU.attach_optimizer` on the backend of `model`.

Cannot be called in direct mode.
"""
function MOIU.attach_optimizer(model::GenericModel)
    error_if_direct_mode(model, :attach_optimizer)
    MOIU.attach_optimizer(backend(model))
    return
end

"""
    set_optimizer(
        model::GenericModel,
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

## Example

```jldoctest
julia> import HiGHS

julia> model = Model();

julia> set_optimizer(model, () -> HiGHS.Optimizer())

julia> set_optimizer(model, HiGHS.Optimizer; add_bridges = false)
```
"""
function set_optimizer(
    model::GenericModel{T},
    @nospecialize(optimizer_constructor);
    add_bridges::Bool = true,
) where {T}
    error_if_direct_mode(model, :set_optimizer)
    if add_bridges
        optimizer = MOI.instantiate(optimizer_constructor; with_bridge_type = T)
        for BT in model.bridge_types
            _moi_call_bridge_function(MOI.Bridges.add_bridge, optimizer, BT)
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
        model::GenericModel;
        ignore_optimize_hook = (model.optimize_hook === nothing),
        kwargs...,
    )

Optimize the model.

If an optimizer has not been set yet (see [`set_optimizer`](@ref)), a
[`NoOptimizer`](@ref) error is thrown.

If `ignore_optimize_hook == true`, the optimize hook is ignored and the model is
solved as if the hook was not set. Keyword arguments `kwargs` are passed to the
`optimize_hook`. An error is thrown if `optimize_hook` is `nothing` and keyword
arguments are provided.

## Example

```jldoctest
julia> import HiGHS

julia> model = Model(HiGHS.Optimizer);

julia> set_silent(model)

julia> function my_optimize_hook(model; foo)
           println("Hook called with foo = ", foo)
           return optimize!(model; ignore_optimize_hook = true)
       end
my_optimize_hook (generic function with 1 method)

julia> set_optimize_hook(model, my_optimize_hook)
my_optimize_hook (generic function with 1 method)

julia> optimize!(model; foo = 2)
Hook called with foo = 2
```
"""
function optimize!(
    model::GenericModel;
    ignore_optimize_hook = (model.optimize_hook === nothing),
    # _differentiation_backend is deprecated. Remove in JuMP v2.0
    _differentiation_backend::MOI.Nonlinear.AbstractAutomaticDifferentiation = MOI.Nonlinear.SparseReverseMode(),
    kwargs...,
)
    # The nlp_model is not kept in sync, so re-set it here.
    # TODO: Consider how to handle incremental solves.
    nlp = nonlinear_model(model)
    if nlp !== nothing
        if _uses_new_nonlinear_interface(model)
            error(
                "Cannot optimize a model which contains the features from " *
                "both the legacy (macros beginning with `@NL`) and new " *
                "(`NonlinearExpr`) nonlinear interfaces. You must use one or " *
                "the other.",
            )
        end
        evaluator = MOI.Nonlinear.Evaluator(
            nlp,
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
    optimizer = unsafe_backend(model)
    if !(optimizer isa MOI.AbstractOptimizer)
        error(
            "Cannot call `optimize!` because the provided optimizer is not " *
            "a subtype of `MOI.AbstractOptimizer`.\n\nThe optimizer is:\n\n" *
            sprint(show, optimizer) *
            "\n",
        )
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
                "(that is, NLobjective and NLconstraint).",
            )
        else
            rethrow(err)
        end
    end
    model.is_model_dirty = false
    return
end

function _uses_new_nonlinear_interface(model)
    if objective_function_type(model) <: GenericNonlinearExpr
        return true
    end
    for (F, S) in list_of_constraint_types(model)
        if F <: GenericNonlinearExpr
            return true
        end
    end
    return false
end

"""
    compute_conflict!(model::GenericModel)

Compute a conflict if the model is infeasible.

The conflict is also called the Irreducible Infeasible Subsystem (IIS).

If an optimizer has not been set yet (see [`set_optimizer`](@ref)), a
[`NoOptimizer`](@ref) error is thrown.

The status of the conflict can be checked with the [`MOI.ConflictStatus`](@ref)
model attribute. Then, the status for each constraint can be queried with
the [`MOI.ConstraintConflictStatus`](@ref) attribute.

See also: [`copy_conflict`](@ref)

## Example

```julia
julia> using JuMP

julia> model = Model(Gurobi.Optimizer);

julia> set_silent(model)

julia> @variable(model, x >= 0);

julia> @constraint(model, c1, x >= 2);

julia> @constraint(model, c2, x <= 1);

julia> optimize!(model)

julia> compute_conflict!(model)

julia> get_attribute(model, MOI.ConflictStatus())
CONFLICT_FOUND::ConflictStatusCode = 3
```
"""
function compute_conflict!(model::GenericModel)
    if mode(model) != DIRECT && MOIU.state(backend(model)) == MOIU.NO_OPTIMIZER
        throw(NoOptimizer())
    end
    MOI.compute_conflict!(backend(model))
    return
end

"""
    termination_status(model::GenericModel)

Return a [`MOI.TerminationStatusCode`](@ref) describing why the solver stopped
(that is, the [`MOI.TerminationStatus`](@ref) attribute).

## Example

```jldoctest
julia> import Ipopt

julia> model = Model(Ipopt.Optimizer);

julia> termination_status(model)
OPTIMIZE_NOT_CALLED::TerminationStatusCode = 0
```
"""
function termination_status(model::GenericModel)
    return MOI.get(model, MOI.TerminationStatus())::MOI.TerminationStatusCode
end

function MOI.get(model::GenericModel, attr::MOI.TerminationStatus)
    if model.is_model_dirty && mode(model) != DIRECT
        return MOI.OPTIMIZE_NOT_CALLED
    end
    return MOI.get(backend(model), attr)
end

"""
    result_count(model::GenericModel)

Return the number of results available to query after a call to
[`optimize!`](@ref).

## Example

```jldoctest
julia> import Ipopt

julia> model = Model(Ipopt.Optimizer);

julia> result_count(model)
0
```
"""
function result_count(model::GenericModel)::Int
    if termination_status(model) == MOI.OPTIMIZE_NOT_CALLED
        return 0
    end
    return MOI.get(model, MOI.ResultCount())
end

"""
    raw_status(model::GenericModel)

Return the reason why the solver stopped in its own words (that is, the
MathOptInterface model attribute [`MOI.RawStatusString`](@ref)).

## Example

```jldoctest
julia> import Ipopt

julia> model = Model(Ipopt.Optimizer);

julia> raw_status(model)
"optimize not called"
```
"""
function raw_status(model::GenericModel)
    if MOI.get(model, MOI.TerminationStatus()) == MOI.OPTIMIZE_NOT_CALLED
        return "optimize not called"
    end
    return MOI.get(model, MOI.RawStatusString())
end

function MOI.get(
    model::GenericModel,
    attr::Union{MOI.PrimalStatus,MOI.DualStatus},
)
    if model.is_model_dirty && mode(model) != DIRECT
        return MOI.NO_SOLUTION
    end
    return MOI.get(backend(model), attr)
end

"""
    primal_status(model::GenericModel; result::Int = 1)

Return a [`MOI.ResultStatusCode`](@ref) describing the status of the most recent
primal solution of the solver (that is, the [`MOI.PrimalStatus`](@ref) attribute)
associated with the result index `result`.

See also: [`result_count`](@ref).

## Example

```jldoctest
julia> import Ipopt

julia> model = Model(Ipopt.Optimizer);

julia> primal_status(model; result = 2)
NO_SOLUTION::ResultStatusCode = 0
```
"""
function primal_status(model::GenericModel; result::Int = 1)
    return MOI.get(model, MOI.PrimalStatus(result))::MOI.ResultStatusCode
end

"""
    dual_status(model::GenericModel; result::Int = 1)

Return a [`MOI.ResultStatusCode`](@ref) describing the status of the most recent
dual solution of the solver (that is, the [`MOI.DualStatus`](@ref) attribute)
associated with the result index `result`.

See also: [`result_count`](@ref).

## Example

```jldoctest
julia> import Ipopt

julia> model = Model(Ipopt.Optimizer);

julia> dual_status(model; result = 2)
NO_SOLUTION::ResultStatusCode = 0
```
"""
function dual_status(model::GenericModel; result::Int = 1)
    return MOI.get(model, MOI.DualStatus(result))::MOI.ResultStatusCode
end

"""
    is_solved_and_feasible(
        model::GenericModel;
        allow_local::Bool = true,
        allow_almost::Bool = false,
        dual::Bool = false,
        result::Int = 1,
    )

Return `true` if the model has a feasible primal solution associated with result
index `result` and the [`termination_status`](@ref) is [`OPTIMAL`](@ref) (the
solver found a global optimum) or [`LOCALLY_SOLVED`](@ref) (the solver found a
local optimum, which may also be the global optimum, but the solver could not
prove so).

If `allow_local = false`, then this function returns `true` only if the
[`termination_status`](@ref) is [`OPTIMAL`](@ref).

If `allow_almost = true`, then the [`termination_status`](@ref) may additionally
be [`ALMOST_OPTIMAL`](@ref) or [`ALMOST_LOCALLY_SOLVED`](@ref) (if `allow_local`),
and the [`primal_status`](@ref) and [`dual_status`](@ref) may additionally be
[`NEARLY_FEASIBLE_POINT`](@ref).

If `dual`, additionally check that an optimal dual solution is available.

If this function returns `false`, use [`termination_status`](@ref),
[`result_count`](@ref), [`primal_status`](@ref) and [`dual_status`](@ref) to
understand what solutions are available (if any).

## Example

```jldoctest
julia> import Ipopt

julia> model = Model(Ipopt.Optimizer);

julia> is_solved_and_feasible(model)
false
```
"""
function is_solved_and_feasible(
    model::GenericModel;
    dual::Bool = false,
    allow_local::Bool = true,
    allow_almost::Bool = false,
    result::Int = 1,
)
    status = termination_status(model)
    ret =
        (status == OPTIMAL) ||
        (allow_local && (status == LOCALLY_SOLVED)) ||
        (allow_almost && (status == ALMOST_OPTIMAL)) ||
        (allow_almost && allow_local && (status == ALMOST_LOCALLY_SOLVED))
    if ret
        primal = primal_status(model; result)
        ret &=
            (primal == FEASIBLE_POINT) ||
            (allow_almost && (primal == NEARLY_FEASIBLE_POINT))
    end
    if ret && dual
        dual_stat = dual_status(model; result)
        ret &=
            (dual_stat == FEASIBLE_POINT) ||
            (allow_almost && (dual_stat == NEARLY_FEASIBLE_POINT))
    end
    return ret
end

"""
    solve_time(model::GenericModel)

If available, returns the solve time in wall-clock seconds reported by the
solver (the [`MOI.SolveTimeSec`](@ref) attribute).

Throws a `MOI.GetAttributeNotAllowed` error if the attribute is not implemented
by the solver.

## Example

```jldoctest; filter=r"[0-9].+"
julia> import HiGHS

julia> model = Model(HiGHS.Optimizer);

julia> set_silent(model)

julia> optimize!(model)

julia> solve_time(model)
1.0488089174032211e-5
```
"""
function solve_time(model::GenericModel)
    return MOI.get(model, MOI.SolveTimeSec())
end

"""
    simplex_iterations(model::GenericModel)

If available, returns the cumulative number of simplex iterations during the
most-recent optimization (the [`MOI.SimplexIterations`](@ref) attribute).

Throws a `MOI.GetAttributeNotAllowed` error if the attribute is not implemented
by the solver.

## Example

```jldoctest; filter=r"[0-9].+"
julia> import HiGHS

julia> model = Model(HiGHS.Optimizer);

julia> set_silent(model)

julia> optimize!(model)

julia> simplex_iterations(model)
0
```
"""
function simplex_iterations(model::GenericModel)
    return MOI.get(model, MOI.SimplexIterations())
end

"""
    barrier_iterations(model::GenericModel)

If available, returns the cumulative number of barrier iterations during the
most-recent optimization (the [`MOI.BarrierIterations`](@ref) attribute).

Throws a `MOI.GetAttributeNotAllowed` error if the attribute is not implemented
by the solver.

## Example

```jldoctest; filter=r"[0-9].+"
julia> import HiGHS

julia> model = Model(HiGHS.Optimizer);

julia> set_silent(model)

julia> optimize!(model)

julia> barrier_iterations(model)
0
```
"""
function barrier_iterations(model::GenericModel)
    return MOI.get(model, MOI.BarrierIterations())
end

"""
    node_count(model::GenericModel)

If available, returns the total number of branch-and-bound nodes explored during
the most recent optimization in a Mixed Integer Program (the
[`MOI.NodeCount`](@ref) attribute).

Throws a `MOI.GetAttributeNotAllowed` error if the attribute is not implemented
by the solver.

## Example

```jldoctest; filter=r"[0-9].+"
julia> import HiGHS

julia> model = Model(HiGHS.Optimizer);

julia> set_silent(model)

julia> optimize!(model)

julia> node_count(model)
0
```
"""
function node_count(model::GenericModel)
    return MOI.get(model, MOI.NodeCount())
end

"""
    get(model::GenericModel, attr::MathOptInterface.AbstractOptimizerAttribute)

Return the value of the attribute `attr` from the model's MOI backend.
"""
function MOI.get(model::GenericModel, attr::MOI.AbstractOptimizerAttribute)
    return MOI.get(backend(model), attr)
end

"""
    struct OptimizeNotCalled <: Exception end

An error thrown when a result attribute cannot be queried before
[`optimize!`](@ref) is called.

## Example

```jldoctest
julia> import Ipopt

julia> model = Model(Ipopt.Optimizer);

julia> objective_value(model)
ERROR: OptimizeNotCalled()
Stacktrace:
[...]
```
"""
struct OptimizeNotCalled <: Exception end

"""
    struct NoOptimizer <: Exception end

An error thrown when no optimizer is set and one is required.

The optimizer can be provided to the [`Model`](@ref) constructor or by calling
[`set_optimizer`](@ref).

## Example

```jldoctest
julia> model = Model();

julia> optimize!(model)
ERROR: NoOptimizer()
Stacktrace:
[...]
```
"""
struct NoOptimizer <: Exception end

# Throws an error if `optimize!` has not been called, that is, if there is no
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
    get(model::GenericModel, attr::MathOptInterface.AbstractModelAttribute)

Return the value of the attribute `attr` from the model's MOI backend.
"""
function MOI.get(model::GenericModel, attr::MOI.AbstractModelAttribute)
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
    model::GenericModel,
    attr::MOI.AbstractVariableAttribute,
    v::GenericVariableRef,
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
    model::GenericModel,
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

function MOI.set(m::GenericModel, attr::MOI.AbstractOptimizerAttribute, value)
    m.is_model_dirty = true
    MOI.set(backend(m), attr, value)
    return
end

function MOI.set(m::GenericModel, attr::MOI.AbstractModelAttribute, value)
    m.is_model_dirty = true
    MOI.set(backend(m), attr, value)
    return
end

function MOI.set(
    model::GenericModel,
    attr::MOI.AbstractVariableAttribute,
    v::GenericVariableRef,
    value,
)
    check_belongs_to_model(v, model)
    model.is_model_dirty = true
    return MOI.set(backend(model), attr, index(v), value)
end

function MOI.set(
    model::GenericModel,
    attr::MOI.AbstractConstraintAttribute,
    cr::ConstraintRef,
    value,
)
    check_belongs_to_model(cr, model)
    model.is_model_dirty = true
    return MOI.set(backend(model), attr, index(cr), value)
end

"""
    get_attribute(model::GenericModel, attr::MOI.AbstractModelAttribute)
    get_attribute(x::GenericVariableRef, attr::MOI.AbstractVariableAttribute)
    get_attribute(cr::ConstraintRef, attr::MOI.AbstractConstraintAttribute)

Get the value of a solver-specifc attribute `attr`.

This is equivalent to calling [`MOI.get`](@ref) with the associated MOI model
and, for variables and constraints, with the associated [`MOI.VariableIndex`](@ref)
or [`MOI.ConstraintIndex`](@ref).

## Example

```jldoctest
julia> model = Model();

julia> @variable(model, x)
x

julia> @constraint(model, c, 2 * x <= 1)
c : 2 x ≤ 1

julia> get_attribute(model, MOI.Name())
""

julia> get_attribute(x, MOI.VariableName())
"x"

julia> get_attribute(c, MOI.ConstraintName())
"c"
```
"""
function get_attribute(model::GenericModel, attr::MOI.AbstractModelAttribute)
    return MOI.get(model, attr)
end

function get_attribute(
    x::GenericVariableRef,
    attr::MOI.AbstractVariableAttribute,
)
    return MOI.get(owner_model(x), attr, x)
end

function get_attribute(cr::ConstraintRef, attr::MOI.AbstractConstraintAttribute)
    return MOI.get(owner_model(cr), attr, cr)
end

"""
    get_attribute(
        model::Union{GenericModel,MOI.OptimizerWithAttributes},
        attr::Union{AbstractString,MOI.AbstractOptimizerAttribute},
    )

Get the value of a solver-specifc attribute `attr`.

This is equivalent to calling [`MOI.get`](@ref) with the associated MOI model.

If `attr` is an `AbstractString`, it is converted to
[`MOI.RawOptimizerAttribute`](@ref).

## Example

```jldoctest
julia> import HiGHS

julia> opt = optimizer_with_attributes(HiGHS.Optimizer, "output_flag" => true);

julia> model = Model(opt);

julia> get_attribute(model, "output_flag")
true

julia> get_attribute(model, MOI.RawOptimizerAttribute("output_flag"))
true

julia> get_attribute(opt, "output_flag")
true

julia> get_attribute(opt, MOI.RawOptimizerAttribute("output_flag"))
true
```
"""
function get_attribute(
    model::Union{GenericModel,MOI.OptimizerWithAttributes},
    attr::MOI.AbstractOptimizerAttribute,
)
    return MOI.get(model, attr)
end

function get_attribute(
    model::Union{GenericModel,MOI.OptimizerWithAttributes},
    name::String,
)
    return get_attribute(model, MOI.RawOptimizerAttribute(name))
end

# This method is needed for string types like String15 coming from a DataFrame.
function get_attribute(
    model::Union{GenericModel,MOI.OptimizerWithAttributes},
    name::AbstractString,
)
    return get_attribute(model, String(name))
end

# Some MOI attributes have a strict value type, like ::String or ::Float64, but
# users may pass other generic subtypes, like SubString instead of String.
# Rather than throw obtuse errors, we can catch and fix some common cases. We
# shouldn't fix every case (for example, AbstractString -> String) because
# users might intentionally want the other subtype.
#
# The default case is to do nothing:
_to_concrete_value_type_if_needed(::MOI.AnyAttribute, value) = value

# A common case is passing an AbstractString instead of a String.
function _to_concrete_value_type_if_needed(
    attr::MOI.AnyAttribute,
    value::AbstractString,
)
    if !(value isa String) && MOI.attribute_value_type(attr) === String
        return String(value)
    end
    return value
end

"""
    set_attribute(model::GenericModel, attr::MOI.AbstractModelAttribute, value)
    set_attribute(x::GenericVariableRef, attr::MOI.AbstractVariableAttribute, value)
    set_attribute(cr::ConstraintRef, attr::MOI.AbstractConstraintAttribute, value)

Set the value of a solver-specifc attribute `attr` to `value`.

This is equivalent to calling [`MOI.set`](@ref) with the associated MOI model
and, for variables and constraints, with the associated [`MOI.VariableIndex`](@ref)
or [`MOI.ConstraintIndex`](@ref).

## Example

```jldoctest
julia> model = Model();

julia> @variable(model, x)
x

julia> @constraint(model, c, 2 * x <= 1)
c : 2 x ≤ 1

julia> set_attribute(model, MOI.Name(), "model_new")

julia> set_attribute(x, MOI.VariableName(), "x_new")

julia> set_attribute(c, MOI.ConstraintName(), "c_new")
```
"""
function set_attribute(
    model::GenericModel,
    attr::MOI.AbstractModelAttribute,
    value,
)
    MOI.set(model, attr, _to_concrete_value_type_if_needed(attr, value))
    return
end

function set_attribute(
    x::GenericVariableRef,
    attr::MOI.AbstractVariableAttribute,
    value,
)
    model = owner_model(x)
    MOI.set(model, attr, x, _to_concrete_value_type_if_needed(attr, value))
    return
end

function set_attribute(
    cr::ConstraintRef,
    attr::MOI.AbstractConstraintAttribute,
    value,
)
    model = owner_model(cr)
    MOI.set(model, attr, cr, _to_concrete_value_type_if_needed(attr, value))
    return
end

"""
    set_attribute(
        model::Union{GenericModel,MOI.OptimizerWithAttributes},
        attr::Union{AbstractString,MOI.AbstractOptimizerAttribute},
        value,
    )

Set the value of a solver-specifc attribute `attr` to `value`.

This is equivalent to calling [`MOI.set`](@ref) with the associated MOI model.

If `attr` is an `AbstractString`, it is converted to
[`MOI.RawOptimizerAttribute`](@ref).

## Example

```jldoctest
julia> import HiGHS

julia> opt = optimizer_with_attributes(HiGHS.Optimizer, "output_flag" => false);

julia> model = Model(opt);

julia> set_attribute(model, "output_flag", false)

julia> set_attribute(model, MOI.RawOptimizerAttribute("output_flag"), true)

julia> set_attribute(opt, "output_flag", true)

julia> set_attribute(opt, MOI.RawOptimizerAttribute("output_flag"), false)
```
"""
function set_attribute(
    model::Union{GenericModel,MOI.OptimizerWithAttributes},
    attr::MOI.AbstractOptimizerAttribute,
    value,
)
    MOI.set(model, attr, _to_concrete_value_type_if_needed(attr, value))
    return
end

function set_attribute(
    model::Union{GenericModel,MOI.OptimizerWithAttributes},
    name::String,
    value,
)
    set_attribute(model, MOI.RawOptimizerAttribute(name), value)
    return
end

# This method is needed for string types like String15 coming from a DataFrame.
function set_attribute(
    model::Union{GenericModel,MOI.OptimizerWithAttributes},
    name::AbstractString,
    value,
)
    set_attribute(model, String(name), value)
    return
end

"""
    set_attributes(
        destination::Union{
            GenericModel,
            MOI.OptimizerWithAttributes,
            GenericVariableRef,
            ConstraintRef,
        },
        pairs::Pair...,
    )

Given a list of `attribute => value` pairs, calls
`set_attribute(destination, attribute, value)` for each pair.

See also: [`set_attribute`](@ref), [`get_attribute`](@ref).

## Example

```jldoctest
julia> import Ipopt

julia> model = Model(Ipopt.Optimizer);

julia> set_attributes(model, "tol" => 1e-4, "max_iter" => 100)
```
is equivalent to:
```jldoctest
julia> import Ipopt

julia> model = Model(Ipopt.Optimizer);

julia> set_attribute(model, "tol", 1e-4)

julia> set_attribute(model, "max_iter", 100)
```
"""
function set_attributes(
    destination::Union{
        GenericModel,
        MOI.OptimizerWithAttributes,
        GenericVariableRef,
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
    optimizer_index(x::GenericVariableRef)::MOI.VariableIndex
    optimizer_index(x::ConstraintRef{<:GenericModel})::MOI.ConstraintIndex

Return the variable or constraint index that corresponds to `x` in the
associated model `unsafe_backend(owner_model(x))`.

This function should be used with [`unsafe_backend`](@ref).

As a safer alternative, use [`backend`](@ref) and [`index`](@ref). See the
docstrings of [`backend`](@ref) and [`unsafe_backend`](@ref) for more details.

## Throws

 * Throws [`NoOptimizer`](@ref) if no optimizer is set.
 * Throws an `ErrorException` if the optimizer is set but is not attached.
 * Throws an `ErrorException` if the index is bridged.

## Example

```jldoctest
julia> import HiGHS

julia> model = Model(HiGHS.Optimizer);

julia> set_silent(model)

julia> @variable(model, x >= 0)
x

julia> MOI.Utilities.attach_optimizer(model)

julia> highs = unsafe_backend(model)
A HiGHS model with 1 columns and 0 rows.

julia> optimizer_index(x)
MOI.VariableIndex(1)
```
"""
function optimizer_index(
    x::Union{GenericVariableRef,ConstraintRef{<:GenericModel}},
)
    return _moi_optimizer_index(backend(owner_model(x)), index(x))
end

"""
    set_start_values(
        model::GenericModel;
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

If it is a function, it must have the form `nonlinear_dual_start(model::GenericModel)`
that returns a vector corresponding to the dual start of the constraints.

The default is [`nonlinear_dual_start_value`](@ref).
"""
function set_start_values(
    model::GenericModel{T};
    variable_primal_start::Union{Nothing,Function} = value,
    constraint_primal_start::Union{Nothing,Function} = value,
    constraint_dual_start::Union{Nothing,Function} = dual,
    nonlinear_dual_start::Union{Nothing,Function} = nonlinear_dual_start_value,
) where {T}
    variable_primal = Dict{GenericVariableRef{T},T}()
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
    constraint_primal = Dict{ConstraintRef,T}()
    constraint_dual = Dict{ConstraintRef,T}()
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
    # It's important that we set the variables first, before setting the
    # constraint starts. This is because some bridges, like Constraint.QuadtoSOC,
    # make use of the variable starts when transforming the constraint starts.
    for (x, primal_start) in variable_primal
        set_start_value(x, primal_start)
    end
    for (ci, primal_start) in constraint_primal
        set_start_value(ci, primal_start)
    end
    for (ci, dual_start) in constraint_dual
        set_dual_start_value(ci, dual_start)
    end
    # Needed for models which bridge `min f(x)` into `min t; t >= f(x)`.
    MOI.set(model, MOI.Bridges.Objective.SlackBridgePrimalDualStart(), nothing)
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

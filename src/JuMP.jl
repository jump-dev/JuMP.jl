#  Copyright 2017, Iain Dunning, Joey Huchette, Miles Lubin, and contributors
#  This Source Code Form is subject to the terms of the Mozilla Public
#  License, v. 2.0. If a copy of the MPL was not distributed with this
#  file, You can obtain one at https://mozilla.org/MPL/2.0/.
#############################################################################
# JuMP
# An algebraic modeling language for Julia
# See https://github.com/jump-dev/JuMP.jl
#############################################################################

"""
    JuMP

An algebraic modeling language for Julia.

For more information, go to https://jump.dev.
"""
module JuMP

using LinearAlgebra
using SparseArrays

import MutableArithmetics
const _MA = MutableArithmetics

import MathOptInterface

"""
    MOI

Shorthand for the MathOptInterface package.
"""
const MOI = MathOptInterface

"""
    MOIU

Shorthand for the MathOptInterface.Utilities package.
"""
const MOIU = MOI.Utilities

"""
    MOIB

Shorthand for the MathOptInterface.Bridges package.
"""
const MOIB = MOI.Bridges

import Calculus
import OrderedCollections.OrderedDict
import ForwardDiff

include("Nonlinear/Nonlinear.jl")
include("Containers/Containers.jl")

# Exports are at the end of the file.

include("utils.jl")

const _MOIVAR = MOI.VariableIndex
const _MOICON{F,S} = MOI.ConstraintIndex{F,S}

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

include("shapes.jl")

# Model

"""
    ModelMode

An enum to describe the state of the CachingOptimizer inside a JuMP model.
"""
@enum(ModelMode, AUTOMATIC, MANUAL, DIRECT)
@doc(
    "`moi_backend` field holds a CachingOptimizer in AUTOMATIC mode.",
    AUTOMATIC
)
@doc("`moi_backend` field holds a CachingOptimizer in MANUAL mode.", MANUAL)
@doc(
    "`moi_backend` field holds an AbstractOptimizer. No extra copy of the " *
    "model is stored. The `moi_backend` must support `add_constraint` etc.",
    DIRECT,
)

"""
    AbstractModel

An abstract type that should be subtyped for users creating JuMP extensions.
"""
abstract type AbstractModel end
# All `AbstractModel`s must define methods for these functions:
# num_variables, object_dictionary

"""
    Model

A mathematical model of an optimization problem.
"""
mutable struct Model <: AbstractModel
    # In MANUAL and AUTOMATIC modes, CachingOptimizer.
    # In DIRECT mode, will hold an AbstractOptimizer.
    moi_backend::MOI.AbstractOptimizer
    # List of shapes of constraints that are not `ScalarShape` or `VectorShape`.
    shapes::Dict{_MOICON,AbstractShape}
    # List of bridges to add in addition to the ones added in
    # `MOI.Bridges.full_bridge_optimizer`. With `BridgeableConstraint`, the
    # same bridge may be added many times so we store them in a `Set` instead
    # of, e.g., a `Vector`.
    bridge_types::Set{Any}
    # Hook into a solve call...function of the form f(m::Model; kwargs...),
    # where kwargs get passed along to subsequent solve calls.
    optimize_hook::Any
    # TODO: Document.
    nlp_data::Union{Nothing,Nonlinear.NonlinearData}
    # Dictionary from variable and constraint names to objects.
    obj_dict::Dict{Symbol,Any}
    # Number of times we add large expressions. Incremented and checked by
    # the `operator_warn` method.
    operator_counter::Int
    # A flag to track whether we have modified the model after calling
    # optimize!.
    is_model_dirty::Bool
    # Enable extensions to attach arbitrary information to a JuMP model by
    # using an extension-specific symbol as a key.
    ext::Dict{Symbol,Any}
end

"""
    Model()

Return a new JuMP model without any optimizer; the model is stored in
a cache.

Use [`set_optimizer`](@ref) to set the optimizer before calling
[`optimize!`](@ref).
"""
function Model()
    caching_opt = MOIU.CachingOptimizer(
        MOIU.UniversalFallback(MOIU.Model{Float64}()),
        MOIU.AUTOMATIC,
    )
    return direct_model(caching_opt)
end

"""
    Model(optimizer_factory; add_bridges::Bool = true)

Return a new JuMP model with the provided optimizer and bridge settings.

See [`set_optimizer`](@ref) for the description of the `optimizer_factory` and
`add_bridges` arguments.

## Examples

Create a model with the optimizer set to `Ipopt`:
```julia
model = Model(Ipopt.Optimizer)
```

Pass an anonymous function which creates a `Gurobi.Optimizer`, and disable
bridges:
```julia
env = Gurobi.Env()
model = Model(() -> Gurobi.Optimizer(env); add_bridges = false)
```
"""
function Model(optimizer_factory; add_bridges::Bool = true)
    model = Model()
    set_optimizer(model, optimizer_factory, add_bridges = add_bridges)
    return model
end

"""
    direct_model(backend::MOI.ModelLike)

Return a new JuMP model using [`backend`](@ref) to store the model and solve it.

As opposed to the [`Model`](@ref) constructor, no cache of the model is stored
outside of [`backend`](@ref) and no bridges are automatically applied to
[`backend`](@ref).

## Notes

The absence of a cache reduces the memory footprint but, it is important to bear
in mind the following implications of creating models using this *direct* mode:

* When [`backend`](@ref) does not support an operation, such as modifying
  constraints or adding variables/constraints after solving, an error is
  thrown. For models created using the [`Model`](@ref) constructor, such
  situations can be dealt with by storing the modifications in a cache and
  loading them into the optimizer when `optimize!` is called.
* No constraint bridging is supported by default.
* The optimizer used cannot be changed the model is constructed.
* The model created cannot be copied.
"""
function direct_model(backend::MOI.ModelLike)
    @assert MOI.is_empty(backend)
    return Model(
        backend,
        Dict{_MOICON,AbstractShape}(),
        Set{Any}(),
        nothing,
        nothing,
        Dict{Symbol,Any}(),
        0,
        false,
        Dict{Symbol,Any}(),
    )
end

"""
    direct_model(factory::MOI.OptimizerWithAttributes)

Create a [`direct_model`](@ref) using `factory`, a `MOI.OptimizerWithAttributes`
object created by [`optimizer_with_attributes`](@ref).

## Example

```julia
model = direct_model(
    optimizer_with_attributes(
        Gurobi.Optimizer,
        "Presolve" => 0,
        "OutputFlag" => 1,
    )
)
```
is equivalent to:
```julia
model = direct_model(Gurobi.Optimizer())
set_optimizer_attribute(model, "Presolve", 0)
set_optimizer_attribute(model, "OutputFlag", 1)
```
"""
function direct_model(factory::MOI.OptimizerWithAttributes)
    optimizer = MOI.instantiate(factory)
    return direct_model(optimizer)
end

Base.broadcastable(model::Model) = Ref(model)

"""
    backend(model::Model)

Return the lower-level MathOptInterface model that sits underneath JuMP. This
model depends on which operating mode JuMP is in (see [`mode`](@ref)).

 * If JuMP is in `DIRECT` mode (i.e., the model was created using
   [`direct_model`](@ref)), the backend will be the optimizer passed to
   [`direct_model`](@ref).
 * If JuMP is in `MANUAL` or `AUTOMATIC` mode, the backend is a
   `MOI.Utilities.CachingOptimizer`.

**This function should only be used by advanced users looking to access
low-level MathOptInterface or solver-specific functionality.**

## Notes

If JuMP is not in `DIRECT` mode, the type returned by `backend` may change
between any JuMP releases. Therefore, only use the public API exposed by
MathOptInterface, and do not access internal fields. If you require access to
the innermost optimizer, see [`unsafe_backend`](@ref). Alternatively, use
[`direct_model`](@ref) to create a JuMP model in `DIRECT` mode.

See also: [`unsafe_backend`](@ref).
"""
backend(model::Model) = model.moi_backend

"""
    unsafe_backend(model::Model)

Return the innermost optimizer associated with the JuMP model `model`.

**This function should only be used by advanced users looking to access
low-level solver-specific functionality. It has a high-risk of incorrect usage.
We strongly suggest you use the alternative suggested below.**

See also: [`backend`](@ref).

## Unsafe behavior

This function is unsafe for two main reasons.

First, the formulation and order of variables and constraints in the unsafe
backend may be different to the variables and constraints in `model`. This
can happen because of bridges, or because the solver requires the variables or
constraints in a specific order. In addition, the variable or constraint index
returned by [`index`](@ref) at the JuMP level may be different to the index of
the corresponding variable or constraint in the `unsafe_backend`. There is no
solution to this. Use the alternative suggested below instead.

Second, the `unsafe_backend` may be empty, or lack some modifications made to
the JuMP model. Thus, before calling `unsafe_backend` you should first call
[`MOI.Utilities.attach_optimizer`](@ref) to ensure that the backend is
synchronized with the JuMP model.
```julia
MOI.Utilities.attach_optimizer(model)
inner = unsafe_backend(model)
```

Moreover, if you modify the JuMP model, the reference you have to the backend
(i.e., `inner` in the example above) may be out-dated, and you should call
[`MOI.Utilities.attach_optimizer`](@ref) again.

This function is also unsafe in the reverse direction: if you modify the unsafe
backend, e.g., by adding a new constraint to `inner`, the changes may be
silently discarded by JuMP when the JuMP `model` is modified or solved.

## Alternative

Instead of `unsafe_backend`, create a model using [`direct_model`](@ref) and
call [`backend`](@ref) instead.

For example, instead of:
```julia
model = Model(HiGHS.Optimizer)
@variable(model, x >= 0)
MOI.Utilities.attach_optimizer(model)
highs = unsafe_backend(model)
```
Use:
```julia
model = direct_model(HiGHS.Optimizer())
@variable(model, x >= 0)
highs = backend(model)  # No need to call `attach_optimizer`.
```
"""
unsafe_backend(model::Model) = unsafe_backend(backend(model))

function unsafe_backend(model::MOIU.CachingOptimizer)
    if MOIU.state(model) == MOIU.NO_OPTIMIZER
        error(
            "Unable to get backend optimizer because CachingOptimizer is " *
            "in state `NO_OPTIMIZER`. Call [`set_optimizer`](@ref) first.",
        )
    end
    return unsafe_backend(model.optimizer)
end

unsafe_backend(model::MOIB.LazyBridgeOptimizer) = unsafe_backend(model.model)
unsafe_backend(model::MOI.ModelLike) = model

_moi_mode(::MOI.ModelLike) = DIRECT

function _moi_mode(model::MOIU.CachingOptimizer)
    return model.mode == MOIU.AUTOMATIC ? AUTOMATIC : MANUAL
end

"""
    mode(model::Model)

Return the [`ModelMode`](@ref) ([`DIRECT`](@ref), [`AUTOMATIC`](@ref), or
[`MANUAL`](@ref)) of `model`.
"""
function mode(model::Model)
    # The type of `backend(model)` is not type-stable, so we use a function
    # barrier (`_moi_mode`) to improve performance.
    return _moi_mode(backend(model))
end

# Internal function.
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
    else
        return _try_get_solver_name(backend(model))
    end
end

_moi_bridge_constraints(::MOI.ModelLike) = false

function _moi_bridge_constraints(model::MOIU.CachingOptimizer)
    return model.optimizer isa MOI.Bridges.LazyBridgeOptimizer
end

"""
    bridge_constraints(model::Model)

When in direct mode, return `false`.
When in manual or automatic mode, return a `Bool` indicating whether the
optimizer is set and unsupported constraints are automatically bridged
to equivalent supported constraints when an appropriate transformation is
available.
"""
function bridge_constraints(model::Model)
    # The type of `backend(model)` is not type-stable, so we use a function
    # barrier (`_moi_bridge_constraints`) to improve performance.
    return _moi_bridge_constraints(backend(model))
end

function _moi_add_bridge(
    model::Nothing,
    BridgeType::Type{<:MOI.Bridges.AbstractBridge},
)
    # No optimizer is attached, the bridge will be added when one is attached
    return
end

function _moi_add_bridge(::MOI.ModelLike, ::Type{<:MOI.Bridges.AbstractBridge})
    return error(
        "Cannot add bridge if `add_bridges` was set to `false` in the `Model` ",
        "constructor.",
    )
end

function _moi_add_bridge(
    bridge_opt::MOI.Bridges.LazyBridgeOptimizer,
    BridgeType::Type{<:MOI.Bridges.AbstractBridge},
)
    MOI.Bridges.add_bridge(bridge_opt, BridgeType{Float64})
    return
end

function _moi_add_bridge(
    caching_opt::MOIU.CachingOptimizer,
    BridgeType::Type{<:MOI.Bridges.AbstractBridge},
)
    _moi_add_bridge(caching_opt.optimizer, BridgeType)
    return
end

"""
     add_bridge(model::Model,
                BridgeType::Type{<:MOI.Bridges.AbstractBridge})

Add `BridgeType` to the list of bridges that can be used to transform
unsupported constraints into an equivalent formulation using only constraints
supported by the optimizer.
"""
function add_bridge(
    model::Model,
    BridgeType::Type{<:MOI.Bridges.AbstractBridge},
)
    push!(model.bridge_types, BridgeType)
    # The type of `backend(model)` is not type-stable, so we use a function
    # barrier (`_moi_add_bridge`) to improve performance.
    _moi_add_bridge(JuMP.backend(model), BridgeType)
    return
end

"""
     print_bridge_graph([io::IO,] model::Model)

Print the hyper-graph containing all variable, constraint, and objective types
that could be obtained by bridging the variables, constraints, and objectives
that are present in the model.

Each node in the hyper-graph corresponds to a variable, constraint, or objective
type.
  * Variable nodes are indicated by `[ ]`
  * Constraint nodes are indicated by `( )`
  * Objective nodes are indicated by `| |`
The number inside each pair of brackets is an index of the node in the
hyper-graph.

Note that this hyper-graph is the full list of possible transformations. When
the bridged model is created, we select the shortest hyper-path(s) from this
graph, so many nodes may be un-used.

For more information, see Legat, B., Dowson, O., Garcia, J., and Lubin, M.
(2020).  "MathOptInterface: a data structure for mathematical optimization
problems." URL: [https://arxiv.org/abs/2002.03447](https://arxiv.org/abs/2002.03447)
"""
print_bridge_graph(model::Model) = print_bridge_graph(Base.stdout, model)

function print_bridge_graph(io::IO, model::Model)
    # The type of `backend(model)` is not type-stable, so we use a function
    # barrier (`_moi_print_bridge_graph`) to improve performance.
    return _moi_print_bridge_graph(io, backend(model))
end

function _moi_print_bridge_graph(io::IO, model::MOI.Bridges.LazyBridgeOptimizer)
    return MOI.Bridges.print_graph(io, model)
end

function _moi_print_bridge_graph(io::IO, model::MOIU.CachingOptimizer)
    return _moi_print_bridge_graph(io, model.optimizer)
end

function _moi_print_bridge_graph(::IO, ::MOI.ModelLike)
    return error(
        "Cannot print bridge graph if `add_bridges` was set to `false` in " *
        "the `Model` constructor.",
    )
end

"""
    empty!(model::Model)::Model

Empty the model, that is, remove all variables, constraints and model
attributes but not optimizer attributes. Always return the argument.

Note: removes extensions data.
"""
function Base.empty!(model::Model)::Model
    # The method changes the Model object to, basically, the state it was when
    # created (if the optimizer was already pre-configured). The exceptions
    # are:
    # * optimize_hook: it is basically an optimizer attribute and we promise
    #   to leave them alone (as do MOI.empty!).
    # * bridge_types: for consistency with MOI.empty! for
    #   MOI.Bridges.LazyBridgeOptimizer.
    # * operator_counter: it is just a counter for a single-time warning
    #   message (so keeping it helps to discover inefficiencies).
    MOI.empty!(model.moi_backend)
    empty!(model.shapes)
    model.nlp_data = nothing
    empty!(model.obj_dict)
    empty!(model.ext)
    model.is_model_dirty = false
    return model
end

"""
    isempty(model::Model)

Verifies whether the model is empty, that is, whether the MOI backend
is empty and whether the model is in the same state as at its creation
apart from optimizer attributes.
"""
function Base.isempty(model::Model)
    MOI.is_empty(model.moi_backend) || return false
    isempty(model.shapes) || return false
    model.nlp_data === nothing || return false
    isempty(model.obj_dict) && isempty(model.ext) || return false
    return !model.is_model_dirty
end

"""
    num_variables(model::Model)::Int64

Returns number of variables in `model`.
"""
num_variables(model::Model)::Int64 = MOI.get(model, MOI.NumberOfVariables())

"""
    object_dictionary(model::Model)

Return the dictionary that maps the symbol name of a variable, constraint, or
expression to the corresponding object.

Objects are registered to a specific symbol in the macros.
For example, `@variable(model, x[1:2, 1:2])` registers the array of variables
`x` to the symbol `:x`.

This method should be defined for any subtype of `AbstractModel`.
"""
object_dictionary(model::Model) = model.obj_dict

"""
    unregister(model::Model, key::Symbol)

Unregister the name `key` from `model` so that a new variable, constraint, or
expression can be created with the same key.

Note that this will not delete the object `model[key]`; it will just remove the
reference at `model[key]`. To delete the object, use
```julia
delete(model, model[key])
unregister(model, key)
```

See also: [`object_dictionary`](@ref).

## Examples

```jldoctest; setup=:(model = Model())
julia> @variable(model, x)
x

julia> @variable(model, x)
ERROR: An object of name x is already attached to this model. If
this is intended, consider using the anonymous construction syntax,
e.g., `x = @variable(model, [1:N], ...)` where the name of the object
does not appear inside the macro.

Alternatively, use `unregister(model, :x)` to first unregister the
existing name from the model. Note that this will not delete the object;
it will just remove the reference at `model[:x]`.
[...]

julia> num_variables(model)
1

julia> unregister(model, :x)

julia> @variable(model, x)
x

julia> num_variables(model)
2
```
"""
function unregister(model::AbstractModel, key::Symbol)
    delete!(object_dictionary(model), key)
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
    set_optimize_hook(model::Model, f::Union{Function,Nothing})

Set the function `f` as the optimize hook for `model`.

`f` should have a signature `f(model::Model; kwargs...)`, where the `kwargs` are
those passed to [`optimize!`](@ref).

## Notes

 * The optimize hook should generally modify the model, or some external state
   in some way, and then call `optimize!(model; ignore_optimize_hook = true)` to
   optimize the problem, bypassing the hook.
 * Use `set_optimize_hook(model, nothing)` to unset an optimize hook.

## Examples

```julia
model = Model()
function my_hook(model::Model; kwargs...)
    print(kwargs)
    return optimize!(model; ignore_optimize_hook = true)
end
set_optimize_hook(model, my_hook)
optimize!(model; test_arg = true)
```
"""
set_optimize_hook(model::Model, f) = (model.optimize_hook = f)

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
    set_optimizer_attribute(model::Model, name::String, value)

Sets solver-specific attribute identified by `name` to `value`.

Note that this is equivalent to
`set_optimizer_attribute(model, MOI.RawOptimizerAttribute(name), value)`.

## Example

```julia
set_optimizer_attribute(model, "SolverSpecificAttributeName", true)
```

See also: [`set_optimizer_attributes`](@ref), [`get_optimizer_attribute`](@ref).
"""
function set_optimizer_attribute(model::Model, name::String, value)
    set_optimizer_attribute(model, MOI.RawOptimizerAttribute(name), value)
    return
end

"""
    set_optimizer_attribute(
        model::Model,
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
    model::Model,
    attr::MOI.AbstractOptimizerAttribute,
    value,
)
    MOI.set(model, attr, value)
    return
end

"""
    set_optimizer_attributes(model::Model, pairs::Pair...)

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
function set_optimizer_attributes(model::Model, pairs::Pair...)
    for (name, value) in pairs
        set_optimizer_attribute(model, name, value)
    end
    return
end

"""
    get_optimizer_attribute(model, name::String)

Return the value associated with the solver-specific attribute named `name`.

Note that this is equivalent to
`get_optimizer_attribute(model, MOI.RawOptimizerAttribute(name))`.

## Example

```julia
get_optimizer_attribute(model, "SolverSpecificAttributeName")
```

See also: [`set_optimizer_attribute`](@ref), [`set_optimizer_attributes`](@ref).
"""
function get_optimizer_attribute(model::Model, name::String)
    return get_optimizer_attribute(model, MOI.RawOptimizerAttribute(name))
end

"""
    get_optimizer_attribute(
        model::Model, attr::MOI.AbstractOptimizerAttribute
    )

Return the value of the solver-specific attribute `attr` in `model`.

## Example

```julia
get_optimizer_attribute(model, MOI.Silent())
```

See also: [`set_optimizer_attribute`](@ref), [`set_optimizer_attributes`](@ref).
"""
function get_optimizer_attribute(
    model::Model,
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

"""
    simplex_iterations(model::Model)

Gets the cumulative number of simplex iterations during the most-recent optimization.

Solvers must implement `MOI.SimplexIterations()` to use this function.
"""
function simplex_iterations(model::Model)
    return MOI.get(model, MOI.SimplexIterations())
end

"""
    barrier_iterations(model::Model)

Gets the cumulative number of barrier iterations during the most recent optimization.

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
    AbstractJuMPScalar <: MutableArithmetics.AbstractMutable

Abstract base type for all scalar types

The subtyping of `AbstractMutable` will allow calls of some `Base` functions
to be redirected to a method in MA that handles type promotion more carefully
(e.g. the promotion in sparse matrix products in SparseArrays usually does not
work for JuMP types) and exploits the mutability of `AffExpr` and `QuadExpr`.
"""
abstract type AbstractJuMPScalar <: _MA.AbstractMutable end
Base.ndims(::Type{<:AbstractJuMPScalar}) = 0
Base.ndims(::AbstractJuMPScalar) = 0

# These are required to create symmetric containers of AbstractJuMPScalars.
LinearAlgebra.symmetric_type(::Type{T}) where {T<:AbstractJuMPScalar} = T
LinearAlgebra.symmetric(scalar::AbstractJuMPScalar, ::Symbol) = scalar
# This is required for linear algebra operations involving transposes.
LinearAlgebra.adjoint(scalar::AbstractJuMPScalar) = scalar

"""
    owner_model(s::AbstractJuMPScalar)

Return the model owning the scalar `s`.
"""
function owner_model end

Base.iterate(x::AbstractJuMPScalar) = (x, true)
Base.iterate(::AbstractJuMPScalar, state) = nothing
Base.isempty(::AbstractJuMPScalar) = false

# Check if two arrays of AbstractJuMPScalars are equal. Useful for testing.
function isequal_canonical(
    x::AbstractArray{<:JuMP.AbstractJuMPScalar},
    y::AbstractArray{<:JuMP.AbstractJuMPScalar},
)
    return size(x) == size(y) && all(JuMP.isequal_canonical.(x, y))
end

include("constraints.jl")
include("variables.jl")
include("objective.jl")

function Base.zero(::Type{V}) where {V<:AbstractVariableRef}
    return zero(GenericAffExpr{Float64,V})
end
Base.zero(v::AbstractVariableRef) = zero(typeof(v))
function Base.one(::Type{V}) where {V<:AbstractVariableRef}
    return one(GenericAffExpr{Float64,V})
end
Base.one(v::AbstractVariableRef) = one(typeof(v))

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
    else
        @assert MOIU.state(model) == MOIU.ATTACHED_OPTIMIZER
        return _moi_optimizer_index(
            model.optimizer,
            model.model_to_optimizer_map[index],
        )
    end
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
    else
        return _moi_optimizer_index(model.model, index)
    end
end

"""
    optimizer_index(v::VariableRef)::MOI.VariableIndex

Return the index of the variable that corresponds to `v` in the optimizer model.
It throws [`NoOptimizer`](@ref) if no optimizer is set and throws an
`ErrorException` if the optimizer is set but is not attached.
"""
function optimizer_index(v::VariableRef)
    model = owner_model(v)
    if mode(model) == DIRECT
        return index(v)
    else
        return _moi_optimizer_index(backend(model), index(v))
    end
end

"""
    optimizer_index(cr::ConstraintRef{Model})::MOI.ConstraintIndex

Return the index of the constraint that corresponds to `cr` in the optimizer
model. It throws [`NoOptimizer`](@ref) if no optimizer is set and throws an
`ErrorException` if the optimizer is set but is not attached or if the
constraint is bridged.
"""
function optimizer_index(cr::ConstraintRef{Model})
    if mode(cr.model) == DIRECT
        return index(cr)
    else
        return _moi_optimizer_index(backend(cr.model), index(cr))
    end
end

"""
    index(cr::ConstraintRef)::MOI.ConstraintIndex

Return the index of the constraint that corresponds to `cr` in the MOI backend.
"""
index(cr::ConstraintRef) = cr.index

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

function MOI.get(model::Model, attr::MOI.TerminationStatus)
    if model.is_model_dirty && mode(model) != DIRECT
        return MOI.OPTIMIZE_NOT_CALLED
    end
    return MOI.get(backend(model), attr)
end

function MOI.get(model::Model, attr::Union{MOI.PrimalStatus,MOI.DualStatus})
    if model.is_model_dirty && mode(model) != DIRECT
        return MOI.NO_SOLUTION
    end
    return MOI.get(backend(model), attr)
end

"""
    get(model::Model, attr::MathOptInterface.AbstractOptimizerAttribute)

Return the value of the attribute `attr` from the model's MOI backend.
"""
function MOI.get(model::Model, attr::MOI.AbstractOptimizerAttribute)
    return MOI.get(backend(model), attr)
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

const _Constant = Union{Number,UniformScaling}
_constant_to_number(x::Number) = x
_constant_to_number(J::UniformScaling) = J.λ

# GenericAffineExpression, AffExpr, AffExprConstraint
include("aff_expr.jl")

# GenericQuadExpr, QuadExpr
# GenericQuadConstraint, QuadConstraint
include("quad_expr.jl")

include("mutable_arithmetics.jl")

include("sets.jl")

# Indicator constraint
include("indicator.jl")
# Complementarity constraint
include("complement.jl")
# SDConstraint
include("sd.jl")

"""
    Base.getindex(m::JuMP.AbstractModel, name::Symbol)

To allow easy accessing of JuMP Variables and Constraints via `[]` syntax.
Returns the variable, or group of variables, or constraint, or group of constraints, of the given name which were added to the model. This errors if multiple variables or constraints share the same name.
"""
function Base.getindex(m::JuMP.AbstractModel, name::Symbol)
    obj_dict = object_dictionary(m)
    if !haskey(obj_dict, name)
        throw(KeyError(name))
    elseif obj_dict[name] === nothing
        error(
            "There are multiple variables and/or constraints named $name that are already attached to this model. If creating variables programmatically, use the anonymous variable syntax x = @variable(m, [1:N], ...). If creating constraints programmatically, use the anonymous constraint syntax con = @constraint(m, ...).",
        )
    else
        return obj_dict[name]
    end
end

"""
    Base.setindex!(m::JuMP.AbstractModel, value, name::Symbol)

stores the object `value` in the model `m` using so that it can be accessed via `getindex`.  Can be called with `[]` syntax.
"""
function Base.setindex!(model::AbstractModel, value, name::Symbol)
    # if haskey(object_dictionary(model), name)
    #     warn("Overwriting the object $name stored in the model. Consider using anonymous variables and constraints instead")
    # end
    return object_dictionary(model)[name] = value
end

"""
    haskey(model::AbstractModel, name::Symbol)

Determine whether the model has a mapping for a given name.
"""
function Base.haskey(model::AbstractModel, name::Symbol)
    return haskey(object_dictionary(model), name)
end

"""
    operator_warn(model::AbstractModel)
    operator_warn(model::Model)

This function is called on the model whenever two affine expressions are added
together without using `destructive_add!`, and at least one of the two
expressions has more than 50 terms.

For the case of `Model`, if this function is called more than 20,000 times then
a warning is generated once.
"""
function operator_warn(::AbstractModel) end
function operator_warn(model::Model)
    model.operator_counter += 1
    if model.operator_counter > 20000
        @warn(
            "The addition operator has been used on JuMP expressions a large " *
            "number of times. This warning is safe to ignore but may " *
            "indicate that model generation is slower than necessary. For " *
            "performance reasons, you should not add expressions in a loop. " *
            "Instead of x += y, use add_to_expression!(x,y) to modify x in " *
            "place. If y is a single variable, you may also use " *
            "add_to_expression!(x, coef, y) for x += coef*y.",
            maxlog = 1
        )
    end
end

include("copy.jl")
include("nlp.jl")
include("operators.jl")
include("macros.jl")
include("optimizer_interface.jl")
include("print.jl")
include("solution_summary.jl")
include("lp_sensitivity2.jl")
include("callbacks.jl")
include("file_formats.jl")
include("feasibility_checker.jl")

# MOI contains a number of Enums that are often accessed by users such as
# `MOI.OPTIMAL`. This piece of code re-exports them from JuMP so that users can
# use: `MOI.OPTIMAL`, `JuMP.OPTIMAL`, or `using JuMP; OPTIMAL`.

const ResultStatusCode = MOI.ResultStatusCode
for name in instances(ResultStatusCode)
    @eval const $(Symbol(name)) = $(name)
end

const TerminationStatusCode = MOI.TerminationStatusCode
for name in instances(TerminationStatusCode)
    @eval const $(Symbol(name)) = $(name)
end

const OptimizationSense = MOI.OptimizationSense
for name in instances(OptimizationSense)
    @eval const $(Symbol(name)) = $(name)
end

# JuMP exports everything except internal symbols, which are defined as those
# whose name starts with an underscore. Macros whose names start with
# underscores are internal as well. If you don't want all of these symbols
# in your environment, then use `import JuMP` instead of `using JuMP`.

# Do not add JuMP-defined symbols to this exclude list. Instead, rename them
# with an underscore.
const _EXCLUDE_SYMBOLS = [Symbol(@__MODULE__), :eval, :include]

for sym in names(@__MODULE__, all = true)
    sym_string = string(sym)
    if sym in _EXCLUDE_SYMBOLS ||
       startswith(sym_string, "_") ||
       startswith(sym_string, "@_")
        continue
    end
    if !(
        Base.isidentifier(sym) ||
        (startswith(sym_string, "@") && Base.isidentifier(sym_string[2:end]))
    )
        continue
    end
    @eval export $sym
end

include("precompile.jl")
_precompile_()

end

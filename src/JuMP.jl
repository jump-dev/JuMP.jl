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

import Base.Meta: isexpr, quot
import LinearAlgebra
import MathOptInterface as MOI
import MutableArithmetics
import OrderedCollections
import OrderedCollections: OrderedDict
import Printf
import SparseArrays

# We can't use import MutableArithmetics as _MA because of a bug in MA.
# Fixed in MutableArithmetics v1.2.3, but would require bumping the compat
# bound so we can keep it as this until necessary.
const _MA = MutableArithmetics

"""
    MOIU

Shorthand for the MathOptInterface.Utilities package.
"""
const MOIU = MOI.Utilities

# TODO(odow): remove this constant
const MOIB = MOI.Bridges

# Exports are at the end of the file.

# These imports must come before the definition of `Model`:
include("shapes.jl")

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

"""
    value_type(::Type{<:Union{AbstractModel,AbstractVariableRef}})

Return the return type of [`value`](@ref) for variables of that model. It
defaults to `Float64` if it is not implemented.
"""
value_type(::Type{<:AbstractModel}) = Float64

mutable struct GenericModel{T} <: AbstractModel
    # In MANUAL and AUTOMATIC modes, CachingOptimizer.
    # In DIRECT mode, will hold an AbstractOptimizer.
    moi_backend::MOI.AbstractOptimizer
    # List of shapes of constraints that are not `ScalarShape` or `VectorShape`.
    shapes::Dict{MOI.ConstraintIndex,AbstractShape}
    # List of bridges to add in addition to the ones added in
    # `MOI.Bridges.full_bridge_optimizer`. With `BridgeableConstraint`, the
    # same bridge may be added many times so we store them in a `Set` instead
    # of, e.g., a `Vector`.
    bridge_types::Set{Any}
    # Hook into a solve call...function of the form f(m::GenericModel; kwargs...),
    # where kwargs get passed along to subsequent solve calls.
    optimize_hook::Any
    # TODO: Document.
    nlp_model::Union{Nothing,MOI.Nonlinear.Model}
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
    # A model-level option that is used as the default for the set_string_name
    # keyword to @variable and @constraint.
    set_string_names_on_creation::Bool
end

value_type(::Type{GenericModel{T}}) where {T} = T

function Base.getproperty(model::GenericModel, name::Symbol)
    if name == :nlp_data
        error(
            "The internal field `.nlp_data` was removed from `Model` in JuMP " *
            "v.1.2.0. If you encountered this message without going " *
            "`model.nlp_data`, it means you are using a package that is " *
            "incompatible with your installed version of JuMP. As a " *
            "temporary fix, install a compatible version with " *
            "`import Pkg; Pkg.pkg\"add JuMP@1.1\"`, then restart Julia for " *
            "the changes to take effect. In addition, you should open a " *
            "GitHub issue for the package you are using so that the issue " *
            "can be fixed for future users.",
        )
    end
    return getfield(model, name)
end

"""
    GenericModel{T}()

Return a new JuMP model of value type `T` without any optimizer; the model is
stored in a cache.

Use [`set_optimizer`](@ref) to set the optimizer before calling
[`optimize!`](@ref).
"""
function GenericModel{T}() where {T}
    caching_opt = MOIU.CachingOptimizer(
        MOIU.UniversalFallback(MOIU.Model{T}()),
        MOIU.AUTOMATIC,
    )
    return direct_model(caching_opt, T)
end

const Model = GenericModel{Float64}

"""
    GenericModel{T}(optimizer_factory; add_bridges::Bool = true) where {T}

Return a new JuMP model with the provided optimizer, bridge settings and
[`value_type`](@ref) `T`.

See [`set_optimizer`](@ref) for the description of the `optimizer_factory` and
`add_bridges` arguments.

## Example

```jldoctest
julia> import Ipopt

julia> model = Model(Ipopt.Optimizer);

julia> solver_name(model)
"Ipopt"
```

```jldoctest
julia> import HiGHS

julia> import MultiObjectiveAlgorithms as MOA

julia> model = Model(() -> MOA.Optimizer(HiGHS.Optimizer); add_bridges = false);
```
"""
function GenericModel{T}(
    (@nospecialize optimizer_factory);
    add_bridges::Bool = true,
) where {T}
    model = GenericModel{T}()
    set_optimizer(model, optimizer_factory; add_bridges = add_bridges)
    return model
end

"""
    direct_model(backend::MOI.ModelLike, T::Type = Float64)

Return a new JuMP model using [`backend`](@ref) to store the model and solve it.

As opposed to the [`GenericModel`](@ref) constructor, no cache of the model is stored
outside of [`backend`](@ref) and no bridges are automatically applied to
[`backend`](@ref).

## Notes

The absence of a cache reduces the memory footprint but, it is important to bear
in mind the following implications of creating models using this *direct* mode:

* When [`backend`](@ref) does not support an operation, such as modifying
  constraints or adding variables/constraints after solving, an error is
  thrown. For models created using the [`GenericModel`](@ref) constructor, such
  situations can be dealt with by storing the modifications in a cache and
  loading them into the optimizer when `optimize!` is called.
* No constraint bridging is supported by default.
* The optimizer used cannot be changed the model is constructed.
* The model created cannot be copied.
"""
function direct_model(backend::MOI.ModelLike, T::Type = Float64)
    @assert MOI.is_empty(backend)
    return GenericModel{T}(
        backend,
        Dict{MOI.ConstraintIndex,AbstractShape}(),
        Set{Any}(),
        nothing,
        nothing,
        Dict{Symbol,Any}(),
        0,
        false,
        Dict{Symbol,Any}(),
        true,
    )
end

"""
    direct_model(factory::MOI.OptimizerWithAttributes)

Create a [`direct_model`](@ref) using `factory`, a `MOI.OptimizerWithAttributes`
object created by [`optimizer_with_attributes`](@ref).

## Example

```jldoctest
julia> import HiGHS

julia> optimizer = optimizer_with_attributes(
           HiGHS.Optimizer,
           "presolve" => "off",
           MOI.Silent() => true,
       );

julia> model = direct_model(optimizer)
A JuMP Model
Feasibility problem with:
Variables: 0
Model mode: DIRECT
Solver name: HiGHS
```
is equivalent to:
```jldoctest
julia> model = direct_model(HiGHS.Optimizer())
A JuMP Model
Feasibility problem with:
Variables: 0
Model mode: DIRECT
Solver name: HiGHS

julia> set_attribute(model, "presolve", "off")

julia> set_attribute(model, MOI.Silent(), true)
```
"""
function direct_model(factory::MOI.OptimizerWithAttributes)
    optimizer = MOI.instantiate(factory)
    return direct_model(optimizer)
end

Base.broadcastable(model::GenericModel) = Ref(model)

"""
    backend(model::GenericModel)

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
backend(model::GenericModel) = model.moi_backend

"""
    unsafe_backend(model::GenericModel)

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

```jldoctest
julia> import HiGHS

julia> model = Model(HiGHS.Optimizer)
A JuMP Model
Feasibility problem with:
Variables: 0
Model mode: AUTOMATIC
CachingOptimizer state: EMPTY_OPTIMIZER
Solver name: HiGHS

julia> MOI.Utilities.attach_optimizer(model)

julia> inner = unsafe_backend(model)
A HiGHS model with 0 columns and 0 rows.
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

```jldoctest
julia> import HiGHS

julia> model = Model(HiGHS.Optimizer);

julia> set_silent(model)

julia> @variable(model, x >= 0)
x

julia> MOI.Utilities.attach_optimizer(model)

julia> highs = unsafe_backend(model)
A HiGHS model with 1 columns and 0 rows.
```

Use:

```jldoctest
julia> import HiGHS

julia> model = direct_model(HiGHS.Optimizer());

julia> set_silent(model)

julia> @variable(model, x >= 0)
x

julia> highs = backend(model)  # No need to call `attach_optimizer`.
A HiGHS model with 1 columns and 0 rows.
```
"""
unsafe_backend(model::GenericModel) = unsafe_backend(backend(model))

function unsafe_backend(model::MOIU.CachingOptimizer)
    if MOIU.state(model) == MOIU.NO_OPTIMIZER
        error(
            "Unable to get backend optimizer because CachingOptimizer is " *
            "in state `NO_OPTIMIZER`. Call [`set_optimizer`](@ref) first.",
        )
    end
    return unsafe_backend(model.optimizer)
end

function unsafe_backend(model::MOI.Bridges.LazyBridgeOptimizer)
    return unsafe_backend(model.model)
end

unsafe_backend(model::MOI.ModelLike) = model

_moi_mode(::MOI.ModelLike) = DIRECT

function _moi_mode(model::MOIU.CachingOptimizer)
    return model.mode == MOIU.AUTOMATIC ? AUTOMATIC : MANUAL
end

"""
    mode(model::GenericModel)

Return the [`ModelMode`](@ref) ([`DIRECT`](@ref), [`AUTOMATIC`](@ref), or
[`MANUAL`](@ref)) of `model`.
"""
function mode(model::GenericModel)
    # The type of `backend(model)` is not type-stable, so we use a function
    # barrier (`_moi_mode`) to improve performance.
    return _moi_mode(backend(model))
end

"""
    set_string_names_on_creation(model::GenericModel, value::Bool)

Set the default argument of the `set_string_name` keyword in the
[`@variable`](@ref) and [`@constraint`](@ref) macros to `value`. This is used to
determine whether to assign `String` names to all variables and constraints in
`model`.

By default, `value` is `true`. However, for larger models calling
`set_string_names_on_creation(model, false)` can improve performance at the cost
of reducing the readability of printing and solver log messages.
"""
function set_string_names_on_creation(model::GenericModel, value::Bool)
    model.set_string_names_on_creation = value
    return
end

function set_string_names_on_creation(model::GenericModel)
    return model.set_string_names_on_creation
end

set_string_names_on_creation(::AbstractModel) = true

_moi_bridge_constraints(::MOI.ModelLike) = false

function _moi_bridge_constraints(model::MOIU.CachingOptimizer)
    return model.optimizer isa MOI.Bridges.LazyBridgeOptimizer
end

"""
    bridge_constraints(model::GenericModel)

When in direct mode, return `false`.

When in manual or automatic mode, return a `Bool` indicating whether the
optimizer is set and unsupported constraints are automatically bridged
to equivalent supported constraints when an appropriate transformation is
available.
"""
function bridge_constraints(model::GenericModel)
    # The type of `backend(model)` is not type-stable, so we use a function
    # barrier (`_moi_bridge_constraints`) to improve performance.
    return _moi_bridge_constraints(backend(model))
end

# No optimizer is attached.
_moi_call_bridge_function(::Function, ::Nothing, args...) = nothing

function _moi_call_bridge_function(::Function, ::MOI.ModelLike, args...)
    return error(
        "Cannot use bridge if `add_bridges` was set to `false` in the `Model` ",
        "constructor.",
    )
end

function _moi_call_bridge_function(
    f::Function,
    model::MOI.Bridges.LazyBridgeOptimizer,
    args...,
)
    return f(model, args...)
end

function _moi_call_bridge_function(
    f::Function,
    model::MOI.Utilities.CachingOptimizer,
    args...,
)
    return _moi_call_bridge_function(f, model.optimizer, args...)
end

"""
     add_bridge(model::AbstractModel, BT::Type{<:MOI.Bridges.AbstractBridge})

Add `BT` to the list of bridges that can be used to transform unsupported
constraints into an equivalent formulation using only constraints supported by
the optimizer.

See also: [`remove_bridge`](@ref).
"""
function add_bridge(model::GenericModel, BT::Type{<:MOI.Bridges.AbstractBridge})
    push!(model.bridge_types, BT)
    _moi_call_bridge_function(
        MOI.Bridges.add_bridge,
        backend(model),
        BT{Float64},
    )
    return
end

"""
     remove_bridge(model::GenericModel, BT::Type{<:MOI.Bridges.AbstractBridge})

Remove `BT` to the list of bridges that can be used to transform unsupported
constraints into an equivalent formulation using only constraints supported by
the optimizer.

See also: [`add_bridge`](@ref).
"""
function remove_bridge(
    model::GenericModel{T},
    BT::Type{<:MOI.Bridges.AbstractBridge},
) where {T}
    delete!(model.bridge_types, BT)
    _moi_call_bridge_function(MOI.Bridges.remove_bridge, backend(model), BT{T})
    if mode(model) != DIRECT
        MOI.Utilities.reset_optimizer(model)
    end
    return
end

"""
     print_bridge_graph([io::IO,] model::GenericModel)

Print the hyper-graph containing all variable, constraint, and objective types
that could be obtained by bridging the variables, constraints, and objectives
that are present in the model.

!!! warning
    This function is intended for advanced users. If you want to see only the
    bridges that are currently used, use [`print_active_bridges`](@ref) instead.

## Explanation of output

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
function print_bridge_graph(io::IO, model::GenericModel)
    return _moi_call_bridge_function(backend(model)) do m
        return MOI.Bridges.print_graph(io, m)
    end
end

print_bridge_graph(model::GenericModel) = print_bridge_graph(Base.stdout, model)

"""
    print_active_bridges([io::IO = stdout,] model::GenericModel)

Print a list of the variable, constraint, and objective bridges that are
currently used in the model.
"""
function print_active_bridges(io::IO, model::GenericModel)
    return _moi_call_bridge_function(backend(model)) do m
        return MOI.Bridges.print_active_bridges(io, m)
    end
end

"""
    print_active_bridges([io::IO = stdout,] model::GenericModel, ::Type{F}) where {F}

Print a list of bridges required for an objective function of type `F`.
"""
function print_active_bridges(io::IO, model::GenericModel, ::Type{F}) where {F}
    return _moi_call_bridge_function(backend(model)) do m
        return MOI.Bridges.print_active_bridges(io, m, moi_function_type(F))
    end
end

"""
    print_active_bridges(
        [io::IO = stdout,]
        model::GenericModel,
        F::Type,
        S::Type{<:MOI.AbstractSet},
    )

Print a list of bridges required for a constraint of type `F`-in-`S`.
"""
function print_active_bridges(
    io::IO,
    model::GenericModel,
    F::Type,
    S::Type{<:MOI.AbstractSet},
)
    return _moi_call_bridge_function(backend(model)) do m
        return MOI.Bridges.print_active_bridges(io, m, moi_function_type(F), S)
    end
end

"""
    print_active_bridges(
        [io::IO = stdout,]
        model::GenericModel,
        S::Type{<:MOI.AbstractSet},
    )

Print a list of bridges required to add a variable constrained to the set `S`.
"""
function print_active_bridges(
    io::IO,
    model::GenericModel,
    S::Type{<:MOI.AbstractSet},
)
    return _moi_call_bridge_function(backend(model)) do m
        return MOI.Bridges.print_active_bridges(io, m, S)
    end
end

function print_active_bridges(model::GenericModel, args...)
    return print_active_bridges(Base.stdout, model, args...)
end

"""
    empty!(model::GenericModel)::GenericModel

Empty the model, that is, remove all variables, constraints and model
attributes but not optimizer attributes. Always return the argument.

Note: removes extensions data.
"""
function Base.empty!(model::GenericModel)::GenericModel
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
    model.nlp_model = nothing
    empty!(model.obj_dict)
    empty!(model.ext)
    model.is_model_dirty = false
    return model
end

"""
    isempty(model::GenericModel)

Verifies whether the model is empty, that is, whether the MOI backend
is empty and whether the model is in the same state as at its creation
apart from optimizer attributes.
"""
function Base.isempty(model::GenericModel)
    return MOI.is_empty(model.moi_backend) &&
           isempty(model.shapes) &&
           model.nlp_model === nothing &&
           isempty(model.obj_dict) &&
           isempty(model.ext) &&
           !model.is_model_dirty
end

"""
    object_dictionary(model::GenericModel)

Return the dictionary that maps the symbol name of a variable, constraint, or
expression to the corresponding object.

Objects are registered to a specific symbol in the macros.
For example, `@variable(model, x[1:2, 1:2])` registers the array of variables
`x` to the symbol `:x`.

This method should be defined for any subtype of `AbstractModel`.
"""
object_dictionary(model::GenericModel) = model.obj_dict

"""
    unregister(model::GenericModel, key::Symbol)

Unregister the name `key` from `model` so that a new variable, constraint, or
expression can be created with the same key.

Note that this will not delete the object `model[key]`; it will just remove the
reference at `model[key]`. To delete the object, use [`delete`](@ref) as well.

See also: [`delete`](@ref), [`object_dictionary`](@ref).

## Example

```jldoctest
julia> model = Model();

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
    Base.getindex(m::JuMP.AbstractModel, name::Symbol)

To allow easy accessing of JuMP Variables and Constraints via `[]` syntax.

Returns the variable, or group of variables, or constraint, or group of
constraints, of the given name which were added to the model. This errors if
multiple variables or constraints share the same name.
"""
function Base.getindex(m::AbstractModel, name::Symbol)
    obj_dict = object_dictionary(m)
    if !haskey(obj_dict, name)
        throw(KeyError(name))
    elseif obj_dict[name] === nothing
        error(
            "There are multiple variables and/or constraints named $name " *
            "that are already attached to this model. If creating variables " *
            "programmatically, use the anonymous variable syntax " *
            "`x = @variable(m, [1:N], ...)`. If creating constraints " *
            "programmatically, use the anonymous constraint syntax " *
            "`con = @constraint(m, ...)`.",
        )
    end
    return obj_dict[name]
end

"""
    Base.setindex!(m::JuMP.AbstractModel, value, name::Symbol)

Stores the object `value` in the model `m` using so that it can be accessed via
`getindex`.  Can be called with `[]` syntax.
"""
function Base.setindex!(model::AbstractModel, value, name::Symbol)
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
    set_optimize_hook(model::GenericModel, f::Union{Function,Nothing})

Set the function `f` as the optimize hook for `model`.

`f` should have a signature `f(model::GenericModel; kwargs...)`, where the `kwargs` are
those passed to [`optimize!`](@ref).

## Notes

 * The optimize hook should generally modify the model, or some external state
   in some way, and then call `optimize!(model; ignore_optimize_hook = true)` to
   optimize the problem, bypassing the hook.
 * Use `set_optimize_hook(model, nothing)` to unset an optimize hook.

## Example

```jldoctest
julia> model = Model();

julia> function my_hook(model::Model; kwargs...)
           println(kwargs)
           println("Calling with `ignore_optimize_hook = true`")
           optimize!(model; ignore_optimize_hook = true)
           return
       end
my_hook (generic function with 1 method)

julia> set_optimize_hook(model, my_hook)
my_hook (generic function with 1 method)

julia> optimize!(model; test_arg = true)
Base.Iterators.Pairs{Symbol, Bool, Tuple{Symbol}, NamedTuple{(:test_arg,), Tuple{Bool}}}(:test_arg => 1)
Calling with `ignore_optimize_hook = true`
ERROR: NoOptimizer()
Stacktrace:
[...]
```
"""
set_optimize_hook(model::GenericModel, f) = (model.optimize_hook = f)

"""
    AbstractJuMPScalar <: MutableArithmetics.AbstractMutable

Abstract base type for all scalar types

The subtyping of `AbstractMutable` will allow calls of some `Base` functions
to be redirected to a method in MA that handles type promotion more carefully
(e.g. the promotion in sparse matrix products in SparseArrays usually does not
work for JuMP types) and exploits the mutability of `AffExpr` and `QuadExpr`.
"""
abstract type AbstractJuMPScalar <: _MA.AbstractMutable end

function _get_index_keyword_indexing_error()
    return ErrorException(
        """
        Indexing an `Array` with keyword arguments is not supported.

        Indexing with keyword arguments _is_ supported for `Container.DenseAxisArray`
        (and `Container.SparseAxisArray`).

        Force the container type by passing `container = DenseAxisArray` to any of
        the JuMP macros. For example:

        ```julia
        julia> model = Model();

        julia> @variable(model, x[i=1:3], container = DenseAxisArray)
        1-dimensional DenseAxisArray{VariableRef,1,...} with index sets:
            Dimension 1, Base.OneTo(3)
        And data, a 3-element Vector{VariableRef}:
        x[1]
        x[2]
        x[3]

        julia> x[i=2]
        x[2]
        ```
        """,
    )
end

function Base.getindex(x::Array{<:AbstractJuMPScalar}; kwargs...)
    if isempty(kwargs)
        if length(x) == 1
            return first(x)
        end
        throw(BoundsError(x, tuple()))
    end
    return throw(_get_index_keyword_indexing_error())
end

"""
    owner_model(s::AbstractJuMPScalar)

Return the model owning the scalar `s`.
"""
function owner_model end

Base.ndims(::Type{<:AbstractJuMPScalar}) = 0
Base.ndims(::AbstractJuMPScalar) = 0

# These are required to create symmetric containers of AbstractJuMPScalars.
LinearAlgebra.symmetric_type(::Type{T}) where {T<:AbstractJuMPScalar} = T
LinearAlgebra.hermitian_type(::Type{T}) where {T<:AbstractJuMPScalar} = T
LinearAlgebra.symmetric(scalar::AbstractJuMPScalar, ::Symbol) = scalar
LinearAlgebra.hermitian(scalar::AbstractJuMPScalar, ::Symbol) = adjoint(scalar)
LinearAlgebra.adjoint(scalar::AbstractJuMPScalar) = conj(scalar)

Base.iterate(x::AbstractJuMPScalar) = (x, true)
Base.iterate(::AbstractJuMPScalar, state) = nothing
Base.isempty(::AbstractJuMPScalar) = false

# Check if two arrays of AbstractJuMPScalars are equal. Useful for testing.
function isequal_canonical(
    x::AbstractArray{<:AbstractJuMPScalar},
    y::AbstractArray{<:AbstractJuMPScalar},
)
    return size(x) == size(y) && all(isequal_canonical.(x, y))
end

const _Constant = Union{Number,LinearAlgebra.UniformScaling}
_constant_to_number(x::Number) = x
_constant_to_number(J::LinearAlgebra.UniformScaling) = J.λ

# These includes are inter-dependent, and _must_ come in this particular order.
include("Containers/Containers.jl")
include("constraints.jl")
include("variables.jl")
include("objective.jl")
include("aff_expr.jl")
include("quad_expr.jl")
include("nlp.jl")
include("macros.jl")
include("optimizer_interface.jl")

# These includes are self-contained and can go in any order, except that
# operators.jl must come after mutable_arithmetics.jl
include("callbacks.jl")
include("complement.jl")
include("copy.jl")
include("feasibility_checker.jl")
include("file_formats.jl")
include("lp_sensitivity2.jl")
include("indicator.jl")
include("reified.jl")
include("mutable_arithmetics.jl")
include("operators.jl")
include("sd.jl")
include("sets.jl")
include("solution_summary.jl")

# print.jl must come last, because it uses types defined in earlier files.
include("print.jl")

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

for sym in names(@__MODULE__; all = true)
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

# We need to explicitly export `MOI` because it was exported in previous
# versions, and it is not now because we use `import MathOptInterface as MOI`.
export MOI

# !!! note
#
#     The next codeblock is a SnoopPrecompile workflow that gets run when
#     JuMP.jl is precompiled. Methods that are called here will be cached so
#     that they do not need to be compiled in every session.
#
#     One downside of the design of JuMP (from the perspective of compilation)
#     is that the API has a large surface area, that is, there are many
#     combinations of function-in-set, and they get compiled individually.
#
#     The model below attempts to construct a very diverse model that utilizes
#     most of JuMP's machinery, but if you find a particular model has a bad
#     time-to-first solve (TTFX), then we likely need to add a new type of
#     constraint or objective to the precompile model.

import SnoopPrecompile

SnoopPrecompile.@precompile_all_calls begin
    # Because lots of the work is done by macros, and macros are expanded
    # at lowering time, not much of this would get precompiled without `@eval`
    @eval begin
        let
            model = Model(
                () -> MOI.Utilities.MockOptimizer(
                    MOI.Utilities.UniversalFallback(
                        MOI.Utilities.Model{Float64}(),
                    ),
                ),
            )
            @variables(model, begin
                x1 >= 0
                0 <= x2 <= 1
                x3 == 2
                x4, Bin
                x5, Int
                x6[1:3]
                x7[2:3]
                x8[i = 1:3; isodd(i)], (start = i)
            end)
            @expressions(model, begin
                a, -1 + x1 + x2
                b, 1 + x1^2 + x2
            end)
            @constraints(model, begin
                c1, a >= 0
                c2[i = 2:3, j = ["a"]], a <= i
                c3, a == 0
                c4, 0 <= a <= 1
                c5, b >= 0
                c6, b <= 0
                c7, b == 0
                c8, 0 <= b <= 1
                [x1, x2, x1] in SecondOrderCone()
                [1.0*x1 x2; x2 x1] >= 0, PSDCone()
                1.0 * x1 ⟂ x2
            end)
            @objective(model, Min, x1)
            @objective(model, Max, a)
            @objective(model, Min, b)
            @NLconstraint(model, c9, 1 * sin(x1) + 2.0 * sin(a) + sin(b) <= 1)
            @NLparameter(model, p == 2)
            @NLexpression(model, expr, x1^p)
            @NLobjective(model, Min, 1 + expr)
            optimize!(model)
            # This block could sit in MOI, but it's a public API at the JuMP
            # level, so it can go here.
            #
            # We evaluate with a `view` because it's a common type for solvers
            # like Ipopt to use.
            d = NLPEvaluator(model)
            MOI.features_available(d)
            MOI.initialize(d, [:Grad, :Jac, :Hess])
            g = zeros(num_nonlinear_constraints(model))
            v_g = view(g, 1:length(g))
            x = zeros(num_variables(model))
            MOI.eval_objective(d, x)
            MOI.eval_constraint(d, g, x)
            MOI.eval_constraint(d, v_g, x)
            MOI.eval_objective_gradient(d, g, x)
            J = zeros(length(MOI.jacobian_structure(d)))
            MOI.eval_constraint_jacobian(d, J, x)
            v_J = view(J, 1:length(J))
            MOI.eval_constraint_jacobian(d, v_J, x)
            H = zeros(length(MOI.hessian_lagrangian_structure(d)))
            v_H = view(H, 1:length(H))
            MOI.eval_hessian_lagrangian(d, v_H, x, 1.0, v_g)
        end
    end
end

include("precompile.jl")
_precompile_()

end

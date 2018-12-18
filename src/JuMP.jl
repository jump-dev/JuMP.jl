#  Copyright 2017, Iain Dunning, Joey Huchette, Miles Lubin, and contributors
#  This Source Code Form is subject to the terms of the Mozilla Public
#  License, v. 2.0. If a copy of the MPL was not distributed with this
#  file, You can obtain one at http://mozilla.org/MPL/2.0/.
#############################################################################
# JuMP
# An algebraic modeling language for Julia
# See http://github.com/JuliaOpt/JuMP.jl
#############################################################################

module JuMP

using LinearAlgebra
using SparseArrays

import MathOptInterface
const MOI = MathOptInterface
const MOIU = MOI.Utilities

import Calculus
import DataStructures.OrderedDict
import ForwardDiff
include("Derivatives/Derivatives.jl")
using .Derivatives

export
    Model, VariableRef, AffExpr, QuadExpr,
    with_optimizer,
    NonlinearConstraint,
    ConstraintRef,
    SecondOrderCone, RotatedSecondOrderCone, PSDCone,
    optimize!,
    set_name,
    set_lower_bound, set_upper_bound,
    set_start_value,
    linear_terms,

# Macros and support functions
    @expression, @expressions, @NLexpression, @NLexpressions,
    @variable, @variables, @constraint, @constraints,
    @NLconstraint, @NLconstraints,
    @SDconstraint, @SDconstraints,
    @objective, @NLobjective,
    @NLparameter, @constraintref

# Deprecations for JuMP v0.18 -> JuMP v0.19 transition
Base.@deprecate(getobjectivevalue, JuMP.objective_value)
Base.@deprecate(getobjectivebound, JuMP.objective_bound)
Base.@deprecate(getvalue,          JuMP.value)
Base.@deprecate(getdual,           JuMP.dual)
Base.@deprecate(numvar,            JuMP.num_variables)
Base.@deprecate(numnlconstr,       JuMP.num_nl_constraints)
Base.@deprecate(setlowerbound,     JuMP.set_lower_bound)
Base.@deprecate(setupperbound,     JuMP.set_upper_bound)
Base.@deprecate(linearterms,       JuMP.linear_terms)

include("utils.jl")

const MOIVAR = MOI.VariableIndex
const MOICON{F,S} = MOI.ConstraintIndex{F,S}
const MOILB = MOICON{MOI.SingleVariable,MOI.GreaterThan{Float64}}
const MOIUB = MOICON{MOI.SingleVariable,MOI.LessThan{Float64}}
const MOIFIX = MOICON{MOI.SingleVariable,MOI.EqualTo{Float64}}
const MOIINT = MOICON{MOI.SingleVariable,MOI.Integer}
const MOIBIN = MOICON{MOI.SingleVariable,MOI.ZeroOne}

MOIU.@model(JuMPMOIModel,
            (MOI.ZeroOne, MOI.Integer),
            (MOI.EqualTo, MOI.GreaterThan, MOI.LessThan, MOI.Interval),
            (MOI.Zeros, MOI.Nonnegatives, MOI.Nonpositives, MOI.SecondOrderCone,
             MOI.RotatedSecondOrderCone, MOI.GeometricMeanCone,
             MOI.PositiveSemidefiniteConeTriangle,
             MOI.PositiveSemidefiniteConeSquare,
             MOI.RootDetConeTriangle, MOI.RootDetConeSquare,
             MOI.LogDetConeTriangle, MOI.LogDetConeSquare),
            (),
            (MOI.SingleVariable,),
            (MOI.ScalarAffineFunction, MOI.ScalarQuadraticFunction),
            (MOI.VectorOfVariables,),
            (MOI.VectorAffineFunction,))

"""
    OptimizerFactory

User-friendly closure that creates new MOI models. New `OptimizerFactory`s are
created with [`with_optimizer`](@ref) and new models are created from the
optimizer factory `optimizer_factory` with `optimizer_factory()`.

## Examples

The following construct an optimizer factory and then use it to create two
independent `IpoptOptimizer`s:
```julia
optimizer_factory = with_optimizer(IpoptOptimizer, print_level=0)
optimizer1 = optimizer_factory()
optimizer2 = optimizer_factory()
```
"""
struct OptimizerFactory
    # The constructor can be
    # * `Function`: a function, or
    # * `DataType`: a type, or
    # * `UnionAll`: a type with missing parameters.
    constructor
    args::Tuple
    kwargs # type changes from Julia v0.6 to v0.7 so we leave it untyped for now
end

"""
    with_optimizer(constructor, args...; kwargs...)

Return an `OptimizerFactory` that creates optimizers using the constructor
`constructor` with positional arguments `args` and keyword arguments `kwargs`.

## Examples

The following returns an optimizer factory that creates `IpoptOptimizer`s using
the constructor call `IpoptOptimizer(print_level=0)`:
```julia
with_optimizer(IpoptOptimizer, print_level=0)
```
"""
function with_optimizer(constructor,
                        args...; kwargs...)
    if !applicable(constructor, args...)
        error("$constructor does not have any method with arguments $args.",
              "The first argument of `with_optimizer` should be callable with",
              " the other argument of `with_optimizer`.")
    end
    return OptimizerFactory(constructor, args, kwargs)
end

function (optimizer_factory::OptimizerFactory)()
    return optimizer_factory.constructor(optimizer_factory.args...;
                                         optimizer_factory.kwargs...)
end

###############################################################################
# Model

# Model has three modes:
# 1) AUTOMATIC: moi_backend field holds a CachingOptimizer in AUTOMATIC mode.
# 2) MANUAL: moi_backend field holds a CachingOptimizer in MANUAL mode.
# 3) DIRECT: moi_backend field holds an AbstractOptimizer. No extra copy of the model is stored. The moi_backend must support add_constraint etc.
# Methods to interact with the CachingOptimizer are defined in solverinterface.jl.
@enum ModelMode AUTOMATIC MANUAL DIRECT

abstract type AbstractModel end
# All `AbstractModel`s must define methods for these functions:
# num_variables, object_dictionary

"""
    Model

A mathematical model of an optimization problem.
"""
mutable struct Model <: AbstractModel
    # Special variablewise properties that we keep track of:
    # lower bound, upper bound, fixed, integrality, binary
    variable_to_lower_bound::Dict{MOIVAR, MOILB}
    variable_to_upper_bound::Dict{MOIVAR, MOIUB}
    variable_to_fix::Dict{MOIVAR, MOIFIX}
    variable_to_integrality::Dict{MOIVAR, MOIINT}
    variable_to_zero_one::Dict{MOIVAR, MOIBIN}
    # In MANUAL and AUTOMATIC modes, CachingOptimizer.
    # In DIRECT mode, will hold an AbstractOptimizer.
    moi_backend::MOI.AbstractOptimizer
    # Hook into a solve call...function of the form f(m::Model; kwargs...),
    # where kwargs get passed along to subsequent solve calls.
    optimize_hook
    # TODO: Document.
    nlp_data
    # Dictionary from variable and constraint names to objects.
    obj_dict::Dict{Symbol, Any}
    # Number of times we add large expressions. Incremented and checked by
    # the `operator_warn` method.
    operator_counter::Int
    # Enable extensions to attach arbitrary information to a JuMP model by
    # using an extension-specific symbol as a key.
    ext::Dict{Symbol, Any}
end

"""
    Model(; caching_mode::MOIU.CachingOptimizerMode=MOIU.AUTOMATIC,
            bridge_constraints::Bool=true)

Return a new JuMP model without any optimizer; the model is stored the model in
a cache. The mode of the `CachingOptimizer` storing this cache is
`caching_mode`. The optimizer can be set later in the [`JuMP.optimize!`](@ref)
call. If `bridge_constraints` is true, constraints that are not supported by the
optimizer are automatically bridged to equivalent supported constraints when
an appropriate transformation is defined in the `MathOptInterface.Bridges`
module or is defined in another module and is explicitely added.
"""
function Model(; caching_mode::MOIU.CachingOptimizerMode=MOIU.AUTOMATIC,
                 solver=nothing)
    if solver !== nothing
        error("The solver= keyword is no longer available in JuMP 0.19 and " *
              "later. See the JuMP documentation " *
              "(http://www.juliaopt.org/JuMP.jl/latest/) for latest syntax.")
    end
    universal_fallback = MOIU.UniversalFallback(JuMPMOIModel{Float64}())
    caching_opt = MOIU.CachingOptimizer(universal_fallback,
                                        caching_mode)
    return direct_model(caching_opt)
end

"""
    Model(optimizer_factory::OptimizerFactory;
          caching_mode::MOIU.CachingOptimizerMode=MOIU.AUTOMATIC,
          bridge_constraints::Bool=true)

Return a new JuMP model using the optimizer factory `optimizer_factory` to
create the optimizer. The optimizer factory can be created by the
[`with_optimizer`](@ref) function.

## Examples

The following creates a model using the optimizer
`IpoptOptimizer(print_level=0)`:
```julia
model = JuMP.Model(with_optimizer(IpoptOptimizer, print_level=0))
```
"""
function Model(optimizer_factory::OptimizerFactory;
               bridge_constraints::Bool=true, kwargs...)
    model = Model(; kwargs...)
    set_optimizer(model, optimizer_factory,
                  bridge_constraints=bridge_constraints)
    return model
end

"""
    direct_model(backend::MOI.ModelLike)

Return a new JuMP model using `backend` to store the model and solve it. As
opposed to the [`Model`](@ref) constructor, no cache of the model is stored
outside of `backend` and no bridges are automatically applied to `backend`.
The absence of cache reduces the memory footprint but it is important to bear
in mind the following implications of creating models using this *direct* mode:

* When `backend` does not support an operation such as adding
  variables/constraints after solver or modifying constraints, an error is
  thrown. With models created using the [`Model`](@ref) constructor, such
  situations can be dealt with by storing the modifications in a cache and
  loading them into the optimizer when `JuMP.optimize!` is called.
* No constraint bridging is supported by default.
* The optimizer used cannot be changed the model is constructed.
* The model created cannot be copied.
"""
function direct_model(backend::MOI.ModelLike)
    @assert MOI.is_empty(backend)
    return Model(Dict{MOIVAR, MOILB}(),
                 Dict{MOIVAR, MOIUB}(),
                 Dict{MOIVAR, MOIFIX}(),
                 Dict{MOIVAR, MOIINT}(),
                 Dict{MOIVAR, MOIBIN}(),
                 backend,
                 nothing,
                 nothing,
                 Dict{Symbol, Any}(),
                 0,
                 Dict{Symbol, Any}())
end

Base.broadcastable(model::Model) = Ref(model)


"""
    backend(model::Model)

Return the lower-level MathOptInterface model that sits underneath JuMP. This
model depends on which operating mode JuMP is in (manual, automatic, or direct),
and whether there are any bridges in the model.

If JuMP is in direct mode (i.e., the model was created using [`JuMP.direct_model`](@ref)),
the backend with be the optimizer passed to `direct_model`. If JuMP is in manual
or automatic mode, the backend is a `MOI.Utilities.CachingOptimizer`.

This function should only be used by advanced users looking to access low-level
MathOptInterface or solver-specific functionality.
"""
backend(model::Model) = model.moi_backend

moi_mode(model::MOI.ModelLike) = DIRECT
function moi_mode(model::MOIU.CachingOptimizer)
    if model.mode == MOIU.AUTOMATIC
        return AUTOMATIC
    else
        return MANUAL
    end
end

"""
    mode(model::Model)

Return mode (DIRECT, AUTOMATIC, MANUAL) of model.
"""
mode(model::Model) = moi_mode(backend(model))
# The type of backend(model) is unknown so we directly redirect to another
# function.

# Direct mode
moi_bridge_constraints(model::MOI.ModelLike) = false
function moi_bridge_constraints(model::MOIU.CachingOptimizer)
    return model.optimizer isa MOI.Bridges.LazyBridgeOptimizer
end

# Internal function.
function try_get_solver_name(model_like)
    try
        return MOI.get(model_like, MOI.SolverName())
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
optimizer is attached. Returns
"SolverName() attribute not implemented by the optimizer." if the attribute is
not implemented.
"""
function solver_name(model::Model)
    if mode(model) != DIRECT &&
        MOIU.state(backend(model)) == MOIU.NO_OPTIMIZER
        return "No optimizer attached."
    else
        return try_get_solver_name(backend(model))
    end
end

"""
    bridge_constraints(model::Model)

When in direct mode, return `false`.
When in manual or automatic mode, return a `Bool` indicating whether the
optimizer is set and unsupported constraints are automatically bridged
to equivalent supported constraints when an appropriate transformation is
available.
"""
bridge_constraints(model::Model) = moi_bridge_constraints(backend(model))
# The type of backend(model) is unknown so we directly redirect to another
# function.


"""
    num_variables(model::Model)

Returns number of variables in `model`.
"""
num_variables(model::Model) = MOI.get(model, MOI.NumberOfVariables())

"""
    num_nl_constraints(model::Model)

Returns the number of nonlinear constraints associated with the `model`.
"""
function num_nl_constraints(model::Model)
    return model.nlp_data !== nothing ? length(model.nlp_data.nlconstr) : 0
end

# TODO(IainNZ): Document these too.
object_dictionary(model::Model) = model.obj_dict

"""
    termination_status(model::Model)

Return the reason why the solver stopped (i.e., the MathOptInterface model
attribute `TerminationStatus`).
"""
function termination_status(model::Model)
    return MOI.get(model, MOI.TerminationStatus())
end

"""
    primal_status(model::Model)

Return the status of the most recent primal solution of the solver (i.e., the
MathOptInterface model attribute `PrimalStatus`).
"""
function primal_status(model::Model)
    return MOI.get(model, MOI.PrimalStatus())
end

"""
    dual_status(model::Model)

Return the status of the most recent dual solution of the solver (i.e., the
MathOptInterface model attribute `DualStatus`).
"""
function dual_status(model::Model)
    return MOI.get(model, MOI.DualStatus())
end

set_optimize_hook(model::Model, f) = (model.optimize_hook = f)


# Abstract base type for all scalar types
abstract type AbstractJuMPScalar end


# These are required to create symmetric containers of AbstractJuMPScalars.
LinearAlgebra.symmetric_type(::Type{T}) where T <: AbstractJuMPScalar = T
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
function isequal_canonical(x::AbstractArray{<:JuMP.AbstractJuMPScalar},
                           y::AbstractArray{<:JuMP.AbstractJuMPScalar})
    return size(x) == size(y) && all(JuMP.isequal_canonical.(x, y))
end

include("constraints.jl")
include("variables.jl")
include("objective.jl")

Base.zero(::Type{V}) where V<:AbstractVariableRef = zero(GenericAffExpr{Float64, V})
Base.zero(v::AbstractVariableRef) = zero(typeof(v))
Base.one(::Type{V}) where V<:AbstractVariableRef = one(GenericAffExpr{Float64, V})
Base.one(v::AbstractVariableRef) = one(typeof(v))

mutable struct VariableNotOwnedError <: Exception
    context::String
end
function Base.showerror(io::IO, ex::VariableNotOwnedError)
    print(io, "VariableNotOwnedError: Variable not owned by model present in $(ex.context)")
end

function verify_ownership(m::Model, vec::Vector{VariableRef})
    n = length(vec)
    @inbounds for i in 1:n
        vec[i].m !== m && return false
    end
    return true
end

Base.copy(v::VariableRef, new_model::Model) = VariableRef(new_model, v.index)
Base.copy(x::Nothing, new_model::Model) = nothing
# TODO: Replace with vectorized copy?
Base.copy(v::AbstractArray{VariableRef}, new_model::Model) = (var -> VariableRef(new_model, var.index)).(v)

function optimizer_index(v::VariableRef)
    model = owner_model(v)
    if mode(model) == DIRECT
        return index(v)
    else
        @assert backend(model).state == MOIU.ATTACHED_OPTIMIZER
        return backend(model).model_to_optimizer_map[index(v)]
    end
end

function optimizer_index(cr::ConstraintRef{Model})
    if mode(cr.model) == DIRECT
        return index(cr)
    else
        @assert backend(cr.model).state == MOIU.ATTACHED_OPTIMIZER
        return backend(cr.model).model_to_optimizer_map[index(cr)]
    end
end

index(cr::ConstraintRef) = cr.index

"""
    has_duals(model::Model)

Return true if the solver has a dual solution available to query, otherwise
return false.

See also [`dual`](@ref) and [`shadow_price`](@ref).
"""
has_duals(model::Model) = dual_status(model) != MOI.NO_SOLUTION

"""
    dual(cr::ConstraintRef)

Get the dual value of this constraint in the result returned by a solver.
Use `has_dual` to check if a result exists before asking for values.
See also [`shadow_price`](@ref).
"""
function dual(cr::ConstraintRef{Model, <:MOICON})
    reshape(MOI.get(cr.model, MOI.ConstraintDual(), cr), dual_shape(cr.shape))
end

"""
    get(m::JuMP.Model, attr::MathOptInterface.AbstractModelAttribute)

Return the value of the attribute `attr` from model's MOI backend.
"""
MOI.get(m::Model, attr::MOI.AbstractModelAttribute) = MOI.get(backend(m), attr)
function MOI.get(m::Model, attr::MOI.AbstractVariableAttribute, v::VariableRef)
    @assert m === owner_model(v) # TODO: Improve the error message.
    MOI.get(backend(m), attr, index(v))
end
function MOI.get(m::Model, attr::MOI.AbstractConstraintAttribute, cr::ConstraintRef)
    @assert m === cr.model # TODO: Improve the error message.
    MOI.get(backend(m), attr, index(cr))
end

MOI.set(m::Model, attr::MOI.AbstractModelAttribute, value) = MOI.set(backend(m), attr, value)
function MOI.set(m::Model, attr::MOI.AbstractVariableAttribute, v::VariableRef, value)
    @assert m === owner_model(v) # TODO: Improve the error message.
    MOI.set(backend(m), attr, index(v), value)
end
function MOI.set(m::Model, attr::MOI.AbstractConstraintAttribute, cr::ConstraintRef, value)
    @assert m === cr.model # TODO: Improve the error message.
    MOI.set(backend(m), attr, index(cr), value)
end

###############################################################################
# GenericAffineExpression, AffExpr, AffExprConstraint
include("aff_expr.jl")



###############################################################################
# GenericQuadExpr, QuadExpr
# GenericQuadConstraint, QuadConstraint
include("quad_expr.jl")

##########################################################################
# SOSConstraint  (special ordered set constraints)
# include("sos.jl")

include("sets.jl")

##########################################################################
# SDConstraint
include("sd.jl")

# Internal functions for registering objects in the scope of the model.
function register_object(model::AbstractModel, name::Symbol, value)
    obj_dict = object_dictionary(model)
    if haskey(obj_dict, name)
        error("An object of name $name is already attached to this model. " *
              "If this is intended, consider using the anonymous construction" *
              " syntax, e.g., x = @variable(model, [1:N], ...) where the " *
              "name of the object does not appear inside the macro.")
        obj_dict[name] = nothing
    else
        obj_dict[name] = value
    end
    return value
end

function register_object(model::AbstractModel, name, value)
    error("Invalid name $name.")
end


"""
    Base.getindex(m::JuMP.AbstractModel, name::Symbol)

To allow easy accessing of JuMP tVariables and Constraints via `[]` syntax.
Returns the variable, or group of variables, or constraint, or group of constraints, of the given name which were added to the model. This errors if multiple variables or constraints share the same name.
"""
function Base.getindex(m::JuMP.AbstractModel, name::Symbol)
    obj_dict = object_dictionary(m)
    if !haskey(obj_dict, name)
        throw(KeyError(name))
    elseif obj_dict[name] === nothing
        error("There are multiple variables and/or constraints named $name that are already attached to this model. If creating variables programmatically, use the anonymous variable syntax x = @variable(m, [1:N], ...). If creating constraints programmatically, use the anonymous constraint syntax con = @constraint(m, ...).")
    else
        return obj_dict[name]
    end
end

"""
    Base.setindex!(m::JuMP.Model, value, name::Symbol)

stores the object `value` in the model `m` using so that it can be accessed via `getindex`.  Can be called with `[]` syntax.
"""
function Base.setindex!(m::JuMP.Model, value, name::Symbol)
    # if haskey(m.obj_dict, name)
    #     warn("Overwriting the object $name stored in the model. Consider using anonymous variables and constraints instead")
    # end
    m.obj_dict[name] = value
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
            "add_to_expression!(x, coef, y) for x += coef*y.", maxlog=1)
        # NOTE: On Julia 1.0 (at least), maxlog=1 does not work correctly.
        # See https://github.com/JuliaLang/julia/issues/28786.
    end
end

##########################################################################
# Types used in the nonlinear code
# TODO: rename "m" field to "model" for style compliance
struct NonlinearExpression
    m::Model
    index::Int
end

struct NonlinearParameter <: AbstractJuMPScalar
    m::Model
    index::Int
end

##########################################################################
include("copy.jl")
include("Containers/Containers.jl")
include("operators.jl")
include("macros.jl")
include("optimizer_interface.jl")
include("nlp.jl")
include("print.jl")


##########################################################################
end

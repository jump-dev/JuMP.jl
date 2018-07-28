#  Copyright 2017, Iain Dunning, Joey Huchette, Miles Lubin, and contributors
#  This Source Code Form is subject to the terms of the Mozilla Public
#  License, v. 2.0. If a copy of the MPL was not distributed with this
#  file, You can obtain one at http://mozilla.org/MPL/2.0/.
#############################################################################
# JuMP
# An algebraic modeling language for Julia
# See http://github.com/JuliaOpt/JuMP.jl
#############################################################################

__precompile__()

module JuMP

using Compat
using Compat.LinearAlgebra
using Compat.SparseArrays

import MathOptInterface
const MOI = MathOptInterface
const MOIU = MOI.Utilities

import Calculus
import DataStructures.OrderedDict
using ForwardDiff
include("Derivatives/Derivatives.jl")
using .Derivatives

export
# Objects
    Model, VariableRef, Norm, AffExpr, QuadExpr,
    with_optimizer,
    # LinearConstraint, QuadConstraint, SDConstraint,
    NonlinearConstraint,
    ConstraintRef,
# Cones
    SecondOrderCone, RotatedSecondOrderCone, PSDCone,
# Functions
    # Model related
    setobjectivesense,
    writeLP, writeMPS,
    #addSOS1, addSOS2,
    optimize,
    internalmodel,
    # VariableRef
    setname,
    #getname,
    setlowerbound, setupperbound,
    #getlowerbound, getupperbound,
    #getvalue, setvalue,
    #getdual,
    #setcategory, getcategory,
    setstartvalue,
    linearindex,
    # Expressions and constraints
    linearterms,

# Macros and support functions
    @LinearConstraint, @LinearConstraints, @QuadConstraint, @QuadConstraints,
    @SOCConstraint, @SOCConstraints,
    @expression, @expressions, @NLexpression, @NLexpressions,
    @variable, @variables, @constraint, @constraints,
    @NLconstraint, @NLconstraints,
    @SDconstraint, @SDconstraints,
    @objective, @NLobjective,
    @NLparameter, @constraintref


include("utils.jl")

const MOIVAR = MOI.VariableIndex
const MOICON{F,S} = MOI.ConstraintIndex{F,S}
const MOILB = MOICON{MOI.SingleVariable,MOI.GreaterThan{Float64}}
const MOIUB = MOICON{MOI.SingleVariable,MOI.LessThan{Float64}}
const MOIFIX = MOICON{MOI.SingleVariable,MOI.EqualTo{Float64}}
const MOIINT = MOICON{MOI.SingleVariable,MOI.Integer}
const MOIBIN = MOICON{MOI.SingleVariable,MOI.ZeroOne}

@MOIU.model JuMPMOIModel (ZeroOne, Integer) (EqualTo, GreaterThan, LessThan, Interval) (Zeros, Nonnegatives, Nonpositives, SecondOrderCone, RotatedSecondOrderCone, GeometricMeanCone, PositiveSemidefiniteConeTriangle, PositiveSemidefiniteConeSquare, RootDetConeTriangle, RootDetConeSquare, LogDetConeTriangle, LogDetConeSquare) () (SingleVariable,) (ScalarAffineFunction,ScalarQuadraticFunction) (VectorOfVariables,) (VectorAffineFunction,)

struct Factory
    # The constructor can be
    # * `Function`: a function, or
    # * `DataType`: a type, or
    # * `UnionAll`: a type with missing parameters.
    constructor::Union{Function, DataType, UnionAll}
    args::Tuple
    kwargs # type changes from Julia v0.6 to v0.7 so we leave it untyped for now
end

"""
    with_optimizer(constructor::Type, args...; kwargs...)

Return a factory that creates optimizers using the constructor `constructor`
with positional arguments `args` and keyword arguments `kwargs`.

## Examples

The following returns a factory that creates `IpoptOptimizer`s using the
constructor call `IpoptOptimizer(print_level=0)`:
```julia
with_optimizer(IpoptOptimizer, print_level=0)
```
"""
function with_optimizer(constructor::Type, args...; kwargs...)
    return Factory(constructor, args, kwargs)
end

"""
    create_model(factory::Factory)

Creates a new model with the factory `factory`.
"""
function create_model(factory::Factory)
    return factory.constructor(factory.args...; factory.kwargs...)
end

###############################################################################
# Model

# Model has three modes:
# 1) Automatic: moibackend field holds a LazyBridgeOptimizer{CachingOptimizer} in Automatic mode.
# 2) Manual: moibackend field holds a LazyBridgeOptimizer{CachingOptimizer} in Manual mode.
# 3) Direct: moibackend field holds an AbstractOptimizer. No extra copy of the model is stored. The moibackend must support addconstraint! etc.
# Methods to interact with the CachingOptimizer are defined in solverinterface.jl.
@enum ModelMode Automatic Manual Direct

abstract type AbstractModel end
# All `AbstractModels` must define `num_variables`.

"""
    Model

A mathematical model of an optimization problem.
"""
mutable struct Model <: AbstractModel

    # Special variablewise properties that we keep track of:
    # lower bound, upper bound, fixed, integrality, binary
    variabletolowerbound::Dict{MOIVAR, MOILB}
    variabletoupperbound::Dict{MOIVAR, MOIUB}
    variabletofix::Dict{MOIVAR, MOIFIX}
    variabletointegrality::Dict{MOIVAR, MOIINT}
    variabletozeroone::Dict{MOIVAR, MOIBIN}

    customnames::Vector

    factory::Union{Nothing, Factory}
    # In Manual and Automatic modes, LazyBridgeOptimizer{CachingOptimizer}.
    # In Direct mode, will hold an AbstractOptimizer.
    moibackend::MOI.AbstractOptimizer
    # Hook into a solve call...function of the form f(m::Model; kwargs...),
    # where kwargs get passed along to subsequent solve calls.
    optimizehook
    # TODO: Document.
    nlpdata#::NLPData
    # Dictionary from variable and constraint names to objects.
    objdict::Dict{Symbol, Any}
    # Number of times we add large expressions. Incremented and checked by
    # the `operator_warn` method.
    operator_counter::Int
    # Enable extensions to attach arbitrary information to a JuMP model by
    # using an extension-specific symbol as a key.
    ext::Dict{Symbol, Any}

    # Inner constructor
    function Model(factory::Union{Nothing, Factory}, moibackend::MOI.ModelLike)
        @assert MOI.isempty(moibackend)
        model = new()
        model.variabletolowerbound = Dict{MOIVAR, MOILB}()
        model.variabletoupperbound = Dict{MOIVAR, MOIUB}()
        model.variabletofix = Dict{MOIVAR, MOIFIX}()
        model.variabletointegrality = Dict{MOIVAR, MOIINT}()
        model.variabletozeroone = Dict{MOIVAR, MOIBIN}()
        model.customnames = VariableRef[]
        model.factory = factory
        model.moibackend = moibackend
        model.optimizehook = nothing
        model.nlpdata = nothing
        model.objdict = Dict{Symbol, Any}()
        model.operator_counter = 0
        model.ext = Dict{Symbol, Any}()
        return model
    end
end

"""
    Model(; caching_mode::MOIU.CachingOptimizerMode=MOIU.Automatic,
            bridge_constraints::Bool=true)

Return a new JuMP model without any optimizer storing the model in a cache.
The mode of the `CachingOptimizer` storing this cache is `caching_mode`.
The optimizer can be set later with [`setoptimizer`](@ref). If
`bridge_constraints` is true, constraints that are not supported by the
optimizer are automatically bridged to equivalent supported constraints when
an appropriate is defined in the `MathOptInterface.Bridges` module or is
defined in another module and is explicitely added.
"""
function Model(; caching_mode::MOIU.CachingOptimizerMode=MOIU.Automatic,
                 bridge_constraints::Bool=true)
    universal_fallback = MOIU.UniversalFallback(JuMPMOIModel{Float64}())
    caching_opt = MOIU.CachingOptimizer(universal_fallback,
                                        caching_mode)
    if bridge_constraints
        backend = MOI.Bridges.fullbridgeoptimizer(caching_opt,
                                                  Float64)
    else
        backend = caching_opt
    end
    return Model(nothing, backend)
end

"""
    Model(factory::Factory;
          caching_mode::MOIU.CachingOptimizerMode=MOIU.Automatic,
          bridge_constraints::Bool=true)

Return a new JuMP model using the factory `factory` to create the optimizer.
This is equivalent to calling `Model` with the same keyword arguments and then
calling [`setoptimizer`](@ref) on the created model with the `factory`. The
factory can be created by the [`with_optimizer`](@ref) function.

## Examples

The following creates a model using the optimizer
`IpoptOptimizer(print_level=0)`:
```julia
model = JuMP.Model(with_optimizer(IpoptOptimizer, print_level=0))
```
"""
function Model(factory::Factory; kwargs...)
    model = Model(; kwargs...)
    setoptimizer(model, factory)
    return model
end

"""
    direct_model(backend::MOI.ModelLike)

Return a new JuMP model using `backend` to store the model and solve it. As
opposed to the [`Model`](@ref) constructor, no cache of the model is stored
outside of `backend` and no bridges are automatically applied to `backend`.
The absence of cache reduces the memory footprint but it is importnat to bear
in mind the following implications of creating models using this *direct* mode:

* When `backend` does not support an operation such as adding
  variables/constraints after solver or modifying constraints, an error is
  thrown. With models created using the [`Model`](@ref) constructor, such
  situations can be dealt with modifying the cache only and copying the model
  cache once `JuMP.optimize` is called.
* When `backend` does not support a constraint type, the constraint is not
  automatically bridged to constraints supported by `backend`.
* The optimizer used cannot be changed. With models created using the
  [`Model`](@ref) constuctor, the variable and constraint indices used
  are the indices corresponding to the cached model so the optimizer can be
  changed but in direct mode, changing the backend would render all variable
  and constraint references invalid are their internal indices corresponds to
  the previous backend.
* The model created cannot be copied.
"""
function direct_model(backend::MOI.ModelLike)
    return Model(nothing, backend)
end

# In Automatic and Manual mode, `model.moibackend` is either directly the
# `CachingOptimizer` if `bridge_constraints=false` was passed in the constructor
# or it is a `LazyBridgeOptimizer` and the `CachingOptimizer` is stored in the
# `model` field
function caching_optimizer(model::Model)
    if model.moibackend isa MOIU.CachingOptimizer
        return model.moibackend
    elseif (model.moibackend isa
            MOI.Bridges.LazyBridgeOptimizer{<:MOIU.CachingOptimizer})
        return model.moibackend.model
    else
        error("The function `caching_optimizer` cannot be called on a model " *
              "in `Direct` mode.")
    end
end

"""
    mode(model::Model)

Return mode (Direct, Automatic, Manual) of model.
"""
function mode(model::Model)
    if !(model.moibackend isa MOI.Bridges.LazyBridgeOptimizer{<:MOIU.CachingOptimizer} ||
         model.moibackend isa MOIU.CachingOptimizer)
        return Direct
    elseif caching_optimizer(model).mode == MOIU.Automatic
        return Automatic
    else
        return Manual
    end
end

"""
    num_variables(model::Model)

Returns number of variables in `model`.
"""
num_variables(model::Model) = MOI.get(model, MOI.NumberOfVariables())

"""
    numnlconstr(model::Model)

Returns the number of nonlinear constraints associated with the `model`.
"""
function numnlconstr(model::Model)
    return model.nlpdata !== nothing ? length(model.nlpdata.nlconstr) : 0
end

"""
    objectivebound(model::Model)

Return the best known bound on the optimal objective value after a call to
`optimize(model)`.
"""
objectivebound(model::Model) = MOI.get(model, MOI.ObjectiveBound())

"""
    objectivevalue(model::Model)

Return the objective value after a call to `optimize(model)`.
"""
objectivevalue(model::Model) = MOI.get(model, MOI.ObjectiveValue())

"""
    objectivesense(model::Model)

Return the objective sense, `:Min`, `:Max`, or `:Feasibility`.
"""
function objectivesense(model::Model)
    moisense = MOI.get(model, MOI.ObjectiveSense())
    if moisense == MOI.MinSense
        return :Min
    elseif moisense == MOI.MaxSense
        return :Max
    else
        @assert moisense == MOI.FeasibilitySense
        return :Feasibility
    end
end

# TODO(IainNZ): Document these too.
# TODO(#1381): Implement Base.copy for Model.
terminationstatus(m::Model) = MOI.get(m, MOI.TerminationStatus())
primalstatus(m::Model) = MOI.get(m, MOI.PrimalStatus())
dualstatus(m::Model) = MOI.get(m, MOI.DualStatus())
setoptimizehook(m::Model, f) = (m.optimizehook = f)


#############################################################################
# AbstractConstraint
# Abstract base type for all constraint types
abstract type AbstractConstraint end
# Abstract base type for all scalar types
# In JuMP, used only for VariableRef. Useful primarily for extensions
abstract type AbstractJuMPScalar end


"""
    owner_model(s::AbstractJuMPScalar)

Return the model owning the scalar `s`.
"""
function owner_model end

Base.start(::AbstractJuMPScalar) = false
Base.next(x::AbstractJuMPScalar, state) = (x, true)
Base.done(::AbstractJuMPScalar, state) = state
Base.isempty(::AbstractJuMPScalar) = false

##########################################################################
# Constraint
# Holds the index of a constraint in a Model.
struct ConstraintRef{M<:AbstractModel,C}
    m::M
    index::C
end

# TODO: should model be a parameter here?
function MOI.delete!(m::Model, cr::ConstraintRef{Model})
    @assert m === cr.m
    MOI.delete!(m.moibackend, index(cr))
end

MOI.isvalid(m::Model, cr::ConstraintRef{Model}) = cr.m === m && MOI.isvalid(m.moibackend, cr.index)

"""
    addconstraint(m::Model, c::AbstractConstraint, name::String="")

Add a constraint `c` to `Model m` and sets its name.
"""
function addconstraint(m::Model, c::AbstractConstraint, name::String="")
    cindex = MOI.addconstraint!(m.moibackend, moi_function_and_set(c)...)
    cref = ConstraintRef(m, cindex)
    if !isempty(name)
        setname(cref, name)
    end
    return cref
end

include("variables.jl")

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

function optimizerindex(v::VariableRef)
    if mode(v.m) == Direct
        return index(v)
    else
        @assert caching_optimizer(v.m).state == MOIU.AttachedOptimizer
        return caching_optimizer(v.m).model_to_optimizer_map[index(v)]
    end
end

function optimizerindex(cr::ConstraintRef{Model})
    if mode(cr.m) == Direct
        return index(cr)
    else
        @assert caching_optimizer(cr.m).state == MOIU.AttachedOptimizer
        return caching_optimizer(cr.m).model_to_optimizer_map[index(cr)]
    end
end

index(cr::ConstraintRef) = cr.index

function hasresultdual(m::Model, REF::Type{<:ConstraintRef{Model, T}}) where {T <: MOICON}
    MOI.canget(m, MOI.ConstraintDual(), REF)
end

"""
    resultdual(cr::ConstraintRef)

Get the dual value of this constraint in the result returned by a solver.
Use `hasresultdual` to check if a result exists before asking for values.
Replaces `getdual` for most use cases.
"""
function resultdual(cr::ConstraintRef{Model, <:MOICON})
    MOI.get(cr.m, MOI.ConstraintDual(), cr)
end

"""
    name(v::ConstraintRef)

Get a constraint's name.
"""
name(cr::ConstraintRef{Model,<:MOICON}) = MOI.get(cr.m, MOI.ConstraintName(), cr)

setname(cr::ConstraintRef{Model,<:MOICON}, s::String) = MOI.set!(cr.m, MOI.ConstraintName(), cr, s)

"""
    canget(m::JuMP.Model, attr::MathOptInterface.AbstractModelAttribute)::Bool

Return `true` if one may query the attribute `attr` from the model's MOI backend.
false if not.
"""
MOI.canget(m::Model, attr::MOI.AbstractModelAttribute) = MOI.canget(m.moibackend, attr)
MOI.canget(m::Model, attr::MOI.AbstractVariableAttribute, ::Type{VariableRef}) = MOI.canget(m.moibackend, attr, MOIVAR)
MOI.canget(m::Model, attr::MOI.AbstractConstraintAttribute, ::Type{ConstraintRef{Model,T}}) where {T <: MOICON} = MOI.canget(m.moibackend, attr, T)

"""
    get(m::JuMP.Model, attr::MathOptInterface.AbstractModelAttribute)

Return the value of the attribute `attr` from model's MOI backend.
"""
MOI.get(m::Model, attr::MOI.AbstractModelAttribute) = MOI.get(m.moibackend, attr)
function MOI.get(m::Model, attr::MOI.AbstractVariableAttribute, v::VariableRef)
    @assert m === v.m
    MOI.get(m.moibackend, attr, index(v))
end
function MOI.get(m::Model, attr::MOI.AbstractConstraintAttribute, cr::ConstraintRef)
    @assert m === cr.m
    MOI.get(m.moibackend, attr, index(cr))
end

MOI.set!(m::Model, attr::MOI.AbstractModelAttribute, value) = MOI.set!(m.moibackend, attr, value)
function MOI.set!(m::Model, attr::MOI.AbstractVariableAttribute, v::VariableRef, value)
    @assert m === v.m
    MOI.set!(m.moibackend, attr, index(v), value)
end
function MOI.set!(m::Model, attr::MOI.AbstractConstraintAttribute, cr::ConstraintRef, value)
    @assert m === cr.m
    MOI.set!(m.moibackend, attr, index(cr), value)
end

###############################################################################
# GenericAffineExpression, AffExpr, AffExprConstraint
include("affexpr.jl")



###############################################################################
# GenericQuadExpr, QuadExpr
# GenericQuadConstraint, QuadConstraint
include("quadexpr.jl")

##########################################################################
# SOSConstraint  (special ordered set constraints)
# include("sos.jl")

include("sets.jl")

##########################################################################
# SDConstraint
include("sd.jl")

# handle dictionary of variables
function registervar(m::AbstractModel, varname::Symbol, value)
    registerobject(m, varname, value, "A variable or constraint named $varname is already attached to this model. If creating variables programmatically, use the anonymous variable syntax x = @variable(m, [1:N], ...).")
end
registervar(m::AbstractModel, varname, value) = error("Invalid variable name $varname")

function registercon(m::AbstractModel, conname::Symbol, value)
    registerobject(m, conname, value, "A variable or constraint named $conname is already attached to this model. If creating constraints programmatically, use the anonymous constraint syntax con = @constraint(m, ...).")
end
registercon(m::AbstractModel, conname, value) = error("Invalid constraint name $conname")

# This function needs to be implemented by all `AbstractModel`s
object_dictionary(m::Model) = m.objdict

function registerobject(m::AbstractModel, name::Symbol, value, errorstring::String)
    objdict = object_dictionary(m)
    if haskey(objdict, name)
        error(errorstring)
        objdict[name] = nothing
    else
        objdict[name] = value
    end
    return value
end


"""
    Base.getindex(m::JuMP.AbstractModel, name::Symbol)

To allow easy accessing of JuMP tVariables and Constraints via `[]` syntax.
Returns the variable, or group of variables, or constraint, or group of constraints, of the given name which were added to the model. This errors if multiple variables or constraints share the same name.
"""
function Base.getindex(m::JuMP.AbstractModel, name::Symbol)
    objdict = object_dictionary(m)
    if !haskey(objdict, name)
        throw(KeyError("No object with name $name"))
    elseif objdict[name] === nothing
        error("There are multiple variables and/or constraints named $name that are already attached to this model. If creating variables programmatically, use the anonymous variable syntax x = @variable(m, [1:N], ...). If creating constraints programmatically, use the anonymous constraint syntax con = @constraint(m, ...).")
    else
        return objdict[name]
    end
end

"""
    Base.setindex!(m::JuMP.Model, value, name::Symbol)

stores the object `value` in the model `m` using so that it can be accessed via `getindex`.  Can be called with `[]` syntax.
"""
function Base.setindex!(m::JuMP.Model, value, name::Symbol)
    # if haskey(m.objdict, name)
    #     warn("Overwriting the object $name stored in the model. Consider using anonymous variables and constraints instead")
    # end
    m.objdict[name] = value
end

"""
    operator_warn(m::AbstractModel)

Everytime two expressions are summed not using `destructive_add!` and one of
the two expressions have more than 50 terms, this function is called on the model.

## Notes for extensions

By default this method does nothing so every new model type must implement this
function in order to print a warning.
"""
function operator_warn(::AbstractModel) end
function operator_warn(m::Model)
    m.operator_counter += 1
    if m.operator_counter > 20000
        Base.warn_once("The addition operator has been used on JuMP expressions a large number of times. This warning is safe to ignore but may indicate that model generation is slower than necessary. For performance reasons, you should not add expressions in a loop. Instead of x += y, use append!(x,y) to modify x in place. If y is a single variable, you may also use push!(x, coef, y) in place of x += coef*y.")
    end
end
function operator_warn(lhs::GenericAffExpr,rhs::GenericAffExpr)
    if length(linearterms(lhs)) > 50 || length(linearterms(rhs)) > 50
        if length(linearterms(lhs)) > 1
            operator_warn(owner_model(first(linearterms(lhs))[2]))
        end
    end
    return
end
operator_warn(lhs,rhs) = nothing

##########################################################################
# Types used in the nonlinear code
struct NonlinearExpression
    m::Model
    index::Int
end

struct NonlinearParameter <: AbstractJuMPScalar
    m::Model
    index::Int
end


##########################################################################
# Behavior that's uniform across all JuMP "scalar" objects

# TODO why do we need this?
const JuMPTypes = Union{AbstractJuMPScalar,
                        NonlinearExpression}
                    #    Norm,
                    #    QuadExpr,
                    #    SOCExpr}
const JuMPScalars = Union{Number,JuMPTypes}

# would really want to do this on ::Type{T}, but doesn't work on v0.4
Base.eltype(::T) where {T<:JuMPTypes} = T
Base.size(::JuMPTypes) = ()
Base.size(x::JuMPTypes,d::Int) = 1
Base.ndims(::JuMPTypes) = 0


##########################################################################
include("containers.jl")
include("operators.jl")
# include("writers.jl")
include("macros.jl")
include("optimizerinterface.jl")
# include("callbacks.jl")
include("nlp.jl")
include("print.jl")


##########################################################################
end

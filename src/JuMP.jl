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
using ForwardDiff
include("Derivatives/Derivatives.jl")
using .Derivatives

export
# Objects
    Model, VariableRef, Norm, AffExpr, QuadExpr, SOCExpr,
    # LinearConstraint, QuadConstraint, SDConstraint, SOCConstraint,
    NonlinearConstraint,
    ConstraintRef,
# Cones
    PSDCone,
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

###############################################################################
# Model

# Model has three modes:
# 1) Automatic: moibackend field holds a CachingOptimizer in Automatic mode.
# 2) Manual: moibackend field holds a CachingOptimizer in Manual mode.
# 3) Direct: moibackend field holds an AbstractOptimizer. No extra copy of the model is stored. The moibackend must support adddconstraint! etc.
# Methods to interact with the CachingOptimizer are defined in solverinterface.jl.
@enum ModelMode Automatic Manual Direct

abstract type AbstractModel end
mutable struct Model <: AbstractModel

    # special variablewise properties that we keep track of:
    # lower bound, upper bound, fixed, integrality, binary
    variabletolowerbound::Dict{MOIVAR,MOILB}
    variabletoupperbound::Dict{MOIVAR,MOIUB}
    variabletofix::Dict{MOIVAR,MOIFIX}
    variabletointegrality::Dict{MOIVAR,MOIINT}
    variabletozeroone::Dict{MOIVAR,MOIBIN}

    customnames::Vector

    # # Variable cones of the form, e.g. (:SDP, 1:9)
    # varCones::Vector{Tuple{Symbol,Any}}

    # Solution data
    objbound
    objval
    # colVal::Vector{Float64}
    # redCosts::Vector{Float64}
    # linconstrDuals::Vector{Float64}
    # conicconstrDuals::Vector{Float64}
    # constr_to_row::Vector{Vector{Int}}
    # # Vector of the same length as sdpconstr.
    # # sdpconstrSym[c] is the list of pairs (i,j), i > j
    # # such that a symmetry-enforcing constraint has been created
    # # between sdpconstr[c].terms[i,j] and sdpconstr[c].terms[j,i]
    # sdpconstrSym::Vector{Vector{Tuple{Int,Int}}}
    moibackend::Union{MOI.AbstractOptimizer,MOIU.CachingOptimizer}
    # callbacks
    callbacks
    # lazycallback
    # cutcallback
    # heurcallback

    # hook into a solve call...function of the form f(m::Model; kwargs...),
    # where kwargs get passed along to subsequent solve calls
    optimizehook
    # # ditto for a print hook
    # printhook


    # # storage vector for merging duplicate terms
    # indexedVector::IndexedVector{Float64}

    nlpdata#::NLPData

    objdict::Dict{Symbol,Any} # dictionary from variable and constraint names to objects

    operator_counter::Int # number of times we add large expressions

    # Extension dictionary - e.g. for robust
    # Extensions should define a type to hold information particular to
    # their functionality, and store an instance of the type in this
    # dictionary keyed on an extension-specific symbol
    ext::Dict{Symbol,Any}
    # Default constructor
    function Model(; mode::ModelMode=Automatic, backend=nothing, optimizer=nothing)
        m = new()
        # TODO make pretty
        m.variabletolowerbound = Dict{MOIVAR,MOILB}()
        m.variabletoupperbound = Dict{MOIVAR,MOIUB}()
        m.variabletofix = Dict{MOIVAR,MOIFIX}()
        m.variabletointegrality = Dict{MOIVAR,MOIINT}()
        m.variabletozeroone = Dict{MOIVAR,MOIBIN}()
        m.customnames = VariableRef[]
        m.objbound = 0.0
        m.objval = 0.0
        if backend != nothing
            # TODO: It would make more sense to not force users to specify Direct mode if they also provide a backend.
            @assert mode == Direct
            @assert optimizer === nothing
            @assert MOI.isempty(backend)
            m.moibackend = backend
        else
            @assert mode != Direct
            m.moibackend = MOIU.CachingOptimizer(MOIU.UniversalFallback(JuMPMOIModel{Float64}()), mode == Automatic ? MOIU.Automatic : MOIU.Manual)
            if optimizer !== nothing
                MOIU.resetoptimizer!(m, optimizer)
            end
        end
        m.callbacks = Any[]
        m.optimizehook = nothing
        # m.printhook = nothing
        m.nlpdata = nothing
        m.objdict = Dict{Symbol,Any}()
        m.operator_counter = 0
        m.ext = Dict{Symbol,Any}()

        return m
    end
end

# Getters/setters

function mode(m::Model)
    if !(m.moibackend isa MOIU.CachingOptimizer)
        return Direct
    elseif m.moibackend.mode == MOIU.Automatic
        return Automatic
    else
        return Manual
    end
end

# temporary name
numvar(m::Model) = MOI.get(m, MOI.NumberOfVariables())

"""
    numnlconstr(m::Model)

returns the number of nonlinear constraints associated with the `Model m`
"""
numnlconstr(m::Model) = m.nlpdata !== nothing ? length(m.nlpdata.nlconstr) : 0


"""
    objectivebound(m::Model)

Return the best known bound on the optimal objective value after a call to `solve`.
"""
objectivebound(m::Model) = MOI.get(m, MOI.ObjectiveBound())

"""
    objectivevalue(m::Model)

Return the objective value after a call to `solve`.
"""
objectivevalue(m::Model) = MOI.get(m, MOI.ObjectiveValue())

"""
    objectivesense(m::Model)

Return the objective sense, `:Min`, `:Max`, or `:Feasibility`.
"""
function objectivesense(m::Model)
    moisense = MOI.get(m, MOI.ObjectiveSense())
    if moisense == MOI.MinSense
        return :Min
    elseif moisense == MOI.MaxSense
        return :Max
    else
        @assert moisense == MOI.FeasibilitySense
        return :Feasibility
    end
end

terminationstatus(m::Model) = MOI.get(m, MOI.TerminationStatus())
primalstatus(m::Model) = MOI.get(m, MOI.PrimalStatus())
dualstatus(m::Model) = MOI.get(m, MOI.DualStatus())

# TODO: Implement Base.copy.

setoptimizehook(m::Model, f) = (m.optimizehook = f)
setprinthook(m::Model, f) = (m.printhook = f)


#############################################################################
# AbstractConstraint
# Abstract base type for all constraint types
abstract type AbstractConstraint end
# Abstract base type for all scalar types
# In JuMP, used only for VariableRef. Useful primarily for extensions
abstract type AbstractJuMPScalar end

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

Base.copy(v::VariableRef, new_model::Model) = VariableRef(new_model, v.col)
Base.copy(x::Nothing, new_model::Model) = nothing
Base.copy(v::AbstractArray{VariableRef}, new_model::Model) = (var -> VariableRef(new_model, var.col)).(v)


function optimizerindex(v::VariableRef)
    if mode(v.m) == Direct
        return index(v)
    else
        @assert v.m.moibackend.state == MOIU.AttachedOptimizer
        return v.m.moibackend.model_to_optimizer_map[index(v)]
    end
end

function optimizerindex(cr::ConstraintRef{Model})
    if mode(cr.m) == Direct
        return index(cr)
    else
        @assert cr.m.moibackend.state == MOIU.AttachedOptimizer
        return cr.m.moibackend.model_to_optimizer_map[index(cr)]
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

##########################################################################
# SDConstraint
include("sd.jl")

# handle dictionary of variables
function registervar(m::Model, varname::Symbol, value)
    registerobject(m, varname, value, "A variable or constraint named $varname is already attached to this model. If creating variables programmatically, use the anonymous variable syntax x = @variable(m, [1:N], ...).")
end
registervar(m::Model, varname, value) = error("Invalid variable name $varname")

function registercon(m::Model, conname::Symbol, value)
    registerobject(m, conname, value, "A variable or constraint named $conname is already attached to this model. If creating constraints programmatically, use the anonymous constraint syntax con = @constraint(m, ...).")
end
registercon(m::Model, conname, value) = error("Invalid constraint name $conname")

function registerobject(m::Model, name::Symbol, value, errorstring::String)
    if haskey(m.objdict, name)
        error(errorstring)
        m.objdict[name] = nothing
    else
        m.objdict[name] = value
    end
    return value
end


"""
    Base.getindex(m::JuMP.Model, name::Symbol)

To allow easy accessing of JuMP tVariables and Constraints via `[]` syntax.
Returns the variable, or group of variables, or constraint, or group of constraints, of the given name which were added to the model. This errors if multiple variables or constraints share the same name.
"""
function Base.getindex(m::JuMP.Model, name::Symbol)
    if !haskey(m.objdict, name)
        throw(KeyError("No object with name $name"))
    elseif m.objdict[name] === nothing
        error("There are multiple variables and/or constraints named $name that are already attached to this model. If creating variables programmatically, use the anonymous variable syntax x = @variable(m, [1:N], ...). If creating constraints programmatically, use the anonymous constraint syntax con = @constraint(m, ...).")
    else
        return m.objdict[name]
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

# usage warnings
function operator_warn(lhs::GenericAffExpr,rhs::GenericAffExpr)
    if length(lhs.vars) > 50 || length(rhs.vars) > 50
        if length(lhs.vars) > 1
            m = lhs.vars[1].m
            m.operator_counter += 1
            if m.operator_counter > 20000
                Base.warn_once("The addition operator has been used on JuMP expressions a large number of times. This warning is safe to ignore but may indicate that model generation is slower than necessary. For performance reasons, you should not add expressions in a loop. Instead of x += y, use append!(x,y) to modify x in place. If y is a single variable, you may also use push!(x, coef, y) in place of x += coef*y.")
            end
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

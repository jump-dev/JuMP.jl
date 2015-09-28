#  Copyright 2015, Iain Dunning, Joey Huchette, Miles Lubin, and contributors
#  This Source Code Form is subject to the terms of the Mozilla Public
#  License, v. 2.0. If a copy of the MPL was not distributed with this
#  file, You can obtain one at http://mozilla.org/MPL/2.0/.
#############################################################################
# JuMP
# An algebraic modeling language for Julia
# See http://github.com/JuliaOpt/JuMP.jl
#############################################################################

isdefined(Base, :__precompile__) && __precompile__()

module JuMP

importall Base.Operators

import MathProgBase

using ReverseDiffSparse, Calculus
import ArrayViews
const subarr = ArrayViews.view

using Compat

export
# Objects
    Model, Variable, Norm, AffExpr, QuadExpr, SOCExpr, AbstractJuMPScalar,
    LinearConstraint, QuadConstraint, SDPConstraint, SOCConstraint,
    ConstraintRef, LinConstrRef,
    JuMPNLPEvaluator,
# Functions
    # Model related
	getCost, getObjectiveValue, getObjective,
    getObjectiveSense, setObjectiveSense, setSolver,
    writeLP, writeMPS, setObjective,
    addConstraint, addSOS1, addSOS2, solve,
    getInternalModel, buildInternalModel, setSolveHook, setPrintHook,
    # Variable
    setName, getName, setLower, setUpper, getLower, getUpper,
    getValue, setValue, getDual, setCategory, getCategory,
    getVar,
    getLinearIndex,
    # Expressions and constraints
    affToStr, quadToStr, exprToStr, conToStr, chgConstrRHS,

# Macros and support functions
    @addConstraint, @addConstraints, @addSDPConstraint,
    @LinearConstraint, @LinearConstraints, @QuadConstraint, @QuadConstraints,
    @SOCConstraint, @SOCConstraints,
    @defVar, @defConstrRef, @setObjective, addToExpression, @defExpr,
    @setNLObjective, @addNLConstraint, @addNLConstraints,
    @defNLExpr

include("JuMPContainer.jl")
#include("JuMPArray.jl")
include("utils.jl")

###############################################################################
# Model class
# Keeps track of all model and column info
type Model
    obj#::QuadExpr
    objSense::Symbol

    linconstr#::Vector{LinearConstraint}
    quadconstr
    sosconstr
    socconstr
    sdpconstr

    # Column data
    numCols::Int
    colNames::Vector{UTF8String}
    colNamesIJulia::Vector{UTF8String}
    colLower::Vector{Float64}
    colUpper::Vector{Float64}
    colCat::Vector{Symbol}

    # Variable cones of the form, e.g. (:SDP, 1:9)
    varCones::Vector{@compat Tuple{Symbol,Any}}

    # Solution data
    objVal
    colVal::Vector{Float64}
    redCosts::Vector{Float64}
    linconstrDuals::Vector{Float64}
    # internal solver model object
    internalModel
    # Solver+option object from MathProgBase
    solver::MathProgBase.AbstractMathProgSolver
    internalModelLoaded::Bool
    # callbacks
    callbacks
    # lazycallback
    # cutcallback
    # heurcallback

    # hook into a solve call...function of the form f(m::Model; kwargs...),
    # where kwargs get passed along to subsequent solve calls
    solvehook
    # ditto for a print hook
    printhook

    # List of JuMPContainer{Variables} associated with model
    dictList::Vector

    # storage vector for merging duplicate terms
    indexedVector::IndexedVector{Float64}

    nlpdata#::NLPData

    varDict::Dict{Symbol,Any} # dictionary from variable names to variable objects
    varData::ObjectIdDict

    getvalue_counter::Int # number of times we call getValue on a JuMPContainer, so that we can print out a warning
    operator_counter::Int # number of times we add large expressions

    # Extension dictionary - e.g. for robust
    # Extensions should define a type to hold information particular to
    # their functionality, and store an instance of the type in this
    # dictionary keyed on an extension-specific symbol
    ext::Dict{Symbol,Any}
end

# dummy solver
type UnsetSolver <: MathProgBase.AbstractMathProgSolver
end

# Default constructor
function Model(;solver=UnsetSolver())
    if !isa(solver,MathProgBase.AbstractMathProgSolver)
        error("solver argument ($solver) must be an AbstractMathProgSolver")
    end
    Model(zero(QuadExpr),              # obj
          :Min,                        # objSense
          LinearConstraint[],          # linconstr
          QuadConstraint[],            # quadconstr
          SOSConstraint[],             # sosconstr
          SOCConstraint[],             # socconstr
          SDPConstraint[],             # sdpconstr
          0,                           # numCols
          UTF8String[],                # colNames
          UTF8String[],                # colNamesIJulia
          Float64[],                   # colLower
          Float64[],                   # colUpper
          Symbol[],                    # colCat
          Vector{@compat Tuple{Symbol,Any}}[], # varCones
          0,                           # objVal
          Float64[],                   # colVal
          Float64[],                   # redCosts
          Float64[],                   # linconstrDuals
          nothing,                     # internalModel
          solver,                      # solver
          false,                       # internalModelLoaded
          Any[],                       # callbacks
          nothing,                     # solvehook
          nothing,                     # printhook
          Any[],                       # dictList
          IndexedVector(Float64,0),    # indexedVector
          nothing,                     # nlpdata
          Dict{Symbol,Any}(),          # varDict
          ObjectIdDict(),              # varData
          0,                           # getvalue_counter
          0,                           # operator_counter
          Dict{Symbol,Any}(),          # ext
    )
end

# Getters/setters
MathProgBase.numvar(m::Model) = m.numCols
MathProgBase.numlinconstr(m::Model) = length(m.linconstr)
MathProgBase.numquadconstr(m::Model) = length(m.quadconstr)
function MathProgBase.numconstr(m::Model)
    c = length(m.linconstr) + length(m.quadconstr) + length(m.sosconstr)
    if m.nlpdata !== nothing
        c += length(m.nlpdata.nlconstr)
    end
    return c
end
function MathProgBase.getsolvetime(m::Model)
    if !m.internalModelLoaded
        error("Model not solved")
    elseif method_exists(MathProgBase.getsolvetime, (typeof(getInternalModel(m)), ))
        return MathProgBase.getsolvetime(getInternalModel(m))
    else
        error("Solve time not implemented for $(typeof(m.solver))")
    end
end
function MathProgBase.getnodecount(m::Model)
    if !m.internalModelLoaded
        error("Model not solved")
    elseif method_exists(MathProgBase.getnodecount, (typeof(getInternalModel(m)), ))
        return MathProgBase.getnodecount(getInternalModel(m))
    else
        error("Node count not implemented for $(typeof(m.solver)).")
    end
end

@Base.deprecate getNumVars(m::Model) MathProgBase.numvar(m)
@Base.deprecate getNumConstraints(m::Model) MathProgBase.numlinconstr(m)

getCost(m::Model) = m.redCosts
getObjective(m::Model) = m.obj
getObjectiveValue(m::Model) = m.objVal
getObjectiveSense(m::Model) = m.objSense
function setObjectiveSense(m::Model, newSense::Symbol)
    if (newSense != :Max && newSense != :Min)
        error("Model sense must be :Max or :Min")
    end
    m.objSense = newSense
end
setObjective(m::Model, something::Any) =
    error("in setObjective: needs three arguments: model, objective sense (:Max or :Min), and expression.")

setObjective(::Model, ::Symbol, x::Array) =
    error("in setObjective: array of size $(size(x)) passed as objective; only scalar objectives are allowed")

function setSolver(m::Model, solver::MathProgBase.AbstractMathProgSolver)
    m.solver = solver
    m.internalModel = nothing
    m.internalModelLoaded = false
end
# Deep copy the model
function Base.copy(source::Model)

    dest = Model()
    dest.solver = source.solver  # The two models are linked by this

    # Objective
    dest.obj = copy(source.obj, dest)
    dest.objSense = source.objSense

    # Constraints
    dest.linconstr  = map(c->copy(c, dest), source.linconstr)
    dest.quadconstr = map(c->copy(c, dest), source.quadconstr)
    dest.sosconstr  = map(c->copy(c, dest), source.sosconstr)
    dest.sdpconstr  = map(c->copy(c, dest), source.sdpconstr)

    # Variables
    dest.numCols = source.numCols
    dest.colNames = source.colNames[:]
    dest.colNamesIJulia = source.colNamesIJulia[:]
    dest.colLower = source.colLower[:]
    dest.colUpper = source.colUpper[:]
    dest.colCat = source.colCat[:]

    # varCones
    dest.varCones = copy(source.varCones)

    # callbacks and hooks
    if !isempty(source.callbacks)
        error("Copying callbacks is not supported")
    end
    if source.solvehook !== nothing
        dest.solvehook = copy(source.solvehook)
    end
    if source.printhook !== nothing
        dest.printhook = copy(source.printhook)
    end

    # variable/extension dicts
    if !isempty(source.ext)
        error("Copying of extension dictionaries is not currently supported")
    end
    dest.varDict = Dict{Symbol,Any}()
    for (symb,v) in source.varDict
        dest.varDict[symb] = copy(v, dest)
    end

    # varData---possibly shouldn't copy

    if source.nlpdata !== nothing
        dest.nlpdata = copy(source.nlpdata)
    end

    return dest
end

getInternalModel(m::Model) = m.internalModel

setSolveHook(m::Model, f) = (m.solvehook = f)
setPrintHook(m::Model, f) = (m.printhook = f)


#############################################################################
# JuMPConstraint
# Abstract base type for all constraint types
abstract JuMPConstraint
# Abstract base type for all scalar types
# In JuMP, used only for Variable. Useful primarily for extensions
abstract AbstractJuMPScalar <: ReverseDiffSparse.Placeholder


#############################################################################
# Variable class
# Doesn't actually do much, just a pointer back to the model
immutable Variable <: AbstractJuMPScalar
    m::Model
    col::Int
end

getLinearIndex(x::Variable) = x.col
ReverseDiffSparse.getplaceindex(x::Variable) = getLinearIndex(x)
Base.isequal(x::Variable,y::Variable) = (x.col == y.col) && (x.m === y.m)

Variable(m::Model, lower, upper, cat::Symbol, name::UTF8String="", value::Number=NaN) =
    error("Attempt to create scalar Variable with lower bound of type $(typeof(lower)) and upper bound of type $(typeof(upper)). Bounds must be scalars in Variable constructor.")

function Variable(m::Model,lower::Number,upper::Number,cat::Symbol,name::UTF8String="",value::Number=NaN)
    m.numCols += 1
    push!(m.colNames, name)
    push!(m.colNamesIJulia, name)
    push!(m.colLower, convert(Float64,lower))
    push!(m.colUpper, convert(Float64,upper))
    push!(m.colCat, cat)
    push!(m.colVal,value)
    if cat == :Fixed
        @assert lower == upper
        m.colVal[end] = lower
    end
    if m.internalModelLoaded
        if method_exists(MathProgBase.addvar!, (typeof(m.internalModel),Vector{Int},Vector{Float64},Float64,Float64,Float64))
            MathProgBase.addvar!(m.internalModel,float(lower),float(upper),0.0)
        else
            Base.warn_once("Solver does not appear to support adding variables to an existing model. Hot-start is disabled.")
            m.internalModelLoaded = false
        end
    end
    return Variable(m, m.numCols)
end

# Name setter/getters
function setName(v::Variable,n::AbstractString)
    v.m.colNames[v.col] = n
    v.m.colNamesIJulia[v.col] = n
end
getName(m::Model, col) = var_str(REPLMode, m, col)
getName(v::Variable) = var_str(REPLMode, v.m, v.col)

# Bound setter/getters
function setLower(v::Variable,lower::Number)
    v.m.colCat[v.col] == :Fixed && error("use setValue for changing the value of a fixed variable")
    v.m.colLower[v.col] = lower
end
function setUpper(v::Variable,upper::Number)
    v.m.colCat[v.col] == :Fixed && error("use setValue for changing the value of a fixed variable")
    v.m.colUpper[v.col] = upper
end
getLower(v::Variable) = v.m.colLower[v.col]
getUpper(v::Variable) = v.m.colUpper[v.col]

# Value setter/getter
function setValue(v::Variable, val::Number)
    v.m.colVal[v.col] = val
    if v.m.colCat[v.col] == :Fixed
        v.m.colLower[v.col] = val
        v.m.colUpper[v.col] = val
    end
end

function setValue(set::Array{Variable}, val::Array)
    promote_shape(size(set), size(val)) # Check dimensions match
    for I in eachindex(set)
        setValue(set[I], val[I])
    end
    nothing
end

# internal method that doesn't print a warning if the value is NaN
_getValue(v::Variable) = v.m.colVal[v.col]

function getValue(v::Variable)
    ret = _getValue(v)
    if isnan(ret)
        Base.warn("Variable value not defined for $(getName(v)). Check that the model was properly solved.")
    end
    ret
end

function getValue(arr::Array{Variable})
    ret = similar(arr, Float64)
    # return immediately for empty array
    if isempty(ret)
        return ret
    end
    # warnedyet is set to true if we've already warned for a component of a JuMPContainer
    warnedyet = false
    m = first(arr).m
    # whether this was constructed via @defVar, essentially
    registered = haskey(m.varData, arr)
    for I in eachindex(arr)
        v = arr[I]
        value = _getValue(v)
        ret[I] = value
        if !warnedyet && isnan(value)
            if registered
                Base.warn("Variable value not defined for component of $(m.varData[arr].name). Check that the model was properly solved.")
                warnedyet = true
            else
                Base.warn("Variable value not defined for $(m.colNames[v.col]). Check that the model was properly solved.")
            end
        end
    end
    # Copy printing data from @defVar for Array{Variable} to corresponding Array{Float64} of values
    if registered
        m.varData[ret] = m.varData[arr]
    end
    ret
end


# Dual value (reduced cost) getter
function getDual(v::Variable)
    if length(v.m.redCosts) < MathProgBase.numvar(v.m)
        error("Variable bound duals (reduced costs) not available. Check that the model was properly solved and no integer variables are present.")
    end
    return v.m.redCosts[v.col]
end

const var_cats = [:Cont, :Int, :Bin, :SemiCont, :SemiInt]
function setCategory(v::Variable, cat::Symbol)
    cat in var_cats || error("Unrecognized variable category $cat. Should be one of:\n    $var_cats")
    v.m.colCat[v.col] = cat
end

getCategory(v::Variable) = v.m.colCat[v.col]

Base.zero(::Type{Variable}) = AffExpr(Variable[],Float64[],0.0)
Base.zero(::Variable) = zero(Variable)
Base.one(::Type{Variable}) = AffExpr(Variable[],Float64[],1.0)
Base.one(::Variable) = one(Variable)

function verify_ownership(m::Model, vec::Vector{Variable})
    n = length(vec)
    @inbounds for i in 1:n
        vec[i].m !== m && return false
    end
    return true
end

Base.copy(v::Variable, new_model::Model) = Variable(new_model, v.col)
function Base.copy(v::Array{Variable}, new_model::Model)
    ret = similar(v, Variable, size(v))
    for I in eachindex(v)
        ret[I] = Variable(new_model, v[I].col)
    end
    ret
end

# Copy methods for variable containers
Base.copy(d::JuMPContainer) = map(copy, d)
Base.copy(d::JuMPContainer, new_model::Model) = map(x -> copy(x, new_model), d)

###############################################################################
# GenericAffineExpression, AffExpr
# GenericRangeConstraint, LinearConstraint
include("affexpr.jl")

###############################################################################
# GenericQuadExpr, QuadExpr
# GenericQuadConstraint, QuadConstraint
include("quadexpr.jl")

##########################################################################
# GenericNorm, Norm
# GenericNormExpr. GenericSOCExpr, SOCExpr
# GenericSOCConstraint, SOCConstraint
include("norms.jl")

##########################################################################
# SOSConstraint  (special ordered set constraints)
include("sos.jl")

##########################################################################
# SDPConstraint is a (dual) semidefinite constraint of the form
# ∑ cᵢ Xᵢ ≥ D, where D is a n×n symmetric data matrix, cᵢ are
# scalars, and Xᵢ are n×n symmetric variable matrices. The inequality
# is taken w.r.t. the psd partial order.
type SDPConstraint <: JuMPConstraint
    terms
end

# Special-case X ≥ 0, which is often convenient
function SDPConstraint(lhs::Matrix, rhs::Number)
    rhs == 0 || error("Cannot construct a semidefinite constraint with nonzero scalar bound $rhs")
    SDPConstraint(lhs)
end

function addConstraint(m::Model, c::SDPConstraint)
    push!(m.sdpconstr,c)
    m.internalModelLoaded = false
    ConstraintRef{SDPConstraint}(m,length(m.sdpconstr))
end

# helper method for mapping going on below
Base.copy(x::Number, new_model::Model) = copy(x)

Base.copy(c::SDPConstraint, new_model::Model) =
    SDPConstraint(map(t -> copy(t, new_model), c.terms))


##########################################################################
# ConstraintRef
# Reference to a constraint for retrieving solution info
immutable ConstraintRef{T<:JuMPConstraint}
    m::Model
    idx::Int
end

typealias LinConstrRef ConstraintRef{LinearConstraint}

getLinearIndex(x::ConstraintRef) = x.idx

function getDual(c::ConstraintRef{LinearConstraint})
    if length(c.m.linconstrDuals) != MathProgBase.numlinconstr(c.m)
        error("Dual solution not available. Check that the model was properly solved and no integer variables are present.")
    end
    return c.m.linconstrDuals[c.idx]
end

function chgConstrRHS(c::ConstraintRef{LinearConstraint}, rhs::Number)
    constr = c.m.linconstr[c.idx]
    sen = sense(constr)
    if sen == :range
        error("Modifying range constraints is currently unsupported.")
    elseif sen == :(==)
        constr.lb = float(rhs)
        constr.ub = float(rhs)
    elseif sen == :>=
        constr.lb = float(rhs)
    else
        @assert sen == :<=
        constr.ub = float(rhs)
    end
end

Variable(m::Model,lower::Number,upper::Number,cat::Symbol,objcoef::Number,
    constraints::JuMPArray,coefficients::Vector{Float64}, name::AbstractString="", value::Number=NaN) =
    Variable(m, lower, upper, cat, objcoef, constraints.innerArray, coefficients, name, value)

# add variable to existing constraints
function Variable(m::Model,lower::Number,upper::Number,cat::Symbol,objcoef::Number,
    constraints::Vector,coefficients::Vector{Float64}, name::AbstractString="", value::Number=NaN)
    for c in constraints
        if !isa(c,ConstraintRef{LinearConstraint})
            error("Unexpected constraint of type $(typeof(c)). Column-wise modeling only supported for linear constraints")
        end
    end
    @assert cat != :Fixed || (lower == upper)
    m.numCols += 1
    push!(m.colNames, name)
    push!(m.colNamesIJulia, name)
    push!(m.colLower, convert(Float64,lower))
    push!(m.colUpper, convert(Float64,upper))
    push!(m.colCat, cat)
    push!(m.colVal,value)
    if cat == :Fixed
        @assert lower == upper
        m.colVal[end] = lower
    end
    v = Variable(m,m.numCols)
    # add to existing constraints
    @assert length(constraints) == length(coefficients)
    for i in 1:length(constraints)
        c::LinearConstraint = m.linconstr[constraints[i].idx]
        coef = coefficients[i]
        push!(c.terms.vars,v)
        push!(c.terms.coeffs,coef)
    end
    push!(m.obj.aff.vars, v)
    push!(m.obj.aff.coeffs,objcoef)

    if m.internalModelLoaded
        if method_exists(MathProgBase.addvar!, (typeof(m.internalModel),Vector{Int},Vector{Float64},Float64,Float64,Float64))
            MathProgBase.addvar!(m.internalModel,Int[c.idx for c in constraints],coefficients,float(lower),float(upper),float(objcoef))
        else
            Base.warn_once("Solver does not appear to support adding variables to an existing model. Hot-start is disabled.")
            m.internalModelLoaded = false
        end
    end

    return v
end

# handle dictionary of variables
function registervar(m::Model, varname::Symbol, value)
    if haskey(m.varDict, varname)
        m.varDict[varname] = nothing # indicate duplicate variable
    else
        m.varDict[varname] = value
    end
    return value
end
registervar(m::Model, varname, value) = value # variable name isn't a simple symbol, ignore

function getVar(m::Model, varname::Symbol)
    if !haskey(m.varDict, varname)
        error("No variable with name $varname")
    elseif m.varDict[varname] === nothing
        error("Multiple variables with name $varname")
    else
        return m.varDict[varname]
    end
end

# usage warnings
function getvalue_warn(x::JuMPContainer{Variable})
    isempty(x) && return
    v = first(values(x))
    m = v.m
    m.getvalue_counter += 1
    if m.getvalue_counter > 400
        Base.warn_once("getValue has been called on a collection of variables a large number of times. For performance reasons, this should be avoided. Instead of getValue(x)[a,b,c], use getValue(x[a,b,c]) to avoid temporary allocations.")
    end
    return
end
getvalue_warn(x::JuMPContainer) = nothing

function operator_warn(lhs::AffExpr,rhs::AffExpr)
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
# Behavior that's uniform across all JuMP "scalar" objects

@compat typealias JuMPTypes Union{AbstractJuMPScalar,
                          Norm,
                          GenericAffExpr,
                          QuadExpr,
                          SOCExpr}
@compat typealias JuMPScalars Union{Number,JuMPTypes}

# would really want to do this on ::Type{T}, but doesn't work on v0.4
Base.eltype{T<:JuMPTypes}(::T) = T
Base.size(::JuMPTypes) = ()
Base.size(x::JuMPTypes,d::Int) = 1
Base.ndims(::JuMPTypes) = 0

##########################################################################
# Operator overloads
include("operators.jl")
# Writers - we support MPS (MILP + QuadObj), LP (MILP)
include("writers.jl")
# Macros - @defVar, sum{}, etc.
include("macros.jl")
# Solvers
include("solvers.jl")
# Callbacks - lazy, cuts, ...
include("callbacks.jl")
# Pretty-printing of JuMP-defined types.
include("print.jl")
# Nonlinear-specific code
include("nlp.jl")

##########################################################################
end

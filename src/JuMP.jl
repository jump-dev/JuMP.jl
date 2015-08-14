#  Copyright 2015, Iain Dunning, Joey Huchette, Miles Lubin, and contributors
#  This Source Code Form is subject to the terms of the Mozilla Public
#  License, v. 2.0. If a copy of the MPL was not distributed with this
#  file, You can obtain one at http://mozilla.org/MPL/2.0/.
#############################################################################
# JuMP
# An algebraic modelling langauge for Julia
# See http://github.com/JuliaOpt/JuMP.jl
#############################################################################

VERSION >= v"0.4.0-dev+6521" && __precompile__()

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
    getObjectiveValue, getObjective,
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
    colNames::Vector{String}
    colNamesIJulia::Vector{String}
    colLower::Vector{Float64}
    colUpper::Vector{Float64}
    colCat::Vector{Symbol}

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
    Model(zero(QuadExpr),:Min,LinearConstraint[],QuadConstraint[],
          SOSConstraint[],SOCConstraint[],SDPConstraint[],
          0,String[],String[],Float64[],Float64[],Symbol[],
          0,Float64[],Float64[],Float64[],nothing,solver,
          false,Any[],nothing,nothing,JuMPContainer[],
          IndexedVector(Float64,0),nothing,Dict{Symbol,Any}(),0,0,Dict{Symbol,Any}())
end

# Getters/setters
MathProgBase.numvar(m::Model) = m.numCols
MathProgBase.numlinconstr(m::Model) = length(m.linconstr)
MathProgBase.numquadconstr(m::Model) = length(m.quadconstr)
function MathProgBase.numconstr(m::Model)
    c = length(m.linconstr) + length(m.quadconstr) + length(m.sosconstr)
    if m.nlpdata != nothing
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

    # callbacks and hooks
    if !isempty(source.callbacks)
        error("Copying callbacks is not supported")
    end
    if source.solvehook != nothing
        dest.solvehook = copy(source.solvehook)
    end
    if source.printhook != nothing
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
    dest.dictList = map(v -> copy(v, dest), source.dictList)

    if source.nlpdata != nothing
        dest.nlpdata = copy(source.nlpdata)
    end

    return dest
end

getInternalModel(m::Model) = m.internalModel

setSolveHook(m::Model, f) = (m.solvehook = f)
setPrintHook(m::Model, f) = (m.printhook = f)

###############################################################################
# Variable class
# Doesn't actually do much, just a pointer back to the model
abstract AbstractJuMPScalar <: ReverseDiffSparse.Placeholder
immutable Variable <: AbstractJuMPScalar
    m::Model
    col::Int
end

getLinearIndex(x::Variable) = x.col
ReverseDiffSparse.getplaceindex(x::Variable) = getLinearIndex(x)
Base.isequal(x::Variable,y::Variable) = isequal(x.col,y.col) && isequal(x.m,y.m)

Variable(m::Model, lower, upper, cat::Symbol, name::String="", value::Number=NaN) =
    error("Attempt to create scalar Variable with lower bound of type $(typeof(lower)) and upper bound of type $(typeof(upper)). Bounds must be scalars in Variable constructor.")

function Variable(m::Model,lower::Number,upper::Number,cat::Symbol,name::String="",value::Number=NaN)
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
function setName(v::Variable,n::String)
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

# internal method that doesn't print a warning if the value is NaN
_getValue(v::Variable) = v.m.colVal[v.col]

function getValue(v::Variable)
    ret = _getValue(v)
    if isnan(ret)
        Base.warn("Variable value not defined for $(getName(v)). Check that the model was properly solved.")
    end
    ret
end

getValue(arr::Array{Variable}) = map(getValue, arr)

# Dual value (reduced cost) getter
function getDual(v::Variable)
    if length(v.m.redCosts) < MathProgBase.numvar(v.m)
        error("Variable bound duals (reduced costs) not available. Check that the model was properly solved and no integer variables are present.")
    end
    return v.m.redCosts[v.col]
end

const var_types = [:Cont, :Int, :Bin, :SemiCont, :SemiInt]
function setCategory(v::Variable, v_type::Symbol)
    v_type in var_types || error("Unrecognized variable type $v_type. Should be one of:\n    $var_types")
    v.m.colCat[v.col] = v_type
end

getCategory(v::Variable) = v.m.colCat[v.col]

Base.zero(::Type{Variable}) = AffExpr(Variable[],Float64[],0.0)
Base.zero(::Variable) = zero(Variable)
Base.one(::Type{Variable}) = AffExpr(Variable[],Float64[],1.0)
Base.one(::Variable) = one(Variable)

verify_ownership(m::Model, vec::Vector{Variable}) = all(v->isequal(v.m,m), vec)

Base.copy(v::Variable, new_model::Model) = Variable(new_model, v.col)

# Copy methods for variable containers
Base.copy(d::JuMPContainer) = map(copy, d)
Base.copy(d::JuMPContainer, new_model::Model) = map(x -> copy(x, new_model), d)
Base.copy{T<:OneIndexedArray}(d::T) = T(map(copy, d.innerArray),
                                       d.name,
                                       d.indexsets,
                                       d.indexexprs)
Base.copy{T<:OneIndexedArray}(d::T, new_model::Model) =
    T(map(v -> copy(v, new_model), d.innerArray),
      d.name,
      d.indexsets,
      d.indexexprs)

###############################################################################
# Generic affine expression class
# Holds a vector of tuples (Var, Coeff)
type GenericAffExpr{CoefType,VarType}
    vars::Array{VarType,1}
    coeffs::Array{CoefType,1}
    constant::CoefType
end

coeftype{CoefType,VarType}(::GenericAffExpr{CoefType,VarType}) = CoefType

Base.zero{CoefType,VarType}(::Type{GenericAffExpr{CoefType,VarType}}) =
    GenericAffExpr{CoefType,VarType}(VarType[],CoefType[],zero(CoefType))
Base.one{CoefType,VarType}(::Type{GenericAffExpr{CoefType,VarType}}) =
    GenericAffExpr{CoefType,VarType}(VarType[],CoefType[],one(CoefType))
Base.zero(a::GenericAffExpr) = zero(typeof(a))
Base.one(a::GenericAffExpr) = one(typeof(a))
Base.start(aff::GenericAffExpr) = 1
Base.done(aff::GenericAffExpr, state::Int) = state > length(aff.vars)
Base.next(aff::GenericAffExpr, state::Int) = ((aff.coeffs[state], aff.vars[state]), state+1)
function Base.in{CoefType,VarType}(x::VarType, aff::GenericAffExpr{CoefType,VarType})
    acc = zero(CoefType)
    for (coef,term) in aff
        if isequal(x, term)
            acc += coef
        end
    end
    return !(acc == zero(CoefType))
end

# More efficient ways to grow an affine expression
# Add a single term to an affine expression
function Base.push!{T,S}(aff::GenericAffExpr{T,S}, new_coeff::T, new_var::S)
    push!(aff.vars, new_var)
    push!(aff.coeffs, new_coeff)
end
# Add an affine expression to an existing affine expression
function Base.append!{T,S}(aff::GenericAffExpr{T,S}, other::GenericAffExpr{T,S})
    append!(aff.vars, other.vars)
    append!(aff.coeffs, other.coeffs)
    aff.constant += other.constant  # Not efficient if CoefType isn't immutable
    aff
end

# Copy an affine expression
Base.copy(aff::GenericAffExpr) = GenericAffExpr(copy(aff.vars),copy(aff.coeffs),copy(aff.constant))

function Base.isequal{T,S}(aff::GenericAffExpr{T,S},other::GenericAffExpr{T,S})
    isequal(aff.constant,other.constant) || return false
    length(aff.vars) == length(other.vars) || return false
    for i in 1:length(aff.vars)
        isequal(aff.vars[i],other.vars[i]) || return false
        isequal(aff.coeffs[i],other.coeffs[i]) || return false
    end
    return true
end

###############################################################################
# Affine expressions, the specific GenericAffExpr used by JuMP
typealias AffExpr GenericAffExpr{Float64,Variable}

AffExpr() = zero(AffExpr)
AffExpr(x::Union(Number,Variable)) = convert(AffExpr, x)

Base.isempty(a::AffExpr) = (length(a.vars) == 0 && a.constant == 0.)
Base.convert(::Type{AffExpr}, v::Variable) = AffExpr([v], [1.], 0.)
Base.convert(::Type{AffExpr}, v::Real) = AffExpr(Variable[], Float64[], v)

function assert_isfinite(a::AffExpr)
    coeffs = a.coeffs
    for i in 1:length(a.vars)
        isfinite(coeffs[i]) || error("Invalid coefficient $(coeffs[i]) on variable $(a.vars[i])")
    end
end

setObjective(m::Model, sense::Symbol, x::Variable) = setObjective(m, sense, convert(AffExpr,x))
function setObjective(m::Model, sense::Symbol, a::AffExpr)
    setObjectiveSense(m, sense)
    m.obj = zero(QuadExpr)
    m.obj.aff = a
end

Base.copy(a::AffExpr, new_model::Model) =
    AffExpr([Variable(new_model, v.col) for v in a.vars], a.coeffs[:], a.constant)

function getValue(a::AffExpr)
    ret = a.constant
    for it in 1:length(a.vars)
        ret += a.coeffs[it] * getValue(a.vars[it])
    end
    return ret
end
getValue(arr::Array{AffExpr}) = map(getValue, arr)

###############################################################################
# QuadExpr class
# Holds a vector of tuples (Var, Var, Coeff), as well as an AffExpr
type GenericQuadExpr{CoefType,VarType}
    qvars1::Vector{VarType}
    qvars2::Vector{VarType}
    qcoeffs::Vector{CoefType}
    aff::GenericAffExpr{CoefType,VarType}
end

coeftype{CoefType,VarType}(::GenericQuadExpr{CoefType,VarType}) = CoefType

typealias QuadExpr GenericQuadExpr{Float64,Variable}

Base.isempty(q::QuadExpr) = (length(q.qvars1) == 0 && isempty(q.aff))
Base.zero{C,V}(::Type{GenericQuadExpr{C,V}}) = GenericQuadExpr(V[], V[], C[], zero(GenericAffExpr{C,V}))
Base.one{C,V}(::Type{GenericQuadExpr{C,V}})  = GenericQuadExpr(V[], V[], C[],  one(GenericAffExpr{C,V}))
Base.zero(q::GenericQuadExpr) = zero(typeof(q))
Base.one(q::GenericQuadExpr)  =  one(typeof(q))

QuadExpr() = zero(QuadExpr)

Base.convert(::Type{QuadExpr}, v::Union(Real,Variable,AffExpr)) = QuadExpr(Variable[], Variable[], Float64[], AffExpr(v))

function Base.append!{T,S}(q::GenericQuadExpr{T,S}, other::GenericQuadExpr{T,S})
    append!(q.qvars1, other.qvars1)
    append!(q.qvars2, other.qvars2)
    append!(q.qcoeffs, other.qcoeffs)
    append!(q.aff, other.aff)
    q
end

function assert_isfinite(q::GenericQuadExpr)
    assert_isfinite(q.aff)
    for i in 1:length(q.qcoeffs)
        isfinite(q.qcoeffs[i]) || error("Invalid coefficient $(q.qcoeffs[i]) on quadratic term $(q.qvars1[i])*$(q.qvars2[i])")
    end
end

function setObjective(m::Model, sense::Symbol, q::QuadExpr)
    m.obj = q
    if m.internalModelLoaded
        if method_exists(MathProgBase.setquadobjterms!, (typeof(m.internalModel), Vector{Cint}, Vector{Cint}, Vector{Float64}))
            verify_ownership(m, m.obj.qvars1)
            verify_ownership(m, m.obj.qvars2)
            MathProgBase.setquadobjterms!(m.internalModel, Cint[v.col for v in m.obj.qvars1], Cint[v.col for v in m.obj.qvars2], m.obj.qcoeffs)
        else
            # we don't (yet) support hot-starting QCQP solutions
            Base.warn_once("JuMP does not yet support adding a quadratic objective to an existing model. Hot-start is disabled.")
            m.internalModelLoaded = false
        end
    end
    setObjectiveSense(m, sense)
end

Base.copy(q::GenericQuadExpr) = GenericQuadExpr(copy(q.qvars1),copy(q.qvars2),copy(q.qcoeffs),copy(q.aff))

# Copy utility function
function Base.copy(q::QuadExpr, new_model::Model)
    return QuadExpr([Variable(new_model, v.col) for v in q.qvars1],
                    [Variable(new_model, v.col) for v in q.qvars2],
                    q.qcoeffs[:], copy(q.aff, new_model))
end

Base.zero(::Type{QuadExpr}) = QuadExpr(Variable[],Variable[],Float64[],zero(AffExpr))
Base.zero(v::QuadExpr) = zero(typeof(v))

function getValue(a::QuadExpr)
    ret = getValue(a.aff)
    for it in 1:length(a.qvars1)
        ret += a.qcoeffs[it] * getValue(a.qvars1[it]) * getValue(a.qvars2[it])
    end
    return ret
end

getValue(arr::Array{QuadExpr}) = map(getValue, arr)

##########################################################################
# GenericNorm, Norm
# Container for √(∑ expr)
type GenericNorm{P,C,V}
    terms::Vector{GenericAffExpr{C,V}}
end
function GenericNorm{C,V}(P, terms::Vector{GenericAffExpr{C,V}})
    if C == Float64 && V == Variable
        # JuMP doesn't support anything else than the L2 Norm
        # Throw an error now before going any further
        P ==   1 && error("JuMP doesn't support L₁ norms.")
        P == Inf && error("JuMP doesn't support L∞ norms.")
        P !=   2 && error("JuMP only supports L₂ norms.")
    end
    GenericNorm{P,C,V}(terms)
end
Base.copy{P,C,V}(x::GenericNorm{P,C,V}) = GenericNorm{P,C,V}(copy(x.terms))

Base.norm{V<:AbstractJuMPScalar}(x::V,           p=2) = vecnorm(x,p)
Base.norm{V<:AbstractJuMPScalar}(x::Array{V},    p=2) = vecnorm(x,p)
Base.norm{V<:AbstractJuMPScalar}(x::JuMPArray{V},p=2) = vecnorm(x,p)
Base.norm{V<:AbstractJuMPScalar}(x::JuMPDict{V}, p=2) = vecnorm(x,p)
Base.norm{C,V}(x::GenericAffExpr{C,V},           p=2) = vecnorm(x,p)
Base.norm{C,V}(x::Array{GenericAffExpr{C,V}},    p=2) = vecnorm(x,p)
Base.norm{C,V}(x::JuMPArray{GenericAffExpr{C,V}},p=2) = vecnorm(x,p)
Base.norm{C,V}(x::JuMPDict{GenericAffExpr{C,V}}, p=2) = vecnorm(x,p)

_vecaff(C,V,x) = map(GenericAffExpr{C,V},vec(x))
Base.vecnorm{V<:AbstractJuMPScalar}(x::V,           p=2) = GenericNorm(p, [GenericAffExpr{Float64,V}(x)] )
Base.vecnorm{V<:AbstractJuMPScalar}(x::Array{V},    p=2) = GenericNorm(p, _vecaff(Float64,V,x) )
Base.vecnorm{V<:AbstractJuMPScalar}(x::JuMPArray{V},p=2) = GenericNorm(p, _vecaff(Float64,V,x.innerArray) )
Base.vecnorm{V<:AbstractJuMPScalar}(x::JuMPDict{V}, p=2) = GenericNorm(p, _vecaff(Float64,V,collect(values(x))) )
Base.vecnorm{C,V}(x::GenericAffExpr{C,V},           p=2) = GenericNorm(p, [x])
Base.vecnorm{C,V}(x::Array{GenericAffExpr{C,V}},    p=2) = GenericNorm(p, vec(x))
Base.vecnorm{C,V}(x::JuMPArray{GenericAffExpr{C,V}},p=2) = GenericNorm(p, vec(x.innerArray))
Base.vecnorm{C,V}(x::JuMPDict{GenericAffExpr{C,V}}, p=2) = GenericNorm(p, collect(values(x)))

# Called by macros
_build_norm{C,V}(P,terms::Vector{GenericAffExpr{C,V}}) = GenericNorm(P,terms)
# The terms vector produced by the macro may not be tightly typed,
# so we need to repackage in a tighter typed vector ourselves.
# This function is needed for performance reasons on only 0.3, as
# the following works just as well on 0.4:
# _build_norm(Lp, terms::Vector{GenericAffExpr}) = _build_norm(Lp, [terms...])
function _build_norm(Lp, terms::Vector{GenericAffExpr})
    if length(terms) == 0
        _build_norm(Lp,AffExpr[])  # Punt
    else
        new_terms = Array(typeof(terms[1]), length(terms))
        for i in 1:length(terms)
            new_terms[i] = terms[i]
        end
        _build_norm(Lp,new_terms)
    end
end

# Short-hand used in operator overloads, etc.
typealias Norm{P} GenericNorm{P,Float64,Variable}

##########################################################################
# GenericNormExpr, SOCExpr
# Container for expressions containing Norms and AffExprs
type GenericNormExpr{P,C,V}
    norm::GenericNorm{P,C,V}
    coeff::C
    aff::GenericAffExpr{C,V}
end

GenericNormExpr{P,C,V}(norm::GenericNorm{P,C,V}) =
    GenericNormExpr{C,V}(norm, one(C), zero(GenericAffExpr{C,V}))
Base.copy{P,C,V}(x::GenericNormExpr{P,C,V}) =
    GenericNormExpr{P,C,V}(copy(x.norm), copy(x.coeff), copy(x.aff))
Base.convert{P,C,V}(::Type{GenericNormExpr{P,C,V}}, x::GenericNorm{P,C,V}) =
    GenericNormExpr{P,C,V}(x, one(C), zero(GenericAffExpr{C,V}))


typealias GenericSOCExpr{C,V} GenericNormExpr{2,C,V}
typealias SOCExpr GenericSOCExpr{Float64,Variable}

validate_soc(socexpr::GenericSOCExpr) = (socexpr.coeff ≥ 0) ||
    error("Invalid second-order cone constraint $(socexpr) ≥ 0")

##########################################################################
# JuMPConstraint
# abstract base for constraint types
abstract JuMPConstraint

##########################################################################
# Generic constraint type with lower and upper bound
type GenericRangeConstraint{TermsType} <: JuMPConstraint
    terms::TermsType
    lb::Float64
    ub::Float64
end

function sense(c::GenericRangeConstraint)
    if c.lb != -Inf
        if c.ub != Inf
            if c.ub == c.lb
                return :(==)
            else
                return :range
            end
        else
                return :>=
        end
    else
        @assert c.ub != Inf
        return :<=
    end
end

function rhs(c::GenericRangeConstraint)
    s = sense(c)
    @assert s != :range
    if s == :<=
        return c.ub
    else
        return c.lb
    end
end

##########################################################################
# LinearConstraint is an affine expression with lower bound (possibly
# -Inf) and upper bound (possibly Inf).
typealias LinearConstraint GenericRangeConstraint{AffExpr}

function addConstraint(m::Model, c::LinearConstraint)
    push!(m.linconstr,c)
    if m.internalModelLoaded
        if method_exists(MathProgBase.addconstr!, (typeof(m.internalModel),Vector{Int},Vector{Float64},Float64,Float64))
            assert_isfinite(c.terms)
            indices, coeffs = merge_duplicates(Cint, c.terms, m.indexedVector, m)
            MathProgBase.addconstr!(m.internalModel,indices,coeffs,c.lb,c.ub)
        else
            Base.warn_once("Solver does not appear to support adding constraints to an existing model. Hot-start is disabled.")
            m.internalModelLoaded = false
        end
    end
    return ConstraintRef{LinearConstraint}(m,length(m.linconstr))
end
addConstraint(m::Model, c::Array{LinearConstraint}) =
    error("The operators <=, >=, and == can only be used to specify scalar constraints. If you are trying to add a vectorized constraint, use the element-wise dot comparison operators (.<=, .>=, or .==) instead")

addVectorizedConstraint(m::Model, v::Array{LinearConstraint}) = map(c->addConstraint(m,c), v)

function Base.copy(c::LinearConstraint, new_model::Model)
    return LinearConstraint(copy(c.terms, new_model), c.lb, c.ub)
end

##########################################################################
# SOSConstraint class
# An SOS constraint.
type SOSConstraint <: JuMPConstraint
    terms::Vector{Variable}
    weights::Vector{Float64}
    sostype::Symbol
end

function constructSOS(m::Model, coll::Vector{AffExpr})
    nvar = length(coll)
    vars = Array(Variable, nvar)
    weight = Array(Float64, nvar)
    for i in 1:length(coll)
        if (length(coll[i].vars) != 1) || (coll[i].constant != 0)
            error("Must specify collection in terms of single variables")
        end
        vars[i] = coll[i].vars[1]
        vars[i].m == m || error("Variable in constraint is not owned by the model")
        weight[i] = coll[i].coeffs[1]
    end
    return vars, weight
end

addSOS1(m::Model, coll) = addSOS1(m, convert(Vector{AffExpr}, coll))

function addSOS1(m::Model, coll::Vector{AffExpr})
    vars, weight = constructSOS(m,coll)
    push!(m.sosconstr, SOSConstraint(vars, weight, :SOS1))
    if m.internalModelLoaded
        idx = Int[v.col for v in vars]
        if applicable(MathProgBase.addsos1!, m.internalModel, idx, weight)
            MathProgBase.addsos1!(m.internalModel, idx, weight)
        else
            Base.warn_once("Solver does not appear to support adding constraints to an existing model. Hot-start is disabled.")
            m.internalModelLoaded = false
        end
    end
    return ConstraintRef{SOSConstraint}(m,length(m.sosconstr))
end

addSOS2(m::Model, coll) = addSOS2(m, convert(Vector{AffExpr}, coll))

function addSOS2(m::Model, coll::Vector{AffExpr})
    vars, weight = constructSOS(m,coll)
    push!(m.sosconstr, SOSConstraint(vars, weight, :SOS2))
    if m.internalModelLoaded
        idx = Int[v.col for v in vars]
        if applicable(MathProgBase.addsos2!, m.internalModel, idx, weight)
            MathProgBase.addsos2!(m.internalModel, idx, weight)
        else
            Base.warn_once("Solver does not appear to support adding constraints to an existing model. Hot-start is disabled.")
            m.internalModelLoaded = false
        end
    end
    return ConstraintRef{SOSConstraint}(m,length(m.sosconstr))
end

Base.copy(sos::SOSConstraint, new_model::Model) =
    SOSConstraint([Variable(new_model,v.col) for v in sos.terms], copy(sos.weights), sos.sostype)


##########################################################################
# SDPConstraint is a (dual) semidefinite constraint of the form
# ∑ cᵢ Xᵢ ≥ D, where D is a n×n symmetric data matrix, cᵢ are
# scalars, and Xᵢ are n×n symmetric variable matrices. The inequality
# is taken w.r.t. the psd partial order.
type SDPConstraint <: JuMPConstraint
    terms # purposely leave this untyped so that we can special-case OneIndexedArray with no additional variables
end

# Special-case X ≥ 0, which is often convenient
function SDPConstraint(lhs::Union(OneIndexedArray,Matrix), rhs::Number)
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
# Generic constraint type for quadratic expressions
# Right-hand side is implicitly taken to be zero, constraint is stored in
# the included QuadExpr.
type GenericQuadConstraint{QuadType} <: JuMPConstraint
    terms::QuadType
    sense::Symbol
end

##########################################################################
# QuadConstraint class
typealias QuadConstraint GenericQuadConstraint{QuadExpr}

function addConstraint(m::Model, c::QuadConstraint)
    push!(m.quadconstr,c)
    if m.internalModelLoaded
        if method_exists(MathProgBase.addquadconstr!, (typeof(m.internalModel),
                                                       Vector{Cint},
                                                       Vector{Float64},
                                                       Vector{Cint},
                                                       Vector{Cint},
                                                       Vector{Float64},
                                                       Char,
                                                       Float64))
            if !((s = string(c.sense)[1]) in ['<', '>', '='])
                error("Invalid sense for quadratic constraint")
            end
            terms = c.terms
            verify_ownership(m, terms.qvars1)
            verify_ownership(m, terms.qvars2)
            MathProgBase.addquadconstr!(m.internalModel,
                                        Cint[v.col for v in c.terms.aff.vars],
                                        c.terms.aff.coeffs,
                                        Cint[v.col for v in c.terms.qvars1],
                                        Cint[v.col for v in c.terms.qvars2],
                                        c.terms.qcoeffs,
                                        s,
                                        -c.terms.aff.constant)
        else
            Base.warn_once("Solver does not appear to support adding quadratic constraints to an existing model. Hot-start is disabled.")
            m.internalModelLoaded = false
        end
    end
    return ConstraintRef{QuadConstraint}(m,length(m.quadconstr))
end
addConstraint(m::Model, c::Array{QuadConstraint}) =
    error("Vectorized constraint added without elementwise comparisons. Try using one of (.<=,.>=,.==).")

addVectorizedConstraint(m::Model, v::Array{QuadConstraint}) = map(c->addConstraint(m,c), v)

# Copy utility function
function Base.copy(c::QuadConstraint, new_model::Model)
    return QuadConstraint(copy(c.terms, new_model), c.sense)
end

##########################################################################
# SOCConstraint is a second-order cone constraint of the form
# α||Ax-b||₂ + cᵀx + d ≤ 0

type GenericSOCConstraint{T<:GenericSOCExpr} <: JuMPConstraint
    normexpr::T

    function GenericSOCConstraint{T}(normexpr::T)
        validate_soc(normexpr)
        new(normexpr)
    end
end

typealias SOCConstraint GenericSOCConstraint{SOCExpr}

function addConstraint{T<:GenericSOCConstraint}(m::Model, c::T)
    push!(m.socconstr,c)
    m.internalModelLoaded = false
    ConstraintRef{T}(m,length(m.socconstr))
end

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
    constraints::JuMPArray,coefficients::Vector{Float64}, name::String="", value::Number=NaN) =
    Variable(m, lower, upper, cat, objcoef, constraints.innerArray, coefficients, name, value)

# add variable to existing constraints
function Variable(m::Model,lower::Number,upper::Number,cat::Symbol,objcoef::Number,
    constraints::Vector,coefficients::Vector{Float64}, name::String="", value::Number=NaN)
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
            if m.operator_counter > 200
                Base.warn_once("The addition operator has been used on JuMP expressions a large number of times. For performance reasons, you should not add expressions in a loop. Instead of x += y, use append!(x,y) to modify x in place. If y is a single variable, you may also use push!(x, coef, y) in place of x += coef*y.")
            end
        end
    end
    return
end
operator_warn(lhs,rhs) = nothing

##########################################################################
# Behavior that's uniform across all JuMP "scalar" objects

typealias JuMPTypes Union(AbstractJuMPScalar,
                          Norm,
                          GenericAffExpr,
                          QuadExpr,
                          SOCExpr)
typealias JuMPScalars Union(Number,JuMPTypes)

# would really want to do this on ::Type{T}, but doesn't work on v0.4
Base.eltype{T<:JuMPTypes}(::T) = T
Base.size(::JuMPTypes) = ()
Base.size(x::JuMPTypes,d::Int) = 1
Base.ndims(::JuMPTypes) = 0

##########################################################################
# Operator overloads

include("operators.jl")
if VERSION > v"0.4-"
    include(joinpath("v0.4","concatenation.jl"))
else
    include(joinpath("v0.3","concatenation.jl"))
end
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

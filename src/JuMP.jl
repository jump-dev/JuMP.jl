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

importall Base.Operators
import Base.map

import MathProgBase
import MathOptInterface
import MathOptInterfaceUtilities
const MOI = MathOptInterface
const MOIU = MathOptInterfaceUtilities

using Calculus
using ReverseDiffSparse
using ForwardDiff

export
# Objects
    Model, Variable, Norm, AffExpr, QuadExpr, SOCExpr,
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
    solve,
    internalmodel,
    # Variable
    setname,
    #getname,
    setlowerbound, setupperbound,
    #getlowerbound, getupperbound,
    #getvalue, setvalue,
    #getdual,
    #setcategory, getcategory,
    setstart,
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
const LBREF = MOICON{MOI.SingleVariable,MOI.GreaterThan{Float64}}
const UBREF = MOICON{MOI.SingleVariable,MOI.LessThan{Float64}}
const FIXREF = MOICON{MOI.SingleVariable,MOI.EqualTo{Float64}}
const INTREF = MOICON{MOI.SingleVariable,MOI.Integer}
const BINREF = MOICON{MOI.SingleVariable,MOI.ZeroOne}

@MOIU.instance JuMPInstance (ZeroOne, Integer) (EqualTo, GreaterThan, LessThan, Interval) (Zeros, Nonnegatives, Nonpositives, SecondOrderCone, RotatedSecondOrderCone, GeometricMeanCone, PositiveSemidefiniteConeTriangle, PositiveSemidefiniteConeSquare, RootDetConeTriangle, RootDetConeSquare, LogDetConeTriangle, LogDetConeSquare) () (SingleVariable,) (ScalarAffineFunction,ScalarQuadraticFunction) (VectorOfVariables,) (VectorAffineFunction,)

###############################################################################
# Model class
# Keeps track of all model and column info
abstract type AbstractModel end
mutable struct Model <: AbstractModel

    instance::JuMPInstance{Float64}
    # special variablewise properties that we keep track of:
    # lower bound, upper bound, fixed, integrality, binary
    variabletolowerbound::Dict{MOIVAR,LBREF}
    variabletoupperbound::Dict{MOIVAR,UBREF}
    variabletofix::Dict{MOIVAR,FIXREF}
    variabletointegrality::Dict{MOIVAR,INTREF}
    variabletozeroone::Dict{MOIVAR,BINREF}
    variabletosolvervariable::Dict{MOIVAR,MOIVAR}

    # convenient solution objects to keep in the model
    variablestart #VariableToValueMap{Float64}
    variableresult #VariableToValueMap{Float64}


    # obj#::QuadExpr
    # objSense::Symbol

    # Mapping from the constraint reference in `instance`
    # and the constraint reference in `solverinstance`
    constrainttosolverconstraint::Dict{MOICON,MOICON}

    # linconstr#::Vector{LinearConstraint}
    # quadconstr
    # sosconstr
    # socconstr
    # sdpconstr

    # Column data
    # numCols::Int
    # colNames::Vector{String}
    # colNamesIJulia::Vector{String}
    # colLower::Vector{Float64}
    # colUpper::Vector{Float64}
    # colCat::Vector{Symbol}

    variablenames # VariableToValueMap{String}

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
    # internal solver instance object
    solverinstance
    # Solver+option object from MPB or MOI
    solverinstanceattached::Bool
    # callbacks
    callbacks
    # lazycallback
    # cutcallback
    # heurcallback

    # hook into a solve call...function of the form f(m::Model; kwargs...),
    # where kwargs get passed along to subsequent solve calls
    solvehook
    # # ditto for a print hook
    # printhook


    # # storage vector for merging duplicate terms
    # indexedVector::IndexedVector{Float64}

    nlpdata#::NLPData
    simplify_nonlinear_expressions::Bool

    objdict::Dict{Symbol,Any} # dictionary from variable and constraint names to objects

    operator_counter::Int # number of times we add large expressions

    # Extension dictionary - e.g. for robust
    # Extensions should define a type to hold information particular to
    # their functionality, and store an instance of the type in this
    # dictionary keyed on an extension-specific symbol
    ext::Dict{Symbol,Any}
    # Default constructor
    function Model(; simplify_nonlinear_expressions::Bool=false)
        # TODO need to support MPB also
        m = new()
        # TODO make pretty
        m.instance = JuMPInstance{Float64}()
        m.variabletolowerbound = Dict{MOIVAR,LBREF}()
        m.variabletoupperbound = Dict{MOIVAR,UBREF}()
        m.variabletofix = Dict{MOIVAR,FIXREF}()
        m.variabletointegrality = Dict{MOIVAR,INTREF}()
        m.variabletozeroone = Dict{MOIVAR,BINREF}()
        m.variablestart = VariableToValueMap{Float64}(m)
        m.variableresult = VariableToValueMap{Float64}(m)
        m.variablenames = VariableToValueMap{String}(m)
        m.customnames = Variable[]
        m.objbound = 0.0
        m.objval = 0.0
        m.solverinstance = nothing
        m.solverinstanceattached = false
        m.callbacks = Any[]
        m.solvehook = nothing
        # m.printhook = nothing
        m.nlpdata = nothing
        m.simplify_nonlinear_expressions = simplify_nonlinear_expressions
        m.objdict = Dict{Symbol,Any}()
        m.operator_counter = 0
        m.ext = Dict{Symbol,Any}()

        return m
    end
end





# Getters/setters

# temporary name
numvar(m::Model) = MOI.get(m.instance, MOI.NumberOfVariables())

# """
#     MathProgBase.numvar(m::Model)
#
# returns the number of variables associated with the `Model m`.
# """
# MathProgBase.numvar(m::Model) = m.numCols
#
# """
#     MathProgBase.numlinconstr(m::Model)
#
# returns the number of linear constraints associated with the `Model m`
# """
# MathProgBase.numlinconstr(m::Model) = length(m.linconstr)
#
# """
#     MathProgBase.numquadconstr(m::Model)
#
# returns the number of quadratic constraints associated with the `Model m`
# """
# MathProgBase.numquadconstr(m::Model) = length(m.quadconstr)

# """
#     numsocconstr(m::Model)
#
# returns the number of second order cone constraints associated with the `Model m`
# """
# numsocconstr(m::Model) = length(m.socconstr)
#
# """
#     numsosconstr(m::Model)
#
# returns the number of sos constraints associated with the `Model m`
# """
# numsosconstr(m::Model) = length(m.sosconstr)

# """
#     numsdconstr(m::Model)
#
# returns the number of semi-definite constraints associated with the `Model m`
# """
# numsdconstr(m::Model) = length(m.sdpconstr)

"""
    numnlconstr(m::Model)

returns the number of nonlinear constraints associated with the `Model m`
"""
numnlconstr(m::Model) = m.nlpdata !== nothing ? length(m.nlpdata.nlconstr) : 0

# """
#     MathProgBase.numconstr(m::Model)
#
# returns the total number of constraints associated with the `Model m`
# """
# function MathProgBase.numconstr(m::Model)
#     c = length(m.linconstr) + length(m.quadconstr) + length(m.socconstr) + length(m.sosconstr) + length(m.sdpconstr)
#     if m.nlpdata !== nothing
#         c += length(m.nlpdata.nlconstr)
#     end
#     return c
# end
#
# for f in MathProgBase.SolverInterface.methods_by_tag[:rewrap]
#     eval(Expr(:import,:MathProgBase,f))
#     @eval function $f(m::Model)
#         # check internal model exists
#         if !m.internalModelLoaded
#             error("Model not solved")
#         else
#             return $f(internalmodel(m))
#         end
#     end
#     eval(Expr(:export,f))
# end


# Doc strings for the auto-wrapped MPB functions above
# it would be preferable to problematically use the docstrings from MPB functions instead

# @doc """
#     getsolvetime(m::Model)
#
# returns the solve time reported by the solver if it is implemented.
# """ getsolvetime(m::Model)
#
# @doc """
#     getnodecount(m::Model)
#
# returns the number of explored branch-and-bound nodes, if it is implemented.
# """ getnodecount(m::Model)
#
# @doc """
#     getobjbound(m::Model)
#
# returns the best known bound on the optimal objective value. This is used, for example, when a branch-and-bound method is stopped before finishing.
# """ getobjbound(m::Model)
#
# @doc """
#     getobjgap(m::Model)
#
# returns the final relative optimality gap as optimization terminated. That is, it returns ``\\frac{|b-f|}{|f|}``, where *b* is the best bound and *f* is the best feasible objective value.
# """ getobjgap(m::Model)
#
# @doc """
#     getrawsolver(m::Model)
#
# returns an object that may be used to access a solver-specific API.
# """ getrawsolver(m::Model)
#
# @doc """
#     getsimplexiter(m::Model)
#
# returns the cumulative number of simplex iterations during the optimization process. In particular, for a MIP it returns the total simplex iterations for all nodes.
# """ getsimplexiter(m::Model)
#
# @doc """
#     getbarrieriter(m::Model)
#
# returns the cumulative number of barrier iterations during the optimization process.
# """ getbarrieriter(m::Model)


# """
#     getobjective(m::Model)
#
# returns the objective function as a `QuadExpr`
# """
# function getobjective(m::Model)
#     traits = ProblemTraits(m)
#     if traits.nlp
#         error("getobjective() not supported for nonlinear models")
#     end
#     return m.obj
# end


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
    moisense = MOI.get(m.instance, MOI.ObjectiveSense())
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

# """
#     setobjectivesense(m::Model, newSense::Symbol)
#
# sets the objective sense (`newSense` is either `:Min` or `:Max`)
# """
# function setobjectivesense(m::Model, newSense::Symbol)
#     if (newSense != :Max && newSense != :Min)
#         error("Model sense must be :Max or :Min")
#     end
#     m.objSense = newSense
# end
# setobjective(m::Model, something::Any) =
#     error("in setobjective: needs three arguments: model, objective sense (:Max or :Min), and expression.")
#
# setobjective(::Model, ::Symbol, x::AbstractArray) =
#     error("in setobjective: array of size $(_size(x)) passed as objective; only scalar objectives are allowed")

# # Deep copy the model
# function Base.copy(source::Model)
#
#     dest = Model()
#     dest.solver = source.solver  # The two models are linked by this
#
#     # Objective
#     dest.obj = copy(source.obj, dest)
#     dest.objSense = source.objSense
#
#     # Constraints
#     dest.linconstr  = map(c->copy(c, dest), source.linconstr)
#     dest.quadconstr = map(c->copy(c, dest), source.quadconstr)
#     dest.sosconstr  = map(c->copy(c, dest), source.sosconstr)
#     dest.sdpconstr  = map(c->copy(c, dest), source.sdpconstr)
#
#     # Variables
#     dest.numCols = source.numCols
#     dest.colNames = source.colNames[:]
#     dest.colNamesIJulia = source.colNamesIJulia[:]
#     dest.colLower = source.colLower[:]
#     dest.colUpper = source.colUpper[:]
#     dest.colCat = source.colCat[:]
#
#     # varCones
#     dest.varCones = copy(source.varCones)
#
#     # callbacks and hooks
#     if !isempty(source.callbacks)
#         error("Copying callbacks is not supported")
#     end
#     if source.solvehook !== nothing
#         dest.solvehook = source.solvehook
#     end
#     if source.printhook !== nothing
#         dest.printhook = source.printhook
#     end
#
#     # variable/extension dicts
#     if !isempty(source.ext)
#         dest.ext = similar(source.ext)
#         for (key, val) in source.ext
#             dest.ext[key] = try
#                 copy(source.ext[key])
#             catch
#                 error("Error copying extension dictionary. Is `copy` defined for all your user types?")
#             end
#         end
#     end
#     dest.objDict = Dict{Symbol,Any}()
#     dest.varData = ObjectIdDict()
#     for (symb,o) in source.objDict
#         newo = copy(o, dest)
#         dest.objDict[symb] = newo
#         if haskey(source.varData, o)
#             dest.varData[newo] = source.varData[o]
#             #dest.varData[newvar] = copy(source.varData[v]) # should we copy this too ? We need to define copy(::JuMPContainerData) too then
#         end
#     end
#
#     if source.nlpdata !== nothing
#         dest.nlpdata = copy(source.nlpdata)
#     end
#
#     return dest
# end

"""
    solverinstance(m::Model)

returns the internal `AbstractSolverInstance` object which can be used to access any functionality that is not exposed by JuMP.
See the MathOptInterface [documentation](XXX).
"""
solverinstance(m::Model) = m.solverinstance

setsolvehook(m::Model, f) = (m.solvehook = f)
setprinthook(m::Model, f) = (m.printhook = f)


#############################################################################
# AbstractConstraint
# Abstract base type for all constraint types
abstract type AbstractConstraint end
# Abstract base type for all scalar types
# In JuMP, used only for Variable. Useful primarily for extensions
abstract type AbstractJuMPScalar end

Base.start(::AbstractJuMPScalar) = false
Base.next(x::AbstractJuMPScalar, state) = (x, true)
Base.done(::AbstractJuMPScalar, state) = state
Base.isempty(::AbstractJuMPScalar) = false

include("variables.jl")

Base.zero(::Type{Variable}) = AffExpr(Variable[],Float64[],0.0)
Base.zero(::Variable) = zero(Variable)
Base.one(::Type{Variable}) = AffExpr(Variable[],Float64[],1.0)
Base.one(::Variable) = one(Variable)

mutable struct VariableNotOwnedError <: Exception
    context::String
end
function Base.showerror(io::IO, ex::VariableNotOwnedError)
    print(io, "VariableNotOwnedError: Variable not owned by model present in $(ex.context)")
end

function verify_ownership(m::Model, vec::Vector{Variable})
    n = length(vec)
    @inbounds for i in 1:n
        vec[i].m !== m && return false
    end
    return true
end

Base.copy(v::Variable, new_model::Model) = Variable(new_model, v.col)
Base.copy(x::Void, new_model::Model) = nothing
Base.copy(v::AbstractArray{Variable}, new_model::Model) = (var -> Variable(new_model, var.col)).(v)

##########################################################################
# ConstraintRef
# Reference to a constraint for retrieving solution info
struct ConstraintRef{M<:AbstractModel,C}
    m::M
    instanceref::C
end
# Base.copy{M,T}(c::ConstraintRef{M,T}, new_model::M) = ConstraintRef{M,T}(new_model, c.idx)

# linearindex(x::ConstraintRef) = x.idx

function solverinstanceref(cr::ConstraintRef{Model, MOICON{F, S}}) where {F, S}
    cr.m.constrainttosolverconstraint[cr.instanceref]::MOICON{F, S}
end

function hasresultdual(cr::ConstraintRef{Model, <:MOICON})
    MOI.get(cr.m.solverinstance, MOI.ConstraintDual(), cr.instanceref)
end

"""
    resultdual(cr::ConstraintRef)

Get the dual value of this constraint in the result returned by a solver.
Use `hasresultdual` to check if a result exists before asking for values.
Replaces `getdual` for most use cases.
"""
function resultdual(cr::ConstraintRef{Model, MOICON{F, S}}) where {F, S}
    MOI.get(cr.m.solverinstance, MOI.ConstraintDual(), solverinstanceref(cr))
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

# # internal method that doesn't print a warning if the value is NaN
# _getDual(c::LinConstrRef) = c.m.linconstrDuals[c.idx]
#
# getdualwarn{T<:Union{ConstraintRef, Int}}(::T) = warn("Dual solution not available. Check that the model was properly solved and no integer variables are present.")
#
# """
#     getdual(c::LinConstrRef)
#
# """
# function getdual(c::LinConstrRef)
#     if length(c.m.linconstrDuals) != MathProgBase.numlinconstr(c.m)
#         getdualwarn(c)
#         NaN
#     else
#         _getDual(c)
#     end
# end

# Returns the number of non-infinity and nonzero bounds on variables
# function getNumBndRows(m::Model)
#     numBounds = 0
#     for i in 1:m.numCols
#         seen = false
#         lb, ub = m.colLower[i], m.colUpper[i]
#         for (_,cone) in m.varCones
#             if i in cone
#                 seen = true
#                 @assert lb == -Inf && ub == Inf
#                 break
#             end
#         end
#
#         if !seen
#             if lb != -Inf && lb != 0
#                 numBounds += 1
#             end
#             if ub != Inf && ub != 0
#                 numBounds += 1
#             end
#         end
#     end
#     return numBounds
# end

# Returns the number of second-order cone constraints
# getNumRows(c::SOCConstraint) = length(c.normexpr.norm.terms) + 1
# getNumSOCRows(m::Model) = sum(getNumRows.(m.socconstr))

# Returns the dual variables corresponding to
# m.sdpconstr[idx] if issdp is true
# m.socconstr[idx] if sdp is not true
# function getconicdualaux(m::Model, idx::Int, issdp::Bool)
#     numLinRows = MathProgBase.numlinconstr(m)
#     numBndRows = getNumBndRows(m)
#     numSOCRows = getNumSOCRows(m)
#     numSDPRows = getNumSDPRows(m)
#     numSymRows = getNumSymRows(m)
#     numRows = numLinRows + numBndRows + numSOCRows + numSDPRows + numSymRows
#     if length(m.conicconstrDuals) != numRows
#         # solve might not have been called so m.constr_to_row might be empty
#         getdualwarn(idx)
#         c = issdp ? m.sdpconstr[idx] : m.socconstr[idx]
#         duals = fill(NaN, getNumRows(c))
#         if issdp
#             duals, Float64[]
#         else
#             duals
#         end
#     else
#         offset = numLinRows + numBndRows
#         if issdp
#             offset += length(m.socconstr)
#         end
#         dual = m.conicconstrDuals[m.constr_to_row[offset + idx]]
#         if issdp
#             offset += length(m.sdpconstr)
#             symdual = m.conicconstrDuals[m.constr_to_row[offset + idx]]
#             dual, symdual
#         else
#             dual
#         end
#     end
# end

# """
#     getdual(c::ConstraintRef{Model,SOCConstraint})
#
#
# """
# function getdual(c::ConstraintRef{Model,SOCConstraint})
#     getconicdualaux(c.m, c.idx, false)
# end

# function setRHS(c::LinConstrRef, rhs::Number)
#     constr = c.m.linconstr[c.idx]
#     sen = sense(constr)
#     if sen == :range
#         error("Modifying range constraints is currently unsupported.")
#     elseif sen == :(==)
#         constr.lb = float(rhs)
#         constr.ub = float(rhs)
#     elseif sen == :>=
#         constr.lb = float(rhs)
#     else
#         @assert sen == :<=
#         constr.ub = float(rhs)
#     end
# end

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

To allow easy accessing of JuMP Variables and Constraints via `[]` syntax.
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
                        GenericAffExpr,
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
include("solverinterface.jl")
# include("callbacks.jl")
include("nlp.jl")
include("print.jl")

# getvalue{T<:JuMPTypes}(arr::Array{T}) = map(getvalue, arr)
#
# function setvalue{T<:AbstractJuMPScalar}(set::Array{T}, val::Array)
#     promote_shape(size(set), size(val)) # Check dimensions match
#     for I in eachindex(set)
#         setvalue(set[I], val[I])
#     end
#     nothing
# end


##########################################################################
end

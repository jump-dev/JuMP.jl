#############################################################################
# JuMP
# An algebraic modelling langauge for Julia
# See http://github.com/JuliaOpt/JuMP.jl
#############################################################################

import Base.getindex
import Base.setindex!
import Base.print
import Base.show

module JuMP

# Use the standard solver interface for LPs and MIPs
using MathProgBase
require(joinpath(Pkg.dir("MathProgBase"),"src","MathProgSolverInterface.jl"))
using MathProgSolverInterface

importall Base

export
# Objects
  Model, Variable, AffExpr, QuadExpr, LinearConstraint, QuadConstraint, MultivarDict,
# Reexported from MathProgBase
  LPSolver, MIPSolver,
# Functions
  # Relevant to all
  print,show,
  # Model related
  getNumVars, getNumConstraints, getObjectiveValue, getObjective,
  getObjectiveSense, setObjectiveSense, writeLP, writeMPS, setObjective,
  addConstraint, addVar, addVars, solve,
  # Variable
  setName, getName, setLower, setUpper, getLower, getUpper, getValue,
  getDual,
  # Expressions and constraints
  affToStr, quadToStr, conToStr,
  
# Macros and support functions
  @addConstraint, @defVar, 
  @defConstrRef, @setObjective, addToExpression

include("JuMPDict.jl")
include("utils.jl")

###############################################################################
# Constants
const CONTINUOUS = 0
const INTEGER = 1
const BINARY = 2
export CONTINUOUS, INTEGER, BINARY

###############################################################################
# Model class
# Keeps track of all model and column info
type Model
  obj#::QuadExpr
  objSense::Symbol
  
  linconstr#::Vector{LinearConstraint}
  quadconstr
  
  # Column data
  numCols::Int
  colNames::Vector{String}
  colLower::Vector{Float64}
  colUpper::Vector{Float64}
  colCat::Vector{Int}
  
  # Solution data
  objVal
  colVal::Vector{Float64}
  redCosts::Vector{Float64}
  linconstrDuals::Vector{Float64}
  # internal solver model object
  internalModel
  # Solver+option objects from MathProgBase
  lpsolver::AbstractMathProgSolver
  mipsolver::AbstractMathProgSolver
end

# Default constructor
function Model(sense::Symbol;lpsolver=MathProgBase.defaultLPsolver,mipsolver=MathProgBase.defaultMIPsolver)
  if (sense != :Max && sense != :Min)
     error("Model sense must be :Max or :Min")
  end
  Model(QuadExpr(),sense,LinearConstraint[], QuadConstraint[],
        0,String[],Float64[],Float64[],Int[],
        0,Float64[],Float64[],Float64[],nothing,lpsolver,mipsolver)
end

# Getters/setters
getNumVars(m::Model) = m.numCols
getNumConstraints(m::Model) = length(m.linconstr)
getObjectiveValue(m::Model) = m.objVal
getObjectiveSense(m::Model) = m.objSense
function setObjectiveSense(m::Model, newSense::Symbol)
  if (newSense != :Max && newSense != :Min)
    error("Model sense must be :Max or :Min")
  end
  m.objSense = newSense
end

# Pretty print
function print(io::IO, m::Model)
  println(io, string(m.objSense," ",quadToStr(m.obj)))
  println(io, "Subject to: ")
  for c in m.linconstr
    println(io, conToStr(c))
  end
  for c in m.quadconstr
    println(io, conToStr(c))
  end
  for i in 1:m.numCols
    print(io, m.colLower[i])
    print(io, " <= ")
    print(io, (m.colNames[i] == "" ? string("_col",i) : m.colNames[i]))
    print(io, " <= ")
    println(io, m.colUpper[i])
  end
end
show(io::IO, m::Model) = print(m.objSense == :Max ? "Maximization problem" :
                                                     "Minimization problem") 
                                                     # What looks good here?

###############################################################################
# Variable class
# Doesn't actually do much, just a pointer back to the model
type Variable
  m::Model
  col::Int
end

function Variable(m::Model,lower::Number,upper::Number,cat::Int,name::String)
  m.numCols += 1
  push!(m.colNames, name)
  push!(m.colLower, convert(Float64,lower))
  push!(m.colUpper, convert(Float64,upper))
  push!(m.colCat, cat)
  return Variable(m, m.numCols)
end

Variable(m::Model,lower::Number,upper::Number,cat::Int) =
  Variable(m,lower,upper,cat,"")

# Name setter/getters
setName(v::Variable,n::String) = (v.m.colNames[v.col] = n)
getName(v::Variable) = (v.m.colNames[v.col] == "" ? string("_col",v.col) : v.m.colNames[v.col])
getName(m::Model, col) = (m.colNames[col] == "" ? string("_col",col) : m.colNames[col])
print(io::IO, v::Variable) = print(io, getName(v))
show(io::IO, v::Variable) = print(io, getName(v))

# Bound setter/getters
setLower(v::Variable,lower::Number) = (v.m.colLower[v.col] = convert(Float64,lower))
setUpper(v::Variable,upper::Number) = (v.m.colUpper[v.col] = convert(Float64,upper))
getLower(v::Variable) = v.m.colLower[v.col]
getUpper(v::Variable) = v.m.colUpper[v.col]

# Value setter/getter
function setValue(v::Variable, val::Number)
    if length(v.m.colVal) < v.col
        resize!(v.m.colVal,v.m.numCols)
    end
    v.m.colVal[v.col] = val
end
getValue(v::Variable) = v.m.colVal[v.col]

# Dual value (reduced cost) getter
getDual(v::Variable) = v.m.redCosts[v.col]

###############################################################################
# Generic affine expression class
# Holds a vector of tuples (Var, Coeff)
type GenericAffExpr{CoefType,VarType}
  vars::Array{VarType,1}
  coeffs::Array{CoefType,1}
  constant::CoefType
end

typealias AffExpr GenericAffExpr{Float64,Variable}

AffExpr() = AffExpr(Variable[],Float64[],0.)

function setObjective(m::Model, a::AffExpr)
  m.obj = QuadExpr()
  m.obj.aff = a
end

print(io::IO, a::AffExpr) = print(io, affToStr(a))
show(io::IO, a::AffExpr) = print(io, affToStr(a))

function affToStr(a::AffExpr, showConstant=true)
  if length(a.vars) == 0
    return string(a.constant)
  end

  # Get reference to model
  m = a.vars[1].m

  # Collect like terms
  indvec = IndexedVector(Float64,m.numCols)
  for ind in 1:length(a.vars)
    addelt(indvec, a.vars[ind].col, a.coeffs[ind])
  end

  # Stringify the terms
  termStrings = Array(UTF8String, 2*length(a.vars)-1)
  if indvec.nnz > 0
    idx = indvec.nzidx[1]
    firstString = "$(indvec.elts[idx]) $(getName(m,idx))"
  end

  numTerms = 0
  for i in 2:indvec.nnz
    numTerms += 1
    idx = indvec.nzidx[i]
    if indvec.elts[idx] < 0
      termStrings[2*numTerms-1] = " - "
    else
      termStrings[2*numTerms-1] = " + "
    end
    termStrings[2*numTerms] = "$(abs(indvec.elts[idx])) $(getName(m,idx))"
  end

  # And then connect them up with +s
  # ret = join(termStrings[1:numTerms], " + ")
  ret = join([firstString, termStrings[1:(2*numTerms)]])
  
  if abs(a.constant) >= 0.000001 && showConstant
    if a.constant < 0
      ret = string(ret, " - ", abs(a.constant))
    else
      ret = string(ret, " + ", a.constant)
    end
  end
  return ret
end

###############################################################################
# QuadExpr class
# Holds a vector of tuples (Var, Var, Coeff), as well as an AffExpr
type QuadExpr
  qvars1::Vector{Variable}
  qvars2::Vector{Variable}
  qcoeffs::Vector{Float64}
  aff::AffExpr
end

QuadExpr() = QuadExpr(Variable[],Variable[],Float64[],AffExpr())

setObjective(m::Model, q::QuadExpr) = (m.obj = q)

print(io::IO, q::QuadExpr) = print(io, quadToStr(q))
show(io::IO, q::QuadExpr) = print(io, quadToStr(q))

function quadToStr(q::QuadExpr)
  if length(q.qvars1) == 0
    return affToStr(q.aff)
  end

  termStrings = Array(UTF8String, 2*length(q.qvars1))
  if length(q.qvars1) > 0
    if q.qcoeffs[1] < 0
      termStrings[1] = "-"
    else
      termStrings[1] = ""
    end
    for ind in 1:length(q.qvars1)
      if ind >= 2
        if q.qcoeffs[ind] < 0
          termStrings[2*ind-1] = " - "
        else 
          termStrings[2*ind-1] = " + "
        end
      end
      if q.qvars1[ind].col == q.qvars2[ind].col
        # Squared term
        termStrings[2*ind] = string(abs(q.qcoeffs[ind])," ",
                                    getName(q.qvars1[ind]),"²")
      else
        # Normal term
        termStrings[2*ind] = string(abs(q.qcoeffs[ind])," ",
                                    getName(q.qvars1[ind]),"*",
                                    getName(q.qvars2[ind]))
      end
    end
  end
  ret = join(termStrings)

  if q.aff.constant == 0 && length(q.aff.vars) == 0
    return ret
  else
    aff = affToStr(q.aff)
    if aff[1] == '-'
      return string(ret, " - ", aff[2:end])
    else
      return string(ret, " + ", aff)
    end
  end
end

##########################################################################
# JuMPConstraint
# abstract base for constraint types
abstract JuMPConstraint

##########################################################################
# LinearConstraint class
# An affine expression with lower bound (possibly -Inf) and upper bound (possibly Inf).
type LinearConstraint <: JuMPConstraint
  terms::AffExpr
  lb::Float64
  ub::Float64
end

LinearConstraint(terms::AffExpr,lb::Number,ub::Number) =
  LinearConstraint(terms,float(lb),float(ub))

function addConstraint(m::Model, c::LinearConstraint)
  push!(m.linconstr,c)
  return ConstraintRef{LinearConstraint}(m,length(m.linconstr))
end

print(io::IO, c::LinearConstraint) = print(io, conToStr(c))
show(io::IO, c::LinearConstraint) = print(io, conToStr(c))

function sense(c::LinearConstraint) 
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

function rhs(c::LinearConstraint)
  s = sense(c)
  @assert s != :range
  if s == :<=
    return c.ub
  else
    return c.lb
  end
end

function conToStr(c::LinearConstraint)
  s = sense(c)
  if s == :range
    return string(c.lb," <= ",affToStr(c.terms,false)," <= ",c.ub)
  else
    return string(affToStr(c.terms,false)," ",s," ",rhs(c))
  end
end

##########################################################################
# QuadConstraint class
# An quadratic constraint. Right-hand side is implicitly taken to be zero; 
# constraint is stored in the included QuadExpr.
type QuadConstraint <: JuMPConstraint
  terms::QuadExpr
  sense::Symbol
end

function addConstraint(m::Model, c::QuadConstraint)
  push!(m.quadconstr,c)
  return ConstraintRef{QuadConstraint}(m,length(m.quadconstr))
end

print(io::IO, c::QuadConstraint) = print(io, conToStr(c))
show(io::IO, c::QuadConstraint)  = print(io, conToStr(c))

conToStr(c::QuadConstraint) = string(quadToStr(c.terms), " ", c.sense, " 0")

##########################################################################
# ConstraintRef
# Reference to a constraint for retrieving solution info
immutable ConstraintRef{T<:JuMPConstraint}
  m::Model
  idx::Int
end

getDual(c::ConstraintRef{LinearConstraint}) = c.m.linconstrDuals[c.idx]

print(io::IO, c::ConstraintRef{LinearConstraint}) = print(io, conToStr(c.m.linconstr[c.idx]))
print(io::IO, c::ConstraintRef{QuadConstraint}) = print(io, conToStr(c.m.quadconstr[c.idx]))
show{T}(io::IO, c::ConstraintRef{T}) = print(io, c)

##########################################################################
# Operator overloads
include("operators.jl")
# Writers - we support MPS (MILP + QuadObj), LP (MILP)
include("writers.jl")
# Solvers
include("solvers.jl")
# Macros - @defVar, sum{}, etc.
include("macros.jl")

##########################################################################
end

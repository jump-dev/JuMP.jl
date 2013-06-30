###############################################################################
# MathProg 
# A MILP+QP modelling langauge for Julia
#
# By Iain Dunning and Miles Lubin
# http://www.github.com/IainNZ/MathProg.jl
###############################################################################

import Base.getindex
import Base.setindex!
import Base.print
import Base.show

module MathProg

# Use the standard solver interface for LPs and MIPs
using MathProgBase
require(joinpath(Pkg.dir("MathProgBase"),"src","LinprogSolverInterface.jl"))
using LinprogSolverInterface

importall Base

export
# Objects
  Model,
  Variable,
  AffExpr,
  QuadExpr,
  LinearConstraint,
  MultivarDict,

# Functions
  # Relevant to all
  print,show,
  # Model related
  getNumVars, getNumConstraints, getObjectiveValue, getObjective,
  setObjectiveSense, writeLP, writeMPS, setObjective, solve,
  addConstraint, addVar, addVars,
  # Variable
  setName, getName, setLower, setUpper, getLower, getUpper, getValue,
  # Expressions and constraints
  affToStr, quadToStr, conToStr,
  
# Macros and support functions
  @addConstraint,
  @defVar,
  @setObjective,
  addToExpression

include("MathProgDict.jl")
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
  
  # Column data
  numCols::Int
  colNames::Vector{String}
  colLower::Vector{Float64}
  colUpper::Vector{Float64}
  colCat::Vector{Int}
  
  # Solution data
  objVal
  colVal::Vector{Float64}
  # internal solver model object
  internalModel
  solverOptions
end

# Default constructor
function Model(sense::Symbol)
  if (sense != :Max && sense != :Min)
     error("Model sense must be :Max or :Min")
  end
  Model(QuadExpr(),sense,LinearConstraint[],
        0,String[],Float64[],Float64[],Int[],
        0,Float64[],nothing,Dict())
end

# Getters/setters
getNumVars(m::Model) = m.numCols
getNumConstraints(m::Model) = length(m.linconstr)
getObjectiveValue(m::Model) = m.objVal
getObjectiveSense(m::Model) = m.objSense
setObjectiveSense(m::Model, newSense:::Symbol) = (m.objSense = newSense)

# Pretty print
function print(io::IO, m::Model)
  println(io, string(m.objSense," ",quadToStr(m.obj)))
  println(io, "Subject to: ")
  for c in m.linconstr
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

# Value getter
getValue(v::Variable) = v.m.colVal[v.col]

###############################################################################
# Affine Expression class
# Holds a vector of tuples (Var, Coeff)
type AffExpr
  vars::Array{Variable,1}
  coeffs::Array{Float64,1}
  constant::Float64
end

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
  termStrings = Array(ASCIIString, length(a.vars))
  numTerms = 0
  for i in 1:indvec.nnz
    idx = indvec.nzidx[i]
    numTerms += 1
    termStrings[numTerms] = "$(indvec.elts[idx]) $(getName(m,idx))"
  end

  # And then connect them up with +s
  ret = join(termStrings[1:numTerms], " + ")
  
  if abs(a.constant) >= 0.000001 && showConstant
    ret = string(ret," + ",a.constant)
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

  termStrings = Array(ASCIIString, length(q.qvars1))
  for ind in 1:length(q.qvars1)
    termStrings[ind] = string(q.qcoeffs[ind]," ",
                              getName(q.qvars1[ind]),"*",
                              getName(q.qvars2[ind]))
  end
  ret = join(termStrings, " + ")

  if q.aff.constant == 0 && length(q.aff.vars) == 0
    return ret
  else
    return string(ret, " + ", affToStr(q.aff))
  end
end

##########################################################################
# LinearConstraint class
# An affine expression with lower bound (possibly -Inf) and upper bound (possibly Inf).
type LinearConstraint
  terms::AffExpr
  lb::Float64
  ub::Float64
end

LinearConstraint(terms::AffExpr,lb::Number,ub::Number) =
  LinearConstraint(terms,float(lb),float(ub))

addConstraint(m::Model, c::LinearConstraint) = push!(m.linconstr,c)

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

########################################################################
# Julp 
# A MILP modelling langauge for Julia
#  Julia
# +  LP
# -----
# =Julp.
#
# By Iain Dunning and Miles Lubin
########################################################################

module Julp

using Base
import Base.(+),Base.(-),Base.(*),Base.(<=),Base.(>=),Base.(==)
import Base.print

export
# Objects
  Model,
  Variable,
  AffExpr,
  Constraint,

# Functions
  print,exprToString,conToString,writeLP,
  setName,getName,setLower,setUpper,getLower,getUpper,
  addConstraint,setObjective,

# Macros
  @SumExpr

macro SumExpr(expr)
  local coefarr = Expr(:comprehension,convert(Vector{Any},[:(convert(Float64,$(expr.args[1].args[2]))),expr.args[2]]),Any)
  local vararr = Expr(:comprehension,convert(Vector{Any},[expr.args[1].args[3],expr.args[2]]),Any)
  esc(:(AffExpr($vararr,$coefarr,0.)))
end

########################################################################
# Model class
# Keeps track of all model and column info
type Model
  objective
  objSense
  
  constraints
  
  # Column data
  numCols
  colNames
  colLower
  colUpper
  colCat
end


# Default constructor
Model(sense::String) = Model(0,sense,Array(Constraint,0),
							0,Array(String,0),Array(Float64,0),Array(Float64,0),Array(Int,0))

# Pretty print
function print(m::Model)
  print(strcat(m.objSense," "))
  print(m.objective)
  println("")
  for c in m.constraints
    print("s.t. ")
    print(c)
    println("")
  end
  for i in 1:m.numCols
    print(m.colLower[i])
    print(" <= ")
    print((m.colNames[i] == "" ? strcat("_col",i) : m.colNames[i]))
    print(" <= ")
    println(m.colUpper[i])
  end
end

########################################################################
# Variable class
# Doesn't actually do much, just a pointer back to the model
type Variable
  m::Model
  col::Integer
end

# Constructor 1 - no name
function Variable(m::Model,lower::Number,upper::Number,cat::Int)
  return Variable(m,lower,upper,cat,"")
end

# Constructor 2 - with name
function Variable(m::Model,lower::Number,upper::Number,cat::Int,name::String)
  m.numCols += 1
  push(m.colNames,name)
  push(m.colLower,convert(Float64,lower))
  push(m.colUpper,convert(Float64,upper))
  push(m.colCat, cat)
  return Variable(m,m.numCols)
end

# Name setter/getters
function getName(v::Variable,n::String)
  v.m.colNames[v.col] = n
end
function getName(v::Variable)
  return (v.m.colNames[v.col] == "" ? strcat("_col",v.col) : v.m.colNames[v.col])
end

function show(io,v::Variable)
  print(io,getName(v))
end

# Bound setter/getters
function setLower(v::Variable,lower::Number)
  v.m.colLower[v.col] = convert(Float64,lower)
end
function setUpper(v::Variable,upper::Number)
  v.m.colUpper[v.col] = convert(Float64,upper)
end
function getLower(v::Variable)
  return v.m.colLower[v.col]
end
function getUpper(v::Variable)
  return v.m.colUpper[v.col]
end

########################################################################
# Affine Expression class
# Holds a vector of tuples (Variable, double)
type AffExpr
  vars::Array{Variable,1}
  coeffs::Array{Float64,1}
  constant::Float64
end

function setObjective(m::Model, a::AffExpr)
  m.objective = a
end

# Pretty printer
function print(a::AffExpr)
  for ind in 1:length(a.vars)
    print(a.coeffs[ind])
    print("*")
    print(getName(a.vars[ind]))
    print(" + ")
  end
  print(a.constant)
end

function exprToString(a::AffExpr)
  @assert length(a.vars) > 0
  seen = zeros(Bool,a.vars[1].m.numCols)
  nSeen = 0
  precomputedStrings = Array(ASCIIString,length(a.vars))
  for ind in 1:length(a.vars)
    thisstr = "$(a.coeffs[ind]) $(getName(a.vars[ind]))"
    precomputedStrings[nSeen+1] = thisstr
    # TODO: check if already seen this variable using seen array
    # if so, go back and update.
    # maybe rows should be cleaned somewhere before this?
    nSeen += 1
  end
  ret = join(precomputedStrings[1:nSeen]," + ")
  
  if abs(a.constant) >= 0.000001
    ret = strcat(ret," + ",a.constant)
  end
  return ret
end

#######################################################################
# Constraint class
# Basically an affine expression with a sense and two sides
# If we don't allow range constraints, we can just keep a
# side and a sense.
type Constraint
  lhs::AffExpr
  sense::String
end

function addConstraint(m::Model, c::Constraint)
  push(m.constraints,c)
end

# Pretty printer
function print(c::Constraint)
  print(c.lhs)
  print(strcat(" ",c.sense," 0"))
end

function conToString(c::Constraint)
  return strcat(exprToString(c.lhs-c.lhs.constant)," ",c.sense," ",-c.lhs.constant)
end

###########################################################
# Overloads
# Variable
# Variable--Variable
function (+)(lhs::Variable, rhs::Variable)
  return AffExpr([lhs,rhs], [1.,1.], 0.)
end
function (-)(lhs::Variable, rhs::Variable)
  return AffExpr([lhs,rhs], [1.,-1.], 0.)
end
# Variable--Number
function (+)(lhs::Variable, rhs::Number)
  return AffExpr([lhs], [1.], rhs)
end
function (-)(lhs::Variable, rhs::Number)
  return AffExpr([lhs], [1.], -rhs)
end
function (*)(lhs::Variable, rhs::Number)
  return AffExpr([lhs],[convert(Float64,rhs)], 0.0)
end
# Variable--AffExpr
function (+)(lhs::Variable, rhs::AffExpr)
  ret = AffExpr(copy(rhs.vars),copy(rhs.coeffs),rhs.constant)
  push(ret.vars, lhs)
  push(ret.coeffs, 1.)
  return ret
end
function (-)(lhs::Variable, rhs::AffExpr)
  ret = AffExpr(copy(rhs.vars),copy(rhs.coeffs),rhs.constant)
  for ind in 1:length(ret.coeffs)
    ret.coeffs *= -1
  end
  ret.constant *= -1
  push(ret.vars, lhs)
  push(ret.coeffs, 1.)
  return ret
end

# Number
# Number--Variable
function (*)(lhs::Number, rhs::Variable)
  return AffExpr([rhs],[convert(Float64,lhs)], 0.0)
end
# Number--AffExpr
function (+)(lhs::Number, rhs::AffExpr)
  ret = AffExpr(copy(rhs.vars),copy(rhs.coeffs),rhs.constant)
  ret.constant += lhs
  return ret
end
function (-)(lhs::Number, rhs::AffExpr)
  ret = AffExpr(copy(rhs.vars),copy(rhs.coeffs),rhs.constant)
  for ind in 1:length(ret.coeffs)
    ret.coeffs *= -1
  end
  ret.constant = lhs - ret.constant
  return ret
end

# AffExpr
# AffExpr--Variable
function (+)(lhs::AffExpr, rhs::Variable) 
  ret = AffExpr(copy(lhs.vars),copy(lhs.coeffs),lhs.constant)
  push(ret.vars, rhs)
  push(ret.coeffs, 1.)
  return ret 
end
function (-)(lhs::AffExpr, rhs::Variable) 
  ret = AffExpr(copy(lhs.vars),copy(lhs.coeffs),lhs.constant)
  push(ret.vars, rhs)
  push(ret.coeffs, -1.)
  return ret 
end
# AffExpr--Number
function (+)(lhs::AffExpr, rhs::Number)
  return AffExpr(copy(lhs.vars),copy(lhs.coeffs), lhs.constant+rhs)
end
function (-)(lhs::AffExpr, rhs::Number)
  return AffExpr(copy(lhs.vars),copy(lhs.coeffs), lhs.constant-rhs)
end
# AffExpr--AffExpr
function (+)(lhs::AffExpr, rhs::AffExpr)
  ret = AffExpr(copy(lhs.vars),copy(lhs.coeffs), lhs.constant)
  ret.vars = cat(1,ret.vars,rhs.vars)
  ret.coeffs = cat(1,ret.coeffs,rhs.coeffs)
  return ret
end

# Constraint
# Constraint--Number
function (<=)(lhs::AffExpr, rhs::Number)
  return Constraint(lhs-rhs,"<=")
end
function (==)(lhs::AffExpr, rhs::Number)
  return Constraint(lhs-rhs,"==")
end
function (>=)(lhs::AffExpr, rhs::Number)
  return Constraint(lhs-rhs,">=")
end

########################################################################
function writeLP(m::Model, fname::String)
  f = open(fname, "w")
  write(f,"\\*Julp-created LP*\\\n")
  
  # Objective
  if m.objSense == "max"
    write(f,"Maximize\n")
  else
    write(f,"Minimize\n")
  end
  objStr = exprToString(m.objective)
  write(f,strcat(" obj: ", objStr, "\n"))
  
  # Constraints
  write(f,"Subject To\n")
  conCount = 0
  tic()
  for c in m.constraints
    conCount += 1
    write(f,strcat(" c",conCount,": ", conToString(c),"\n"))
  end
  toc()
  print("In writing constraints\n")

  # Bounds
  write(f,"Bounds\n")
  for i in 1:m.numCols
    n = (m.colNames[i] == "" ? strcat("_col",i) : m.colNames[i])
    if abs(m.colLower[i]) > 0.00001
      write(f,strcat(" ",m.colLower[i]," <= ",n))
      if abs(m.colUpper[i] - 1e10) > 0.00001 # need a "infinite" number?
        write(f,strcat(" <= ",m.colUpper[i],"\n"))
      else
        write(f,"\n")
      end
    else
      if abs(m.colUpper[i] - 1e10) > 0.00001 # need a "infinite" number?
        write(f,strcat(" ",n," <= ",m.colUpper[i],"\n"))
      end
    end
  end
  
  # Integer - to do


  # Done
  write(f,"End\n")
  close(f)
end

###########################################################
end

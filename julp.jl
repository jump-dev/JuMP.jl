###########################################################
# Julp 
# A MILP modelling langauge for Julia
#  Julia
# +  LP
# -----
# =Julp.
#
# By Iain Dunning and Miles Lubin
###########################################################

macro SumExpr(expr)
  x = Expr(:comprehension,convert(Vector{Any},[:($(expr.args[1].args[3]), convert(Float64,$(expr.args[1].args[2])) ),expr.args[2] ]),Any)
  :(AffExpr($x))
end

module Julp

using Base
import Base.(+),Base.(-),Base.(*),Base.(<=),Base.(>=),Base.(==)

export
# Objects
  Model,
  Variable,
  AffExpr,
  Constraint,

# Functions
  PrintModel,
  SetName,GetName,SetLower,SetUpper,GetLower,GetUpper,
  PrintExpr,ExprToString,
  AddConstraint,PrintCon,ConToString,
  WriteLP
  
# Model class
# Keeps track of all model and column info
type Model
  names
  cols
  objective
  sense
  constraints #::Array{Constraint} - can't figure out how
              # to use before define
  colLower
  colUpper
end

# Default constructor
Model(sense::String) = Model(Array(String,0),0,0,sense,Array(Constraint,0),Array(Float64,0),Array(Float64,0))

# Pretty print
function PrintModel(m::Model)
  print(strcat(m.sense," "))
  PrintExpr(m.objective)
  println("")
  for c in m.constraints
    print("s.t. ")
    PrintCon(c)
    println("")
  end
  for i in 1:m.cols
    print(m.colLower[i])
    print(" <= ")
    print((m.names[i] == "" ? strcat("_col",i) : m.names[i]))
    print(" <= ")
    println(m.colUpper[i])
  end
end

###########################################################
# Variable class
# Doesn't actually do much, just a pointer back to the
# model effectively
type Variable
  m::Model
  col::Integer
end

# Constructor 1 - no name
function Variable(m::Model,lower::Number,upper::Number)
  m.cols += 1
  push(m.names,"")
  push(m.colLower,convert(Float64,lower))
  push(m.colUpper,convert(Float64,upper))
  return Variable(m,m.cols)
end

# Constructor 2 - with name
function Variable(m::Model,n::String,lower::Number,upper::Number)
  m.cols += 1
  push(m.names,n)
  push(m.colLower,convert(Float64,lower))
  push(m.colUpper,convert(Float64,upper))
  return Variable(m,m.cols)
end

# Name setter/getters
function SetName(v::Variable,n::String)
  v.m.names[v.col] = n
end
function GetName(v::Variable)
  #n = v.m.names[v.col]
  #if n==""
  #  n = strcat("_col",v.col)
  #end
  return (v.m.names[v.col] == "" ? strcat("_col",v.col) : v.m.names[v.col])
end

function show(io,v::Variable)
  print(io,GetName(v))
end

# Bound setter/getters
function SetLower(v::Variable,lower::Number)
  v.m.colLower[v.col] = convert(Float64,lower)
end
function SetUpper(v::Variable,upper::Number)
  v.m.colUpper[v.col] = convert(Float64,upper)
end
function GetLower(v::Variable)
  return v.m.colLower[v.col]
end
function GetUpper(v::Variable)
  return v.m.colUpper[v.col]
end

###########################################################
# Affine Expression class
# Holds a vector of tuples (Variable, double)
type AffExpr
  data::Array{(Variable,Float64),1}
  constant::Float64
end

# Basic constructor
function AffExpr(data::Array{(Variable,Float64),1})
  AffExpr(data,0.)
end

# Pretty printer
function PrintExpr(a::AffExpr)
  for pair in a.data
    print(pair[2])
    print("*")
    print(GetName(pair[1]))
    print(" + ")
  end
  print(a.constant)
end

function lpSum(terms::Array{AffExpr})
  numTerms = sum([length(t.data) for t in terms])
  data = Array((Variable,Float64),numTerms)
  pos = 1
  for t in terms
    for v in t.data
      data[pos] = v
      pos+=1
    end
  end
  return AffExpr(data,0.)
end

function ExprToString(a::AffExpr)
  @assert length(a.data) > 0
  seen = zeros(Bool,a.data[1][1].m.cols)
  nSeen = 0
  precomputedStrings = Array(ASCIIString,length(a.data))
  for pair in a.data
    thisstr = "$(pair[2]) $(GetName(pair[1]))"
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

###########################################################
# Constraint class
# Basically an affine expression with a sense and two sides
# If we don't allow range constraints, we can just keep a
# side and a sense.
type Constraint
  lhs::AffExpr
  sense::String
end

function AddConstraint(m::Model, c::Constraint)
  push(m.constraints,c)
end

# Pretty printer
function PrintCon(c::Constraint)
  PrintExpr(c.lhs)
  print(strcat(" ",c.sense," 0"))
  #print(ConToString(c))
end

function ConToString(c::Constraint)
  return strcat(ExprToString(c.lhs-c.lhs.constant)," ",c.sense," ",-c.lhs.constant)
end

###########################################################
# Overloads
# Variable
# Variable--Variable
function (+)(lhs::Variable, rhs::Variable)
  return AffExpr([(lhs,+1.0) (rhs,+1.0)], 0.0)
end
function (-)(lhs::Variable, rhs::Variable)
  return AffExpr([(lhs,+1.0) (rhs,-1.0)], 0.0)
end
# Variable--Number
function (+)(lhs::Variable, rhs::Number)
  return AffExpr([(lhs,+1.0)], +rhs)
end
function (-)(lhs::Variable, rhs::Number)
  return AffExpr([(lhs,+1.0)], -rhs)
end
function (*)(lhs::Variable, rhs::Number)
  return AffExpr([(lhs,convert(Float64,rhs))], 0.0)
end
# Variable--AffExpr
function (+)(lhs::Variable, rhs::AffExpr)
  ret = AffExpr(copy(rhs.data),rhs.constant)
  push(ret.data, (lhs,1.0))
  return ret
end

# Number
# Number--Variable
function (*)(lhs::Number, rhs::Variable)
  return AffExpr([(rhs,convert(Float64,lhs))], 0.0)
end
# Number--Number (LOL!)
# Number--AffExpr
function (+)(lhs::Number, rhs::AffExpr)
  return AffExpr(copy(rhs.data), rhs.constant+lhs)
end

# AffExpr
# AffExpr--Variable
function (+)(lhs::AffExpr, rhs::Variable) 
  ret = AffExpr(copy(lhs.data),lhs.constant)
  push(ret.data, (rhs,1.0))
  return ret
end
# AffExpr--Number
function (+)(lhs::AffExpr, rhs::Number)
  return AffExpr(copy(lhs.data), lhs.constant+rhs)
end
function (-)(lhs::AffExpr, rhs::Number)
  return AffExpr(copy(lhs.data), lhs.constant-rhs)
end
# AffExpr--AffExpr
function (+)(lhs::AffExpr, rhs::AffExpr)
  ret = AffExpr(copy(lhs.data),lhs.constant)
  ret.data = cat(1,ret.data,rhs.data)
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

###########################################################
function WriteLP(m::Model, fname::String)
  f = open(fname, "w")
  write(f,"\\*Julp-created LP*\\\n")
  
  # Objective
  if m.sense == "max"
    write(f,"Maximize\n")
  else
    write(f,"Minimize\n")
  end
  objStr = ExprToString(m.objective)
  write(f,strcat(" obj: ", objStr, "\n"))
  
  # Constraints
  write(f,"Subject To\n")
  conCount = 0
  tic()
  for c in m.constraints
    conCount += 1
    write(f,strcat(" c",conCount,": ", ConToString(c),"\n"))
  end
  toc()
  print("In writing constraints\n")

  # Bounds
  write(f,"Bounds\n")
  for i in 1:m.cols
    n = (m.names[i] == "" ? strcat("_col",i) : m.names[i])
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

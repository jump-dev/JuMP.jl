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

module Julp

import Base.*

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

function ExprToString(a::AffExpr)
  # This is "wrong" because it doesn't collect variables that
  # might appear multiple times e.g. 1x + 1x + 2y = 2x + 2y
  @assert length(a.data) > 0
  ret = "$(a.data[1][2]) $(GetName(a.data[1][1]))"
  if (mod(length(a.data),2) == 0)
    upto = length(a.data)-1
  else
    upto = length(a.data)
  end
  for i in 2:2:upto
    @assert i+1 <= length(a.data)
    pair1 = a.data[i]
    pair2 = a.data[i+1]
    ret = "$(ret) + $(pair1[2]) $(GetName(pair1[1])) + $(pair2[2]) $(GetName(pair2[1]))"
  end
  if upto < length(a.data)
    @assert upto + 1 == length(a.data)
    pair = a.data[length(a.data)]
    ret = "$(ret) + $(pair[2]) $(GetName(pair[1]))"
  end
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

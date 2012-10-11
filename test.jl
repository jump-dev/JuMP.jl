###########################################################
# JULP
###########################################################
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
    print(m.names[i])
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
  return v.m.names[v.col]
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
###########################################################



###########################################################
# Test script
m = Model("max")

x = Variable(m,"x",0,1e10)
y = Variable(m,"y",0,1e10)

m.objective = 5*x + 3*y
AddConstraint(m, 1.0*x + 5*y <= 3.0)
PrintModel(m)

exit()



affexpr = 5*x
PrintExpr(affexpr)
println("Changing variable name")
println(m)
SetName(y,"ydash")
println(m)

println("Adding two terms")
affexpr2 = 5*x + 3*y
PrintExpr(affexpr2)

println("Adding a third term")
z = Variable(m,"z")
affexpr3 = affexpr2 + z
PrintExpr(affexpr2) # Check we haven't modified affexpr2
PrintExpr(affexpr3)

affexpr3 = affexpr2 + 1*z
PrintExpr(affexpr3)


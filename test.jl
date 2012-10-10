#######################################
# Model class
# Keeps track of all column info
type Model
  names
  cols
end

# Default constructor
Model() = Model(Array(String,1),0)

###########################################################
# Variable class
# Doesn't actually do much, just a pointer back to the
# model effectively
type Variable
  m::Model
  col::Integer
end

#######################################
# Methods of the above classes
function Variable(m::Model)
  m.cols += 1
  if m.cols == 1
    m.names[1] = ""
  else
    push(m.names,"")
  end
  return Variable(m,m.cols)
end

function Variable(m::Model,n::String)
  m.cols += 1
  if m.cols == 1
    m.names[1] = n
  else
    push(m.names,"")
  end
  return Variable(m,m.cols)
end

function SetName(v::Variable,n::String)
  v.m.names[v.col] = n
end

function GetName(v::Variable)
  return v.m.names[v.col]
end

function PrintExpr(a::Array)
  for pair in a
    print(pair[2])
    print("*")
    print(GetName(pair[1]))
    print(" + ")
  end
  println("")
end

# Overloads
(*)(lhs::Variable,rhs::Number)   = [(lhs,rhs)]
(*)(lhs::Number,  rhs::Variable) = [(rhs,lhs)]
(+)(lhs::Array,   rhs::Array)    = hcat(lhs,rhs)
function (+)(lhs::Array,   rhs::Variable) 
  # Because it can't seem to cat in the 1st dimension
  # we ended up with a column vector above
  # But then push doesn't work...
  # No idea whats going on
  ret = copy(lhs)
  #show(ret)
  #push(ret,(rhs,1))
  #ret = hcat(ret,[(rhs,1)])
  return ret
end

# Test script
m = Model()

x = Variable(m,"x")
y = Variable(m)

affexpr = 5*x
PrintExpr(affexpr)


println("Changing variable name")
println(m)
SetName(y,"y")
println(m)

println("Adding two terms")
affexpr2 = 5*x + 3*y
PrintExpr(affexpr2)

println("Adding a third term")
z = Variable(m,"z")
affexpr3 = affexpr2 + z
PrintExpr(affexpr3)
# This doesn't work either
#affexpr3 = affexpr2 + 1*z
#PrintExpr(affexpr3)

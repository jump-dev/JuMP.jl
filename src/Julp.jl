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

importall Base

export
# Objects
  Model,
  Variable,
  AffExpr,
  Constraint,

# Functions
  print,exprToString,conToString,writeLP,writeMPS,
  setName,getName,setLower,setUpper,getLower,getUpper,
  addConstraint,setObjective,

# Macros and support functions
  @sumExpr,
  lpSum

include("macros.jl")

macro sumExprOld(expr)
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
  objIsQuad
  
  constraints
  
  # Column data
  numCols::Int
  colNames::Vector{String}
  colLower::Vector{Float64}
  colUpper::Vector{Float64}
  colCat::Vector{Int}
end


# Default constructor
Model(sense::String) = Model(0,sense,false,Array(Constraint,0),
							0,String[],Float64[],Float64[],Int[])

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
  push!(m.colNames,name)
  push!(m.colLower,convert(Float64,lower))
  push!(m.colUpper,convert(Float64,upper))
  push!(m.colCat, cat)
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

function AffExpr()
  return AffExpr(Variable[],Float64[],0.)
end

function setObjective(m::Model, a::AffExpr)
  m.objective = a
  m.objIsQuad = false
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
# QuadExpr class
type QuadExpr
  quadVars1::Vector{Variable}
  quadVars2::Vector{Variable}
  quadCoeffs::Vector{Float64}
  aff::AffExpr
end

function setObjective(m::Model, q::QuadExpr)
  m.objective = q
  m.objIsQuad = true
end

# Pretty printer
function print(q::QuadExpr)
  for ind in 1:length(q.quadVars1)
    print(q.quadCoeffs[ind])
    print("*")
    print(getName(q.quadVars1[ind]))
    print("*")
    print(getName(q.quadVars2[ind]))
    print(" + ")
  end
  print(q.aff)
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
  push!(m.constraints,c)
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
function (*)(lhs::Variable, rhs::Variable)
  return QuadExpr([lhs],[rhs], [1.], AffExpr(Variable[],Float64[],0.))
end
# Variable--Number
function (+)(lhs::Variable, rhs::Number)
  return AffExpr([lhs], [1.], convert(Float64,rhs))
end
function (-)(lhs::Variable, rhs::Number)
  return AffExpr([lhs], [1.], -convert(Float64,rhs))
end
function (*)(lhs::Variable, rhs::Number)
  return AffExpr([lhs],[convert(Float64,rhs)], 0.0)
end
# Variable--AffExpr
function (+)(lhs::Variable, rhs::AffExpr)
  ret = AffExpr(copy(rhs.vars),copy(rhs.coeffs),rhs.constant)
  push!(ret.vars, lhs)
  push!(ret.coeffs, 1.)
  return ret
end
function (-)(lhs::Variable, rhs::AffExpr)
  ret = AffExpr(copy(rhs.vars),copy(rhs.coeffs),rhs.constant)
  for ind in 1:length(ret.coeffs)
    ret.coeffs *= -1
  end
  ret.constant *= -1
  push!(ret.vars, lhs)
  push!(ret.coeffs, 1.)
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
  push!(ret.vars, rhs)
  push!(ret.coeffs, 1.)
  return ret 
end
function (-)(lhs::AffExpr, rhs::Variable) 
  ret = AffExpr(copy(lhs.vars),copy(lhs.coeffs),lhs.constant)
  push!(ret.vars, rhs)
  push!(ret.coeffs, -1.)
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
  ret = AffExpr(copy(lhs.vars), copy(lhs.coeffs), lhs.constant)
  ret.vars = cat(1,ret.vars,rhs.vars)
  ret.coeffs = cat(1,ret.coeffs,rhs.coeffs)
  return ret
end
function (*)(lhs::AffExpr, rhs::AffExpr)
  ret = QuadExpr(Variable[],Variable[],Float64[],AffExpr(Variable[],Float64[],0.))
  
  # Quadratic terms
  for ind1 = 1:length(lhs.coeffs)
    for ind2 = 1:length(rhs.coeffs)
      v1 = lhs.vars[ind1]
      v2 = rhs.vars[ind2]
      c  = lhs.coeffs[ind1]*rhs.coeffs[ind2]
      push!(ret.quadVars1,v1)
      push!(ret.quadVars2,v2)
      push!(ret.quadCoeffs,c)
    end
  end
  
  # Expr 1 constant * Expr 2 terms
  if lhs.constant != 0
    for ind2 = 1:length(rhs.coeffs)
      v2 = rhs.vars[ind2]
      c  = lhs.constant*rhs.coeffs[ind2]
      push!(ret.aff.vars,v2)
      push!(ret.aff.coeffs,c)
    end
    ret.aff.constant += lhs.constant*rhs.constant
  end
  
  # Expr 2 constant * Expr 1 terms
  if rhs.constant != 0
    for ind1 = 1:length(lhs.coeffs)
      v1 = lhs.vars[ind1]
      c  = rhs.constant*lhs.coeffs[ind1]
      push!(ret.aff.vars,v1)
      push!(ret.aff.coeffs,c)
    end
  end
  
  return ret
end

# QuadExpr
# QuadExpr--QuadExpr
function (+)(lhs::QuadExpr, rhs::QuadExpr)
  ret = QuadExpr(Variable[],Variable[],Float64[],AffExpr())
  ret.quadVars1 = cat(1,lhs.quadVars1,rhs.quadVars1)
  ret.quadVars2 = cat(1,lhs.quadVars2,rhs.quadVars2)
  ret.quadCoeffs = cat(1,lhs.quadCoeffs,rhs.quadCoeffs)
  ret.aff = lhs.aff + rhs.aff
  return ret
end

# Constraint
# Constraint--Number
function (<=)(lhs::AffExpr, rhs::Number)
  lhs.constant -= rhs #TODO: this is temporary, need to restructure
  return Constraint(lhs,"<=")
end
function (==)(lhs::AffExpr, rhs::Number)
  lhs.constant -= rhs
  return Constraint(lhs,"==")
end
function (>=)(lhs::AffExpr, rhs::Number)
  return Constraint(lhs-rhs,">=")
end

########################################################################
function writeMPS(m::Model, fname::String)
  f = open(fname, "w")
  
  numRows = length(m.constraints)
  # cache these strings
  #colnames = [ "x$(col)" for col in 1:m.numCols ]
  #colnames = [ convert(ASCIIString,@sprintf("x%d",col)) for col in 1:m.numCols ]
  # 
  #rownames = [ "CON$(i)" for i in 1:(numRows) ]
  #rownames = [ @sprintf("CON%d",i) for i in 1:(numRows) ]
  #push!(rownames, "obj")
  
  # Objective and constraint names
  write(f,"ROWS\n")
  write(f," N  CON$(numRows+1)\n")
  for c in 1:numRows
    senseChar = 'L'
    if m.constraints[c].sense == "=="
      senseChar = 'E'
    elseif m.constraints[c].sense == ">="
      senseChar = 'G'
    end
    #write(f,strcat(" ",senseChar,"  CON",c,"\n"))
    #@printf(f," %c  %s\n",senseChar,rownames[c])
    @printf(f," %c  CON%d\n",senseChar,c)
  end
  #tic() 
  # load rows into SparseMatrixCSC
  rowptr = Array(Int,numRows+2)
  nnz = 0
  for c in 1:numRows
      nnz += length(m.constraints[c].lhs.coeffs)
  end
  objaff = (m.objIsQuad) ? m.objective.aff : m.objective
  nnz += length(objaff.coeffs)
  colval = Array(Int,nnz)
  rownzval = Array(Float64,nnz)
  nnz = 0
  for c in 1:numRows
      rowptr[c] = nnz + 1
      coeffs = m.constraints[c].lhs.coeffs
      vars = m.constraints[c].lhs.vars
      for ind in 1:length(coeffs)
          nnz += 1
          colval[nnz] = vars[ind].col
          rownzval[nnz] = coeffs[ind]
      end
  end
  rowptr[numRows+1] = nnz + 1
  for ind in 1:length(objaff.coeffs)
      nnz += 1
      colval[nnz] = objaff.vars[ind].col
      rownzval[nnz] = objaff.coeffs[ind]
  end
  rowptr[numRows+2] = nnz + 1

  rowmat = SparseMatrixCSC(m.numCols,numRows+1, rowptr, colval, rownzval)
  colmat = rowmat'
  colptr = colmat.colptr
  rowval = colmat.rowval
  nzval = colmat.nzval
  #toc()
  #println("time building column structure")
    
  # Output each column
  write(f,"COLUMNS\n")
  for col in 1:m.numCols
    for ind in colmat.colptr[col]:(colmat.colptr[col+1]-1)
      #write(f,"    $(colnames[col])  $(rownames[rowval[ind]])  $(nzval[ind])\n")
      #@printf(f,"    %s  %s  %f\n",colnames[col],rownames[rowval[ind]],nzval[ind])
      @printf(f,"    x%d  CON%d  %f\n",col,rowval[ind],nzval[ind])
    end
  end
  
  # RHSs
  write(f,"RHS\n")
  for c in 1:numRows
    #write(f,"    rhs    CON$(c)    $(-m.constraints[c].lhs.constant)\n")
    #@printf(f,"    rhs    %s    %f\n",rownames[c],-m.constraints[c].lhs.constant)
    @printf(f,"    rhs    CON%d    %f\n",c,-m.constraints[c].lhs.constant)
  end
  
  # BOUNDS
  write(f,"BOUNDS\n")
  for col in 1:m.numCols
    if m.colLower[col] == 0 && m.colUpper[col] > 0
      # Default lower 0, and an upper
      #write(f,"  UP BOUND x$(col) $(m.colUpper[col])\n")
      #@printf(f,"  UP BOUND %s %f\n", colnames[col], m.colUpper[col])
      @printf(f,"  UP BOUND x%d %f\n", col, m.colUpper[col])
    elseif m.colLower[col] == -Inf && m.colUpper[col] == +Inf
      # Free
      #write(f,"  FR BOUND x$(col)\n")
      @printf(f, "  FR BOUND x%d\n", col)
    elseif m.colLower[col] != -Inf && m.colUpper[col] == +Inf
      # No upper, but a lower
      #write(f,"  PL BOUND x$(col) \n")
      #write(f,"  LO BOUND x$(col) $(m.colLower[col])\n")
      @printf(f, "  PL BOUND x%d\n  LO BOUND x%d %f\n",col,col,m.colLower[col])
    elseif m.colLower[col] == -Inf && m.colUpper[col] != +Inf
      # No lower, but a upper
      #write(f,"  MI BOUND x$(col) \n")
      #write(f,"  UP BOUND x$(col) $(m.colUpper[col])\n")
      @printf(f,"  MI BOUND x%d\n  UP BOUND x%d %f\n",col,col,m.colUpper[col])
    else
      # Lower and upper
      #write(f,"  LO BOUND x$(col) $(m.colLower[col])\n")
      #write(f,"  UP BOUND x$(col) $(m.colUpper[col])\n")
      @printf(f, "  LO BOUND x%d %f\n  UP BOUND x%d %f\n",col,col,m.colLower[col],m.colUpper[col])
    end
  end
  
  # Quadratic objective
  if m.objIsQuad
    write(f,"QMATRIX\n")
    qv1 = m.objective.quadVars1
    qv2 = m.objective.quadVars2
    qc  = m.objective.quadCoeffs
    for ind = 1:length(qv1)
      if qv1[ind].col == qv2[ind].col
        # Diagonal element
        #write(f,"  x$(qv1[ind].col)  x$(qv2[ind].col)  $(2*qc[ind])\n")
        @printf(f,"  x%d x%d  %f\n",qv1[ind].col,qv2[ind].col, 2qc[ind])
      else
        # Off diagonal, and we're gonna assume no duplicates
        #write(f,"  x$(qv1[ind].col)  x$(qv2[ind].col)  $(  qc[ind])\n")
        #write(f,"  x$(qv2[ind].col)  x$(qv1[ind].col)  $(  qc[ind])\n")
        @printf(f, "  x%d x%d %f\n", qv1[ind].col,qv2[ind].col, qc[ind])
        @printf(f, "  x%d x%d %f\n", qv2[ind].col,qv1[ind].col, qc[ind])
      end
    end
  end
  
  write(f,"ENDATA\n")
  close(f)
end

function writeLP(m::Model, fname::String)
  f = open(fname, "w")
  write(f,"NAME Julp-created LP \n")
  
  if m.objIsQuad
    print("Can't handle quad obj for LP yet\n")
    return
  end
  
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
  #tic()
  for c in m.constraints
    conCount += 1
    write(f,strcat(" c",conCount,": ", conToString(c),"\n"))
  end
  #toc()
  #print("In writing constraints\n")

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
function lpSum(expr)
  ret = AffExpr()
  for j in expr
	ret.vars   = cat(1, ret.vars,   j.vars)
	ret.coeffs = cat(1, ret.coeffs, j.coeffs)
  end
  return ret
end

###########################################################
end

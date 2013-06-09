########################################################################
# MathProg 
# A MILP+QP modelling langauge for Julia
#
# By Iain Dunning and Miles Lubin
# http://www.github.com/IainNZ/MathProg.jl
########################################################################

import Base.getindex
import Base.setindex!
import Base.print

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
  Constraint,
  MultivarDict,

# Functions
  print,exprToStr,conToString,writeLP,writeMPS,
  setName,getName,setLower,setUpper,getLower,getUpper,getValue,
  addConstraint,setObjective,solve,addVar,addVars,

# Macros and support functions
  @addConstraint,
  @defVar,
  @setObjective,
  addToExpression,
  lpSum

include("MathProgDict.jl")

########################################################################
# Constants
const CONTINUOUS = 0
const INTEGER = 1
const BINARY = 2
export CONTINUOUS, INTEGER, BINARY

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
  
  # Solution data
  objVal
  colVal::Vector{Float64}
  # internal solver model object
  internalModel
  solverOptions
end


# Default constructor
Model(sense::String) = Model(0,sense,false,Array(Constraint,0),
							0,String[],Float64[],Float64[],Int[],0,Float64[],nothing,Dict())

# Pretty print
function print(m::Model)
  print(string(m.objSense," "))
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
    print((m.colNames[i] == "" ? string("_col",i) : m.colNames[i]))
    print(" <= ")
    println(m.colUpper[i])
  end
end

########################################################################
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
print(io::IO, v::Variable) = print(io, getName(v))

# Bound setter/getters
setLower(v::Variable,lower::Number) = (v.m.colLower[v.col] = convert(Float64,lower))
setUpper(v::Variable,upper::Number) = (v.m.colUpper[v.col] = convert(Float64,upper))
getLower(v::Variable) = v.m.colLower[v.col]
getUpper(v::Variable) = v.m.colUpper[v.col]

# Value getter
getValue(v::Variable) = v.m.colVal[v.col]

########################################################################
# Affine Expression class
# Holds a vector of tuples (Variable, double)
type AffExpr
  vars::Array{Variable,1}
  coeffs::Array{Float64,1}
  constant::Float64
end

AffExpr() = AffExpr(Variable[],Float64[],0.)

setObjective(m::Model, a::AffExpr) = (m.objective = a)

print(io::IO, a::AffExpr) = print(io, exprToStr(a))

function exprToStr(a::AffExpr)
  if length(a.vars) == 0
    return string(a.constant)
  end
  seen = zeros(Bool,a.vars[1].m.numCols)
  nSeen = 0
  precomputedStrings = Array(ASCIIString,length(a.vars))
  for ind in 1:length(a.vars)
    thisstr = "$(a.coeffs[ind]) $(getName(a.vars[ind]))"
    precomputedStrings[nSeen+1] = thisstr
    # TODO: collect like terms
    nSeen += 1
  end
  ret = join(precomputedStrings[1:nSeen]," + ")
  
  if abs(a.constant) >= 0.000001
    ret = string(ret," + ",a.constant)
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

QuadExpr() = QuadExpr(Variable[],Variable[],Float64[],AffExpr())

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
  print(string(" ",c.sense," 0"))
end

function conToString(c::Constraint)
  return string(exprToString(c.lhs-c.lhs.constant)," ",c.sense," ",-c.lhs.constant)
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
(*)(lhs::Number, rhs::AffExpr) = (*)(rhs,lhs)
# Number--QuadExpr
(*)(lhs::Number, rhs::QuadExpr) = QuadExpr(copy(rhs.quadVars1),copy(rhs.quadVars2),lhs*rhs.quadCoeffs,lhs*rhs.aff)

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
function (/)(lhs::AffExpr, rhs::Number)
  return AffExpr(copy(lhs.vars), lhs.coeffs/rhs, lhs.constant/rhs)
end
function (*)(lhs::AffExpr, rhs::Number)
  return AffExpr(copy(lhs.vars), lhs.coeffs*rhs, lhs.constant*rhs)
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

  write(f,"NAME   MathProgModel\n")
  
  numRows = length(m.constraints)
  
  # Objective and constraint names 
  gc_disable()
  write(f,"ROWS\n")
  write(f," N  CON$(numRows+1)\n")
  for c in 1:numRows
    senseChar = 'L'
    if m.constraints[c].sense == "=="
      senseChar = 'E'
    elseif m.constraints[c].sense == ">="
      senseChar = 'G'
    end
    @printf(f," %c  CON%d\n",senseChar,c)
  end
  gc_enable()

  # Load rows into SparseMatrixCSC
  gc_disable()
  rowptr = Array(Int,numRows+2)
  nnz = 0
  for c in 1:numRows
      nnz += length(m.constraints[c].lhs.coeffs)
  end
  objaff::AffExpr = (m.objIsQuad) ? m.objective.aff : m.objective
  nnz += length(objaff.coeffs)
  colval = Array(Int,nnz)
  rownzval = Array(Float64,nnz)
  nnz = 0
  for c in 1:numRows
      rowptr[c] = nnz + 1
      # TODO: type assertion shouldn't be necessary
      constr::Constraint = m.constraints[c]
      coeffs = constr.lhs.coeffs
      vars = constr.lhs.vars
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
  gc_enable()
    
  # Output each column
  gc_disable()
  write(f,"COLUMNS\n")
  for col in 1:m.numCols
    for ind in colmat.colptr[col]:(colmat.colptr[col+1]-1)
      @printf(f,"    VAR%d  CON%d  %f\n",col,rowval[ind],nzval[ind])
    end
  end
  gc_enable()
  
  # RHSs
  gc_disable()
  write(f,"RHS\n")
  for c in 1:numRows
    @printf(f,"    rhs    CON%d    %f\n",c,-m.constraints[c].lhs.constant)
  end
  gc_enable()
  
  # BOUNDS
  gc_disable()
  write(f,"BOUNDS\n")
  for col in 1:m.numCols
    if m.colLower[col] == 0 && m.colUpper[col] > 0
      # Default lower 0, and an upper
      @printf(f,"  UP BOUND VAR%d %f\n", col, m.colUpper[col])
    elseif m.colLower[col] == -Inf && m.colUpper[col] == +Inf
      # Free
      @printf(f, "  FR BOUND VAR%d\n", col)
    elseif m.colLower[col] != -Inf && m.colUpper[col] == +Inf
      # No upper, but a lower
      @printf(f, "  PL BOUND VAR%d\n  LO BOUND VAR%d %f\n",col,col,m.colLower[col])
    elseif m.colLower[col] == -Inf && m.colUpper[col] != +Inf
      # No lower, but a upper
      @printf(f,"  MI BOUND VAR%d\n  UP BOUND VAR%d %f\n",col,col,m.colUpper[col])
    else
      # Lower and upper
      @printf(f, "  LO BOUND x%d %f\n  UP BOUND x%d %f\n",col,col,m.colLower[col],m.colUpper[col])
    end
  end
  gc_enable()
  
  # Quadratic objective
  gc_disable()
  if m.objIsQuad
    write(f,"QMATRIX\n")
    qv1 = m.objective.quadVars1
    qv2 = m.objective.quadVars2
    qc  = m.objective.quadCoeffs
    for ind = 1:length(qv1)
      if qv1[ind].col == qv2[ind].col
        # Diagonal element
        @printf(f,"  x%d x%d  %f\n",qv1[ind].col,qv2[ind].col, 2qc[ind])
      else
        # Off diagonal, and we're gonna assume no duplicates
        @printf(f, "  x%d x%d %f\n", qv1[ind].col,qv2[ind].col, qc[ind])
        @printf(f, "  x%d x%d %f\n", qv2[ind].col,qv1[ind].col, qc[ind])
      end
    end
  end
  
  write(f,"ENDATA\n")
  close(f)
  gc_enable()
end

function writeLP(m::Model, fname::String)
  f = open(fname, "w")
  write(f, "NAME Julp-created LP \n")
  
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
  objaff::AffExpr = (m.objIsQuad) ? m.objective.aff : m.objective
  write(f, " obj: ")
  nnz = length(objaff.coeffs)
  for ind in 1:(nnz-1)
    @printf(f, "%f VAR%d + ", objaff.coeffs[ind], objaff.vars[ind].col)
  end
  if nnz >= 1
    @printf(f, "%f VAR%d\n", objaff.coeffs[nnz], objaff.vars[nnz].col)
  end
  
  # Constraints
  #gc_disable()
  write(f,"Subject To\n")
  for i in 1:length(m.constraints)
    @printf(f, " c%d: ", i)

    c::Constraint = m.constraints[i]
    nnz = length(c.lhs.coeffs)
    for ind in 1:(nnz-1)
      @printf(f, "%f VAR%d + ", c.lhs.coeffs[ind], c.lhs.vars[ind].col)
    end
    if nnz >= 1
      @printf(f, "%f VAR%d", c.lhs.coeffs[nnz], c.lhs.vars[nnz].col)
    end
   
    # Sense and RHS
    if c.sense == "=="
      @printf(f, " = %f\n", -c.lhs.constant)
    elseif c.sense == "<="
      @printf(f, " <= %f\n", -c.lhs.constant)
    else
      @printf(f, " >= %f\n", -c.lhs.constant)
    end
  end
  #gc_enable()

  # Bounds
  #gc_disable()
  write(f,"Bounds\n")
  for i in 1:m.numCols    
    if m.colLower[i] == -Inf
      # No low bound
      if m.colUpper[i] == +Inf
        # Free
        @printf(f, " VAR%d free\n", i)
      else
        # x <= finite
        @printf(f, " -inf <= VAR%d <= %f\n", i, m.colUpper[i])
      end
    else
      # Low bound exists
      if m.colUpper[i] == +Inf
        # x >= finite
        @printf(f, " %f <= VAR%d <= +inf\n", m.colLower[i], i)
      else
        # finite <= x <= finite
        @printf(f, " %f <= VAR%d <= %f\n", m.colLower[i], i, m.colUpper[i])
      end
    end
  end
  #gc_enable()
  

  # Integer - to do


  # Done
  write(f,"End\n")
  close(f)
end

###########################################################
# Deprecated
function lpSum(expr)
  ret = AffExpr()
  for j in expr
	ret.vars   = cat(1, ret.vars,   j.vars)
	ret.coeffs = cat(1, ret.coeffs, j.coeffs)
  end
  return ret
end

###########################################################
# Solvers
function solve(m::Model)
	# Analyze model to see if any integers
	anyInts = false
	for j = 1:m.numCols
		if m.colCat[j] == INTEGER || m.colCat[j] == BINARY
			anyInts = true
			break
		end
	end
	
	if anyInts
		solveMIP(m)
	else
		solveLP(m)
	end
end

# prepare objective, constraint matrix, and row bounds
function prepProblem(m::Model)

    objaff::AffExpr = m.objIsQuad ? m.objective.aff : m.objective
    
    # We already have dense column lower and upper bounds

    # Create dense objective vector
    f = zeros(m.numCols)
    for ind in 1:length(objaff.vars)
        f[objaff.vars[ind].col] = objaff.coeffs[ind]
    end

    # Create row bounds
    numRows = length(m.constraints)
    rowlb = fill(-1e10, numRows)
    rowub = fill(+1e10, numRows)
    for c in 1:numRows
        if m.constraints[c].sense == "<="
            rowub[c] = -m.constraints[c].lhs.constant
        elseif m.constraints[c].sense == ">="
            rowlb[c] = -m.constraints[c].lhs.constant
        else
            rowub[c] = -m.constraints[c].lhs.constant
            rowlb[c] = -m.constraints[c].lhs.constant
        end
    end

    # Create sparse A matrix
    # First we build it row-wise, then use the efficient transpose
    # Theory is, this is faster than us trying to do it ourselves
    # Intialize storage
    rowptr = Array(Int,numRows+1)
    nnz = 0
    for c in 1:numRows
        nnz += length(m.constraints[c].lhs.coeffs)
    end
    colval = Array(Int,nnz)
    rownzval = Array(Float64,nnz)

    # Fill it up
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

    # Build the object
    rowmat = SparseMatrixCSC(m.numCols, numRows, rowptr, colval, rownzval)
    A = rowmat'

    return f, A, rowlb, rowub

end

function solveLP(m::Model)
    if m.objIsQuad
        error("Quadratic objectives are not fully supported yet")
    end

    f, A, rowlb, rowub = prepProblem(m)  

    # Ready to solve
    if MathProgBase.lpsolver == nothing
        error("No LP solver installed. Please run Pkg.add(\"Clp\") and restart Julia.")
    end
    m.internalModel = MathProgBase.lpsolver.model(;m.solverOptions...)
    loadproblem(m.internalModel, A, m.colLower, m.colUpper, f, rowlb, rowub)
    setsense(m.internalModel, m.objSense == "max" ? :Max : :Min)
    optimize(m.internalModel)
    stat = status(m.internalModel)

    if stat != :Optimal
        println("Warning: LP not solved to optimality, status: ", stat)
    else
        # store solution values in model
        m.objVal = getobjval(m.internalModel)
        m.objVal += m.objective.constant
        m.colVal = getsolution(m.internalModel)
    end

    return stat

end

function solveMIP(m::Model)
    if m.objIsQuad && string(MathProgBase.mipsolver) != "Gurobi"
        error("Quadratic objectives are not fully supported yet")
    end

    f, A, rowlb, rowub = prepProblem(m)

    # Build vartype vector
    vartype = zeros(Char,m.numCols)
    for j = 1:m.numCols
        if m.colCat[j] == CONTINUOUS
            vartype[j] = 'C'
        else
            vartype[j] = 'I'
        end
    end

    # Ready to solve
    if MathProgBase.mipsolver == nothing
        error("No MIP solver installed. Please run Pkg.add(\"CoinMP\") and restart Julia.")
    end
    
    m.internalModel = MathProgBase.mipsolver.model(;m.solverOptions...)
    # CoinMP doesn't support obj senses...
    if m.objSense == "max"
        f = -f
    end
    loadproblem(m.internalModel, A, m.colLower, m.colUpper, f, rowlb, rowub)
    setvartype(m.internalModel, vartype)
    # undocumented support for quadratic MIPs with gurobi:
    if m.objIsQuad
        gurobisolver = getrawsolver(m.internalModel)
        MathProgBase.mipsolver.add_qpterms!(gurobisolver, [v.col for v in m.objective.quadVars1], [v.col for v in m.objective.quadVars2], m.objective.quadCoeffs)
    end

    optimize(m.internalModel)
    stat = status(m.internalModel)

    if stat != :Optimal
        println("Warning: MIP not solved to optimality, status: ", stat)
    else
        # store solution values in model
        m.objVal = getobjval(m.internalModel)
        if m.objSense == "max"
            m.objVal = -m.objVal
        end
        m.objVal += !m.objIsQuad ? m.objective.constant : m.objective.aff.constant
        m.colVal = getsolution(m.internalModel)
    end

    return stat

end

include("macros.jl")

###########################################################
end

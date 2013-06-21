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
    qv1 = m.objective.qvars1
    qv2 = m.objective.qvars2
    qc  = m.objective.qcoeffs
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

###############################################################################
# LP File Writer
# We use the formatting defined at:
#   http://lpsolve.sourceforge.net/5.0/CPLEX-format.htm
function writeLP(m::Model, fname::String)

  f = open(fname, "w")

  # Coin's LP reader likes models to have a name
  write(f, "NAME Julp-created LP \n")
  
  if m.objIsQuad
    error("LP writer does not support quadratic objectives.\n")
  end
  
  # Objective
  if m.objSense == :Max
    write(f,"Maximize\n")
  else
    write(f,"Minimize\n")
  end
  objaff::AffExpr = m.objective
  write(f, " obj: ")
  nnz = length(objaff.coeffs)
  for ind in 1:(nnz-1)
    @printf(f, "%f VAR%d + ", objaff.coeffs[ind], objaff.vars[ind].col)
  end
  if nnz >= 1
    @printf(f, "%f VAR%d\n", objaff.coeffs[nnz], objaff.vars[nnz].col)
  end
  
  # Constraints
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

  # Bounds
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

  # Integer - don't handle binaries specially
  write(f,"General\n")
  for i in 1:m.numCols
    if m.colCat[i] != CONTINUOUS
      @printf(f, " VAR%d\n", i)
    end
  end

  # Done
  write(f,"End\n")
  close(f)
end


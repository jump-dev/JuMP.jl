function writeMPS(m::Model, fname::String)
  f = open(fname, "w")

  write(f,"NAME   MathProgModel\n")
  
  numRows = length(m.linconstr)
  
  # Objective and constraint names 
  gc_disable()
  write(f,"ROWS\n")
  write(f," N  CON$(numRows+1)\n")
  hasrange = false
  for c in 1:numRows
    rowsense = sense(m.linconstr[c])
    if rowsense == :(<=)
      senseChar = 'L'
    elseif rowsense == :(==)
      senseChar = 'E'
    elseif rowsense == :(>=)
      senseChar = 'G'
    else
      hasrange = true
      senseChar = 'E'
    end
    @printf(f," %c  CON%d\n",senseChar,c)
  end
  gc_enable()

  # Load rows into SparseMatrixCSC
  gc_disable()
  rowptr = Array(Int,numRows+2)
  nnz = 0
  for c in 1:numRows
      nnz += length(m.linconstr[c].terms.coeffs)
  end
  objaff::AffExpr = m.obj.aff
  objlincoef = objaff.coeffs
  if m.objSense == :Max
      println("Warning, MPS does not support maximization sense. Flipping objective coefficients.")
      objlincoef = -objaff.coeffs
  end


  nnz += length(objaff.coeffs)
  colval = Array(Int,nnz)
  rownzval = Array(Float64,nnz)
  nnz = 0
  for c in 1:numRows
      rowptr[c] = nnz + 1
      # TODO: type assertion shouldn't be necessary
      constr::LinearConstraint = m.linconstr[c]
      coeffs = constr.terms.coeffs
      vars = constr.terms.vars
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
      rownzval[nnz] = objlincoef[ind]
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
  inintegergroup = false
  write(f,"COLUMNS\n")
  for col in 1:m.numCols
    if m.colCat[col] != CONTINUOUS && !inintegergroup
      @printf(f,"    MARKER    'MARKER'                 'INTORG'\n")
      inintegergroup = true
    elseif m.colCat[col] == CONTINUOUS && inintegergroup
      @printf(f,"    MARKER    'MARKER'                 'INTEND'\n")
      inintegergroup = false
    end
    for ind in colmat.colptr[col]:(colmat.colptr[col+1]-1)
      @printf(f,"    VAR%d  CON%d  %f\n",col,rowval[ind],nzval[ind])
    end
  end
  if inintegergroup
    @printf(f,"    MARKER    'MARKER'                 'INTEND'\n")
  end
  gc_enable()
  
  # RHSs
  gc_disable()
  write(f,"RHS\n")
  for c in 1:numRows
    rowsense = sense(m.linconstr[c])
    if rowsense != :range
      @printf(f,"    rhs    CON%d    %f\n",c,rhs(m.linconstr[c]))
    else
      @printf(f,"    rhs    CON%d    %f\n",c,m.linconstr[c].lb)
    end
  end
  gc_enable()

  # RANGES
  if hasrange
    gc_disable()
    write(f,"RANGES\n")
    for c in 1:numRows
      rowsense = sense(m.linconstr[c])
      if rowsense == :range
        @printf(f,"    rhs    CON%d    %f\n",c,m.linconstr[c].ub-m.linconstr[c].lb)
      end
    end
  end

  
  # BOUNDS
  gc_disable()
  write(f,"BOUNDS\n")
  for col in 1:m.numCols
    if m.colLower[col] == 0
      if m.colUpper[col] != Inf
        # Default lower 0, and an upper
        @printf(f,"  UP BOUND VAR%d %f\n", col, m.colUpper[col])
      end
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
      @printf(f, "  LO BOUND VAR%d %f\n  UP BOUND VAR%d %f\n",col,m.colLower[col],col,m.colUpper[col])
    end
  end
  gc_enable()
  
  # Quadratic objective
  gc_disable()
  if length(m.obj.qvars1) != 0
    write(f,"QMATRIX\n")
    qv1 = m.obj.qvars1
    qv2 = m.obj.qvars2
    qc  = m.obj.qcoeffs
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

  if length(m.obj.qvars1) != 0
    error("LP writer does not support quadratic objectives.\n")
  end
  
  # Objective
  if m.objSense == :Max
    write(f,"Maximize\n")
  else
    write(f,"Minimize\n")
  end
  objaff::AffExpr = m.obj.aff
  write(f, " obj: ")
  nnz = length(objaff.coeffs)
  for ind in 1:(nnz-1)
    @printf(f, "%f VAR%d + ", objaff.coeffs[ind], objaff.vars[ind].col)
  end
  if nnz >= 1
    @printf(f, "%f VAR%d\n", objaff.coeffs[nnz], objaff.vars[nnz].col)
  end
  
  # Constraints
  function writeconstrterms(c::LinearConstraint)
    nnz = length(c.terms.coeffs)
    for ind in 1:(nnz-1)
      @printf(f, "%f VAR%d + ", c.terms.coeffs[ind], c.terms.vars[ind].col)
    end
    if nnz >= 1
      @printf(f, "%f VAR%d", c.terms.coeffs[nnz], c.terms.vars[nnz].col)
    end
  end
  write(f,"Subject To\n")
  constrcount = 1
  for i in 1:length(m.linconstr)
    @printf(f, " c%d: ", constrcount)

    c::LinearConstraint = m.linconstr[i]
    rowsense = sense(c)
    if rowsense != :range
      writeconstrterms(c)
      if rowsense == :(==)
        @printf(f, " = %f\n", rhs(c))
      elseif rowsense == :<=
        @printf(f, " <= %f\n", rhs(c))
      else 
        @assert rowsense == :>=
        @printf(f, " >= %f\n", rhs(c))
      end
      constrcount += 1
    else
      writeconstrterms(c)
      @printf(f, " >= %f\n", c.lb)
      @printf(f, " c%d: ", constrcount+1)
      writeconstrterms(c)
      @printf(f, " <= %f\n", c.ub)
      constrcount += 2
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


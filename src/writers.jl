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
        t = m.colCat[col]
        (t == :SemiCont || t == :SemiInt) && error("The MPS file writer does not currently support semicontinuous or semi-integer variables")
        if (t == :Bin || t == :Int) && !inintegergroup
            @printf(f,"    MARKER    'MARKER'                 'INTORG'\n")
            inintegergroup = true
        elseif (t == :Cont || t == :Fixed) && inintegergroup
            @printf(f,"    MARKER    'MARKER'                 'INTEND'\n")
            inintegergroup = false
        end
        for ind in colmat.colptr[col]:(colmat.colptr[col+1]-1)
            @printf(f,"    VAR%d  CON%d  ",col,rowval[ind])
            print_shortest(f,nzval[ind])
            println(f)
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
            @printf(f,"    rhs    CON%d    ",c)
            print_shortest(f,rhs(m.linconstr[c]))
        else
            @printf(f,"    rhs    CON%d    ",c)
            print_shortest(f,m.linconstr[c].lb)
        end
        println(f)
    end
    gc_enable()

    # RANGES
    if hasrange
        gc_disable()
        write(f,"RANGES\n")
        for c in 1:numRows
            rowsense = sense(m.linconstr[c])
            if rowsense == :range
                @printf(f,"    rhs    CON%d    ",c)
                print_shortest(f,m.linconstr[c].ub-m.linconstr[c].lb)
                println(f)
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
                @printf(f,"  UP BOUND VAR%d ", col)
                print_shortest(f, m.colUpper[col])
                println(f)
            end
        elseif m.colLower[col] == -Inf && m.colUpper[col] == +Inf
            # Free
            @printf(f, "  FR BOUND VAR%d\n", col)
        elseif m.colLower[col] != -Inf && m.colUpper[col] == +Inf
            # No upper, but a lower
            @printf(f, "  PL BOUND VAR%d\n  LO BOUND VAR%d ",col,col)
            print_shortest(f,m.colLower[col])
            println(f)
        elseif m.colLower[col] == -Inf && m.colUpper[col] != +Inf
            # No lower, but a upper
            @printf(f,"  MI BOUND VAR%d\n  UP BOUND VAR%d ",col,col)
            print_shortest(f,m.colUpper[col])
            println(f)
        else
            # Lower and upper
            @printf(f, "  LO BOUND VAR%d ",col)
            print_shortest(f,m.colLower[col])
            println(f)
            @printf(f, "  UP BOUND VAR%d ",col)
            print_shortest(f,m.colUpper[col])
            println(f)
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
                @printf(f,"  VAR%d VAR%d  ", qv1[ind].col,qv2[ind].col)
                print_shortest(f,2qc[ind])
                println(f)
            else
                # Off diagonal, and we're gonna assume no duplicates
                @printf(f, "  VAR%d VAR%d ", qv1[ind].col,qv2[ind].col)
                print_shortest(f, qc[ind])
                println(f)
                @printf(f, "  VAR%d VAR%d ", qv2[ind].col,qv1[ind].col)
                print_shortest(f, qc[ind])
                println(f)
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
        print_shortest(f, objaff.coeffs[ind])
        @printf(f, " VAR%d + ", objaff.vars[ind].col)
    end
    if nnz >= 1
        print_shortest(f, objaff.coeffs[nnz])
        @printf(f, " VAR%d\n", objaff.vars[nnz].col)
    end
    
    # Constraints
    function writeconstrterms(c::LinearConstraint)
        nnz = length(c.terms.coeffs)
        for ind in 1:(nnz-1)
            print_shortest(f, c.terms.coeffs[ind])
            @printf(f, " VAR%d + ", c.terms.vars[ind].col)
        end
        if nnz >= 1
            print_shortest(f, c.terms.coeffs[nnz])
            @printf(f, " VAR%d", c.terms.vars[nnz].col)
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
                @printf(f, " = ")
                print_shortest(f, rhs(c))
                println(f)
            elseif rowsense == :<=
                @printf(f, " <= ")
                print_shortest(f, rhs(c))
                println(f)
            else 
                @assert rowsense == :>=
                @printf(f, " >= ")
                print_shortest(f, rhs(c))
                println(f)
            end
            constrcount += 1
        else
            writeconstrterms(c)
            @printf(f, " >= ")
            print_shortest(f, c.lb)
            println(f)
            @printf(f, " c%d: ", constrcount+1)
            writeconstrterms(c)
            @printf(f, " <= ")
            print_shortest(f, c.ub)
            println(f)
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
                @printf(f, " -inf <= VAR%d <= ", i)
                print_shortest(f, m.colUpper[i])
                println(f)
            end
        else
            # Low bound exists
            if m.colUpper[i] == +Inf
                # x >= finite
                @printf(f, " ")
                print_shortest(f, m.colLower[i])
                @printf(f," <= VAR%d <= +inf\n", i)
            else
                # finite <= x <= finite
                @printf(f, " ")
                print_shortest(f, m.colLower[i])
                @printf(f, " <= VAR%d <= ", i)
                print_shortest(f, m.colUpper[i])
                println(f)
            end
        end
    end

    # Integer - don't handle binaries specially
    write(f,"General\n")
    for i in 1:m.numCols
        t = m.colCat[i]
        (t == :SemiCont || t == :SemiInt) && error("The LP file writer does not currently support semicontinuous or semi-integer variables")
        if t == :Bin || t == :Int
            @printf(f, " VAR%d\n", i)
        end
    end

    # Done
    write(f,"End\n")
    close(f)
end


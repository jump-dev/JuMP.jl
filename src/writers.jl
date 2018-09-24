#  Copyright 2017, Iain Dunning, Joey Huchette, Miles Lubin, and contributors
#  This Source Code Form is subject to the terms of the Mozilla Public
#  License, v. 2.0. If a copy of the MPL was not distributed with this
#  file, You can obtain one at http://mozilla.org/MPL/2.0/.
if VERSION ≥ v"0.7-"
    const print_shortest = Base.Grisu.print_shortest
end

function writeMPS(m::Model, fname::AbstractString)
    f = open(fname, "w")

    write(f,"NAME   JuMPModel\n")

    numRows = length(m.linconstr)

    # Objective and constraint names
    Compat.GC.enable(false)
    write(f,"ROWS\n")
    write(f," N  OBJ\n")
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
    Compat.GC.enable(true)

    objlincoef = prepAffObjective(m)
    rowlb, rowub = prepConstrBounds(m)
    if m.objSense == :Max
        println("Warning, MPS does not support maximization sense. Flipping objective coefficients.")
        objlincoef = -objlincoef
    end

    A = prepConstrMatrix(m)
    colptr = A.colptr
    rowval = A.rowval
    nzval = A.nzval

    # Output each column
    Compat.GC.enable(false)
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
        for ind in colptr[col]:(colptr[col+1]-1)
            @printf(f,"    VAR%d  CON%d  ",col,rowval[ind])
            print_shortest(f,nzval[ind])
            println(f)
        end
        @printf(f,"    VAR%d  OBJ  ",col)
        print_shortest(f,objlincoef[col])
        println(f)
    end
    if inintegergroup
        @printf(f,"    MARKER    'MARKER'                 'INTEND'\n")
    end
    Compat.GC.enable(true)

    # RHSs
    Compat.GC.enable(false)
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
    Compat.GC.enable(true)

    # RANGES
    if hasrange
        Compat.GC.enable(false)
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
    Compat.GC.enable(false)
    write(f,"BOUNDS\n")
    for col in 1:m.numCols
        if m.colLower[col] == 0
            if m.colUpper[col] != Inf
                # Default lower 0, and an upper
                @printf(f,"  UP BOUND VAR%d ", col)
                print_shortest(f, m.colUpper[col])
                println(f)
            else
                # Default bounds. Explicitly state for solvers like Gurobi. See issue #792
                @printf(f,"  PL BOUND VAR%d", col)
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
    Compat.GC.enable(true)

    # Quadratic objective
    Compat.GC.enable(false)
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
                #@printf(f, "  VAR%d VAR%d ", qv2[ind].col,qv1[ind].col)
                #print_shortest(f, qc[ind])
                #println(f)
            end
        end
    end

    write(f,"ENDATA\n")
    close(f)
    Compat.GC.enable(true)
end

###############################################################################
# LP File Writer
# We use the formatting defined at:
#   http://lpsolve.sourceforge.net/5.0/CPLEX-format.htm

varname_generic(m::Model, col::Integer) = "VAR$(col)"

function varname_given(m::Model, col::Integer)
    # TODO: deal with non-ascii characters?
    name = getname(m, col)
    for (pat, sub) in [("[", "_"), ("]", ""), (",", "_")]
        name = replace(name, pat => sub)
    end
    name
end

function writeLP(m::Model, fname::AbstractString; genericnames=true)
    varname = genericnames ? varname_generic : varname_given

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
        if ind == 1
            print_shortest(f, objaff.coeffs[ind])
        else
            print_shortest(f, abs(objaff.coeffs[ind]))
        end
        @printf(f, " %s %s ", varname(m, objaff.vars[ind].col), (objaff.coeffs[ind+1] < 0) ? "-" : "+")
    end
    if nnz >= 1
        if nnz == 1
            print_shortest(f, objaff.coeffs[nnz])
        else
            print_shortest(f, abs(objaff.coeffs[nnz]))
        end
        @printf(f, " %s\n", varname(m, objaff.vars[nnz].col))
    end

    # Constraints
    function writeconstrterms(c::LinearConstraint)
        nnz = length(c.terms.coeffs)
        for ind in 1:(nnz-1)
            if ind == 1
                print_shortest(f, c.terms.coeffs[ind])
            else
                print_shortest(f, abs(c.terms.coeffs[ind]))
            end
            @printf(f, " %s %s ", varname(m, c.terms.vars[ind].col), (c.terms.coeffs[ind+1] < 0) ? "-" : "+")
        end
        if nnz >= 1
            if nnz == 1
                print_shortest(f, c.terms.coeffs[nnz])
            else
                print_shortest(f, abs(c.terms.coeffs[nnz]))
            end
            @printf(f, " %s", varname(m, c.terms.vars[nnz].col))
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
    # SOS constraints
    for i in 1:length(m.sosconstr)
        @printf(f, " c%d: ", constrcount)

        c::SOSConstraint = m.sosconstr[i]
        if c.sostype == :SOS1
            @printf(f, "S1::")
        elseif  c.sostype == :SOS2
            @printf(f, "S2::")
        end

        @assert length(c.terms) == length(c.weights)
        for j in 1:length(c.terms)
            @printf(f, " %s:", varname(m, c.terms[j].col))
            print_shortest(f, c.weights[j])
        end

        println(f)
        constrcount += 1
    end

    # Bounds
    write(f,"Bounds\n")
    for i in 1:m.numCols
        if m.colLower[i] == -Inf
            # No low bound
            if m.colUpper[i] == +Inf
                # Free
                @printf(f, " %s free\n", varname(m, i))
            else
                # x <= finite
                @printf(f, " -inf <= %s <= ", varname(m, i))
                print_shortest(f, m.colUpper[i])
                println(f)
            end
        else
            # Low bound exists
            if m.colUpper[i] == +Inf
                # x >= finite
                @printf(f, " ")
                print_shortest(f, m.colLower[i])
                @printf(f," <= %s <= +inf\n", varname(m, i))
            else
                # finite <= x <= finite
                @printf(f, " ")
                print_shortest(f, m.colLower[i])
                @printf(f, " <= %s <= ", varname(m, i))
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
            @printf(f, " %s\n", varname(m, i))
        end
    end

    # Done
    write(f,"End\n")
    close(f)
end

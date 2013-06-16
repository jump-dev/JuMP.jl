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

    objaff::AffExpr = m.objective # TODO check
    
    # We already have dense column lower and upper bounds

    # Create dense objective vector
    f = zeros(m.numCols)
    for ind in 1:length(objaff.vars)
        f[objaff.vars[ind].col] = objaff.coeffs[ind]
    end

    # Create row bounds
    numRows = length(m.constraints)
    rowlb = fill(-Inf, numRows)
    rowub = fill(+Inf, numRows)
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
    tmprow = IndexedVector(Float64,m.numCols)
    tmpelts = tmprow.elts
    tmpnzidx = tmprow.nzidx
    for c in 1:numRows
        rowptr[c] = nnz + 1
        coeffs = m.constraints[c].lhs.coeffs
        vars = m.constraints[c].lhs.vars
        # collect duplicates
        for ind in 1:length(coeffs)
            addelt(tmprow,vars[ind].col,coeffs[ind])
        end
        for i in 1:tmprow.nnz
            nnz += 1
            idx = tmpnzidx[i]
            colval[nnz] = idx
            rownzval[nnz] = tmpelts[idx]
        end
        empty!(tmprow)
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
        MathProgBase.mipsolver.add_qpterms!(gurobisolver, [v.col for v in m.quadobj.qvars1], [v.col for v in m.quadobj.qvars2], m.quadobj.qcoeffs)
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
        m.objVal += m.objective.constant
        m.colVal = getsolution(m.internalModel)
    end

    return stat

end


setLPSolver(s::Symbol) = MathProgBase.setlpsolver(s)
setMIPSolver(s::Symbol) = MathProgBase.setmipsolver(s)

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

function quadraticGurobi(m::Model, solvermodule, ismip = false)
    # ugly hack for now until we get CoinMP to support setting objective senses
    doflip = false
    if ismip && m.objSense == :Max
        doflip = true
    end

    if length(m.obj.qvars1) != 0
        gurobisolver = getrawsolver(m.internalModel)
        solvermodule.add_qpterms!(gurobisolver, [v.col for v in m.obj.qvars1], [v.col for v in m.obj.qvars2], !doflip ? m.obj.qcoeffs : -m.obj.qcoeffs)
    end

# Add quadratic constraint to solver
    for k in 1:length(m.quadconstr)
        qconstr = m.quadconstr[k]
        gurobisolver = getrawsolver(m.internalModel)
        if !((s = string(qconstr.sense)[1]) in ['<', '>', '='])
            error("Invalid sense for quadratic constraint")
        end

        solvermodule.add_qconstr!(gurobisolver, 
                                  [v.col for v in qconstr.terms.aff.vars], 
                                  qconstr.terms.aff.coeffs, 
                                  [v.col for v in qconstr.terms.qvars1], 
                                  [v.col for v in qconstr.terms.qvars2], 
                                  qconstr.terms.qcoeffs, 
                                  s, 
                                  -qconstr.terms.aff.constant)
    end

    if length(m.quadconstr) > 0
        solvermodule.update_model!(gurobisolver)
    end
end

# prepare objective, constraint matrix, and row bounds
function prepProblem(m::Model)

    objaff::AffExpr = m.obj.aff
    
    # We already have dense column lower and upper bounds

    # Create dense objective vector
    f = zeros(m.numCols)
    for ind in 1:length(objaff.vars)
        f[objaff.vars[ind].col] += objaff.coeffs[ind]
    end

    # Create row bounds
    numRows = length(m.linconstr)
    rowlb = fill(-Inf, numRows)
    rowub = fill(+Inf, numRows)
    for c in 1:numRows
        rowlb[c] = m.linconstr[c].lb
        rowub[c] = m.linconstr[c].ub
    end

    # Create sparse A matrix
    # First we build it row-wise, then use the efficient transpose
    # Theory is, this is faster than us trying to do it ourselves
    # Intialize storage
    rowptr = Array(Int,numRows+1)
    nnz = 0
    for c in 1:numRows
        nnz += length(m.linconstr[c].terms.coeffs)
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
        coeffs = m.linconstr[c].terms.coeffs
        vars = m.linconstr[c].terms.vars
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
    f, A, rowlb, rowub = prepProblem(m)  

    # Ready to solve
    if MathProgBase.lpsolver == nothing
        error("No LP solver installed. Please run Pkg.add(\"Clp\") and restart Julia.")
    end

    m.internalModel = MathProgBase.lpsolver.model(;m.solverOptions...)
    loadproblem(m.internalModel, A, m.colLower, m.colUpper, f, rowlb, rowub)

    if length(m.obj.qvars1) != 0 || length(m.quadconstr) != 0
        if string(MathProgBase.lpsolver) != "Gurobi"
            error("Quadratic objectives/constraints are currently only supported using Gurobi")
        end
        quadraticGurobi(m, MathProgBase.lpsolver)
    end

    setsense(m.internalModel, m.objSense)

    optimize(m.internalModel)
    stat = status(m.internalModel)

    if stat != :Optimal
        println("Warning: LP not solved to optimality, status: ", stat)
    else
        # store solution values in model
        m.objVal = getobjval(m.internalModel)
        m.objVal += m.obj.aff.constant
        m.colVal = getsolution(m.internalModel)
        m.redCosts = getreducedcosts(m.internalModel)
        m.linconstrDuals = getconstrduals(m.internalModel)
    end

    return stat

end

function solveMIP(m::Model)
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
        error("No MIP solver installed. Please run Pkg.add(\"Cbc\") and restart Julia.")
    end
    
    m.internalModel = MathProgBase.mipsolver.model(;m.solverOptions...)
    # CoinMP doesn't support obj senses...
    if m.objSense == :Max
        f = -f
    end
    loadproblem(m.internalModel, A, m.colLower, m.colUpper, f, rowlb, rowub)
    setvartype(m.internalModel, vartype)

    if length(m.obj.qvars1) != 0 || length(m.quadconstr) != 0
        if string(MathProgBase.mipsolver) != "Gurobi"
            error("Quadratic objectives/constraints are currently only supported using Gurobi")
        end
        quadraticGurobi(m, MathProgBase.mipsolver, true)
    end

    optimize(m.internalModel)
    stat = status(m.internalModel)

    if stat != :Optimal
        println("Warning: MIP not solved to optimality, status: ", stat)
    end
    # It's possible that we have a feasible solution if we're not optimal
    # TODO: Test this behavior on various solvers
    try
        # store solution values in model
        m.objVal = getobjval(m.internalModel)
        if m.objSense == :Max
            m.objVal = -m.objVal
        end
        m.objVal += m.obj.aff.constant
        m.colVal = getsolution(m.internalModel)
    end

    return stat

end


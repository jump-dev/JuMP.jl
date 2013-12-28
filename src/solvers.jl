if Pkg.installed("Gurobi") != nothing
  eval(Expr(:using,:Gurobi))
end

if Pkg.installed("CPLEXLink") != nothing
  eval(Expr(:using,:CPLEXLink))
end

function solve(m::Model)
  # Analyze model to see if any integers
  anyInts = false
  for j = 1:m.numCols
    if m.colCat[j] == INTEGER
      anyInts = true
      break
    end
  end
	
  if anyInts
    if isa(m.solver,MathProgBase.MissingSolver)
      m.solver = MathProgBase.defaultMIPsolver
      s = solveMIP(m)
      # Clear solver in case we change problem types
      m.solver = MathProgBase.MissingSolver("",Symbol[])
      return s
    else
      solveMIP(m)
    end
  else
    if isa(m.solver,MathProgBase.MissingSolver)
      m.solver = MathProgBase.defaultLPsolver
      s = solveLP(m)
      m.solver = MathProgBase.MissingSolver("",Symbol[])
      return s
    else
      solveLP(m)
    end
  end
end

function quadCheck(m::Model, ismip = false)
    if length(m.obj.qvars1) != 0 || length(m.quadconstr) != 0
        if isdefined(:(GurobiSolver))
            GurobiCheck  = isa(m.solver,GurobiSolver)
        else
            GurobiCheck  = false
        end
        if isdefined(:(CplexSolver))
            CplexCheck  = isa(m.solver,CplexSolver)
        else
            CplexCheck  = false
        end
        if !GurobiCheck && !CplexCheck
            error("Quadratic objectives/constraints are currently only supported using Gurobi or CPLEXLink")
        end
        if !ismip && isa(m.solver,GurobiSolver)
            # Gurobi by default will not compute duals
            # if quadratic constraints are present.
            push!(m.solver.options,(:QCPDual,1))
        end
        return true
    end
    return false
end

function addQuadratics(m::Model)

    if length(m.obj.qvars1) != 0
      setquadobj!(m.internalModel, Cint[v.col for v in m.obj.qvars1], Cint[v.col for v in m.obj.qvars2], m.obj.qcoeffs)
    end

# Add quadratic constraint to solver
    for k in 1:length(m.quadconstr)
        qconstr = m.quadconstr[k]
        if !((s = string(qconstr.sense)[1]) in ['<', '>', '='])
          error("Invalid sense for quadratic constraint")
        end
        addquadconstr!(m.internalModel, Cint[v.col for v in qconstr.terms.aff.vars], qconstr.terms.aff.coeffs, Cint[v.col for v in qconstr.terms.qvars1], Cint[v.col for v in qconstr.terms.qvars2], qconstr.terms.qcoeffs, s, -qconstr.terms.aff.constant)
    end
    writeproblem(m.internalModel, "out.lp")
end

# prepare objective, constraint matrix, and row bounds
function prepProblemBounds(m::Model)

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
    return f, rowlb, rowub
end

# prepare column-wise constraint matrix
function prepConstrMatrix(m::Model)

    # Create sparse A matrix
    # First we build it row-wise, then use the efficient transpose
    # Theory is, this is faster than us trying to do it ourselves
    # Intialize storage
    numRows = length(m.linconstr)
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
end

function solveLP(m::Model)
    f, rowlb, rowub = prepProblemBounds(m)  

    # Ready to solve

    if !m.firstsolve
        try
            setvarLB!(m.internalModel, m.colLower)
            setvarUB!(m.internalModel, m.colUpper)
            setconstrLB!(m.internalModel, rowlb)
            setconstrUB!(m.internalModel, rowub)
            setobj!(m.internalModel, f)
        catch
            warn("LP solver does not appear to support hot-starts. Problem will be solved from scratch.")
            m.firstsolve = true
        end
    end
    if m.firstsolve
        A = prepConstrMatrix(m)
        hasQuad = quadCheck(m)
        m.internalModel = model(m.solver)
        loadproblem!(m.internalModel, A, m.colLower, m.colUpper, f, rowlb, rowub, m.objSense)

        if hasQuad
            addQuadratics(m)
        end
    end 

    optimize!(m.internalModel)
    stat = status(m.internalModel)

    if stat != :Optimal
        warn("LP not solved to optimality, status: $stat")
        if stat == :Infeasible
            try
                m.linconstrDuals = getinfeasibilityray(m.internalModel)
            catch
                println("Infeasibility ray (Farkas proof) not available")
            end
        elseif stat == :Unbounded
            try
                m.colVal = getunboundedray(m.internalModel)
            catch
                println("Unbounded ray not available")
            end
        end
    else
        # store solution values in model
        m.objVal = getobjval(m.internalModel)
        m.objVal += m.obj.aff.constant
        m.colVal = getsolution(m.internalModel)
        try
            m.redCosts = getreducedcosts(m.internalModel)
            m.linconstrDuals = getconstrduals(m.internalModel)
        catch
            warn("Dual solutions not available")
        end
        m.firstsolve = false
    end

    return stat

end

function solveMIP(m::Model)
    f, rowlb, rowub = prepProblemBounds(m)
    A = prepConstrMatrix(m)

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
    
    hasQuad = quadCheck(m, true)
   
    m.internalModel = model(m.solver)
    
    loadproblem!(m.internalModel, A, m.colLower, m.colUpper, f, rowlb, rowub, m.objSense)
    setvartype!(m.internalModel, vartype)

    if !all(isnan(m.colVal))
        try
            setwarmstart!(m.internalModel, m.colVal)
        catch
            Base.warn_once("MIP solver does not appear to support warm start solution.")
        end
    end

    if hasQuad
        addQuadratics(m)
    end
    registercallbacks(m)

    optimize!(m.internalModel)
    stat = status(m.internalModel)

    if stat != :Optimal
        println("Warning: MIP not solved to optimality, status: ", stat)
    end
    # It's possible that we have a feasible solution if we're not optimal
    # TODO: Test this behavior on various solvers
    try
        # store solution values in model
        m.objVal = getobjval(m.internalModel)
        m.objVal += m.obj.aff.constant
        m.colVal = getsolution(m.internalModel)
    end

    return stat

end


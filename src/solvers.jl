function solve(m::Model;IpoptOptions::Dict=Dict(),load_model_only=false, suppress_warnings=false)
    load_model_only == true && warn("load_model_only keyword is deprecated; use the buildInternalModel function instead")
    if m.nlpdata != nothing
        if length(IpoptOptions) > 0
            error("Specifying options by using IpoptOptions is no longer supported. Use \"m = Model(solver=IpoptSolver(option1=value1,option2=value2,...)\" instead, after loading the Ipopt package.")
        end
        if isa(m.solver,UnsetSolver)
            m.solver = MathProgBase.defaultNLPsolver
        end
        s = solvenlp(m, suppress_warnings=suppress_warnings)
        return s
    end
    # Analyze model to see if any integers
    anyInts = (length(m.sosconstr) > 0)
    if !anyInts
        for j = 1:m.numCols
            if m.colCat[j] != :Cont
                anyInts = true
                break
            end
        end
    end

    if isa(m.solver,UnsetSolver) &&
      (length(m.obj.qvars1) > 0 || length(m.quadconstr) > 0)
        m.solver = MathProgBase.defaultQPsolver
    end
    if anyInts
        if isa(m.solver,UnsetSolver)
            m.solver = MathProgBase.defaultMIPsolver
            s = solveMIP(m; load_model_only=load_model_only, suppress_warnings=suppress_warnings)
            # Clear solver in case we change problem types
            m.solver = UnsetSolver()
            m.internalModelLoaded = false
            return s
        else
            solveMIP(m; load_model_only=load_model_only, suppress_warnings=suppress_warnings)
        end
    else
        if isa(m.solver,UnsetSolver)
            m.solver = MathProgBase.defaultLPsolver
            s = solveLP(m, load_model_only=load_model_only, suppress_warnings=suppress_warnings)
            m.solver = UnsetSolver()
            return s
        else
            solveLP(m; load_model_only=load_model_only, suppress_warnings=suppress_warnings)
        end
    end
end

function addQuadratics(m::Model)

    if length(m.obj.qvars1) != 0
        assert_isfinite(m.obj)
        verify_ownership(m, m.obj.qvars1)
        verify_ownership(m, m.obj.qvars2)
        # Check for solver support for quadratic objectives happens in MPB
        MathProgBase.setquadobjterms!(m.internalModel, Cint[v.col for v in m.obj.qvars1], Cint[v.col for v in m.obj.qvars2], m.obj.qcoeffs)
    end

    # Add quadratic constraint to solver
    for k in 1:length(m.quadconstr)
        qconstr = m.quadconstr[k]
        if !((s = string(qconstr.sense)[1]) in ['<', '>', '='])
            error("Invalid sense for quadratic constraint")
        end
        terms::QuadExpr = qconstr.terms
        assert_isfinite(terms)
        for ind in 1:length(terms.qvars1)
            if (terms.qvars1[ind].m != m) || (terms.qvars2[ind].m != m)
                error("Variable not owned by model present in constraints")
            end
        end
        affidx = Cint[v.col for v in qconstr.terms.aff.vars]
        var1idx = Cint[v.col for v in qconstr.terms.qvars1]
        var2idx = Cint[v.col for v in qconstr.terms.qvars2]
        if applicable(MathProgBase.addquadconstr!, m.internalModel, affidx, qconstr.terms.aff.coeffs, var1idx, var2idx, qconstr.terms.qcoeffs, s, -qconstr.terms.aff.constant) 
            MathProgBase.addquadconstr!(m.internalModel, affidx, qconstr.terms.aff.coeffs, var1idx, var2idx, qconstr.terms.qcoeffs, s, -qconstr.terms.aff.constant)
        else
            error("Solver does not support quadratic constraints")
        end
    end
end

function addSOS(m::Model)
    for i in 1:length(m.sosconstr)
        sos = m.sosconstr[i]
        indices = Int[v.col for v in sos.terms]
        if sos.sostype == :SOS1
            if applicable(MathProgBase.addsos1!, m.internalModel, indices, sos.weights)
                MathProgBase.addsos1!(m.internalModel, indices, sos.weights)
            else
                error("Solver does not support SOS constraints")
            end
        elseif sos.sostype == :SOS2
            if applicable(MathProgBase.addsos2!, m.internalModel, indices, sos.weights)
                MathProgBase.addsos2!(m.internalModel, indices, sos.weights)
            else
                error("Solver does not support SOS constraints")
            end
        end
    end
end

# prepare objective, constraint matrix, and row bounds
function prepProblemBounds(m::Model)

    objaff::AffExpr = m.obj.aff
    assert_isfinite(objaff)
    verify_ownership(m, objaff.vars)
        
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
    linconstr = m.linconstr::Vector{LinearConstraint}
    numRows = length(linconstr)
    rowptr = Array(Int,numRows+1)
    nnz = 0
    for c in 1:numRows
        nnz += length(linconstr[c].terms.coeffs)
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
        assert_isfinite(linconstr[c].terms)
        coeffs = linconstr[c].terms.coeffs
        vars = linconstr[c].terms.vars
        # collect duplicates
        for ind in 1:length(coeffs)
            if !is(vars[ind].m, m)
                error("Variable not owned by model present in constraints")
            end
            addelt!(tmprow,vars[ind].col,coeffs[ind])
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
    # Note that rowmat doesn't have sorted indices, so technically doesn't
    # follow SparseMatrixCSC format. But it's safe to take the transpose.
    A = rowmat'
end

function solveLP(m::Model; load_model_only=false, suppress_warnings=false)
    f, rowlb, rowub = prepProblemBounds(m)  

    # Ready to solve
    noQuads = (length(m.quadconstr) == 0) && (length(m.obj.qvars1) == 0)
    if m.internalModelLoaded
        if applicable(MathProgBase.setvarLB!, m.internalModel, m.colLower) &&
           applicable(MathProgBase.setvarUB!, m.internalModel, m.colUpper) &&
           applicable(MathProgBase.setconstrLB!, m.internalModel, rowlb) &&
           applicable(MathProgBase.setconstrUB!, m.internalModel, rowub) &&
           applicable(MathProgBase.setobj!, m.internalModel, f) &&
           applicable(MathProgBase.setsense!, m.internalModel, m.objSense) &&
           applicable(MathProgBase.setvartype!, m.internalModel, [:Cont])            
            MathProgBase.setvarLB!(m.internalModel, m.colLower)
            MathProgBase.setvarUB!(m.internalModel, m.colUpper)
            MathProgBase.setconstrLB!(m.internalModel, rowlb)
            MathProgBase.setconstrUB!(m.internalModel, rowub)
            MathProgBase.setobj!(m.internalModel, f)
            MathProgBase.setsense!(m.internalModel, m.objSense)
            MathProgBase.setvartype!(m.internalModel, fill(:Cont,m.numCols))
        else
            !suppress_warnings && Base.warn_once("Solver does not appear to support hot-starts. Problem will be solved from scratch.")
            m.internalModelLoaded = false
        end
    end
    if !m.internalModelLoaded
        A = prepConstrMatrix(m)
        m.internalModel = MathProgBase.model(m.solver)
        MathProgBase.loadproblem!(m.internalModel, A, m.colLower, m.colUpper, f, rowlb, rowub, m.objSense)
        addQuadratics(m)
        m.internalModelLoaded = true
    end 

    if !load_model_only
        MathProgBase.optimize!(m.internalModel)
        stat = MathProgBase.status(m.internalModel)
    else
        stat = :NotSolved
    end

    if stat == :NotSolved
        # do nothing
    elseif stat != :Optimal
        !suppress_warnings && warn("Not solved to optimality, status: $stat")
        m.colVal = fill(NaN, m.numCols)
        m.objVal = NaN
        if stat == :Infeasible
            if noQuads && applicable(MathProgBase.getinfeasibilityray, m.internalModel)
                m.linconstrDuals = MathProgBase.getinfeasibilityray(m.internalModel)
            else
                noQuads && !suppress_warnings && warn("Infeasibility ray (Farkas proof) not available")
                m.linconstrDuals = fill(NaN, length(m.linconstr))
            end
        elseif stat == :Unbounded
            if noQuads && applicable(MathProgBase.getunboundedray, m.internalModel)
                m.colVal = MathProgBase.getunboundedray(m.internalModel)
            else
                noQuads && !suppress_warnings && warn("Unbounded ray not available")
                m.colVal = fill(NaN, m.numCols)
            end
        else
            try # guess try/catch is necessary because we're not sure what return status we have
                m.colVal = MathProgBase.getsolution(m.internalModel)
            catch
                m.colVal = fill(NaN, m.numCols)
            end
            try
                # store solution values in model
                m.objVal = MathProgBase.getobjval(m.internalModel)
                m.objVal += m.obj.aff.constant
            catch
                m.objVal = NaN
            end
        end
    else
        # store solution values in model
        m.objVal = MathProgBase.getobjval(m.internalModel)
        m.objVal += m.obj.aff.constant
        m.colVal = MathProgBase.getsolution(m.internalModel)
        if noQuads && applicable(MathProgBase.getreducedcosts, m.internalModel) &&
                      applicable(MathProgBase.getconstrduals,  m.internalModel)
            m.redCosts = MathProgBase.getreducedcosts(m.internalModel)
            m.linconstrDuals = MathProgBase.getconstrduals(m.internalModel)
        else
            noQuads && !suppress_warnings && warn("Dual solutions not available")
            m.redCosts = fill(NaN, length(m.linconstr))
            m.linconstrDuals = fill(NaN, length(m.linconstr))
        end
    end

    return stat
end

function solveMIP(m::Model; load_model_only=false, suppress_warnings=false)
    f, rowlb, rowub = prepProblemBounds(m)
    A = prepConstrMatrix(m)

    # Ready to solve

    
    if m.internalModelLoaded
        if applicable(MathProgBase.setvarLB!, m.internalModel, m.colLower) &&
           applicable(MathProgBase.setvarUB!, m.internalModel, m.colUpper) &&
           applicable(MathProgBase.setconstrLB!, m.internalModel, rowlb) &&
           applicable(MathProgBase.setconstrUB!, m.internalModel, rowub) &&
           applicable(MathProgBase.setobj!, m.internalModel, f) &&
           applicable(MathProgBase.setsense!, m.internalModel, m.objSense) &&
           applicable(MathProgBase.setvartype!, m.internalModel, m.colCat)
            MathProgBase.setvarLB!(m.internalModel, m.colLower)
            MathProgBase.setvarUB!(m.internalModel, m.colUpper)
            MathProgBase.setconstrLB!(m.internalModel, rowlb)
            MathProgBase.setconstrUB!(m.internalModel, rowub)
            MathProgBase.setobj!(m.internalModel, f)
            MathProgBase.setsense!(m.internalModel, m.objSense)
            MathProgBase.setvartype!(m.internalModel, m.colCat)
        else
            !suppress_warnings && Base.warn_once("Solver does not appear to support hot-starts. Problem will be solved from scratch.")
            m.internalModelLoaded = false
        end
    end
    if !m.internalModelLoaded
        m.internalModel = MathProgBase.model(m.solver)
        
        MathProgBase.loadproblem!(m.internalModel, A, m.colLower, m.colUpper, f, rowlb, rowub, m.objSense)
        if applicable(MathProgBase.setvartype!, m.internalModel, m.colCat)
            MathProgBase.setvartype!(m.internalModel, m.colCat)
        else
            error("Solver does not support discrete variables")
        end

        addSOS(m)

        addQuadratics(m)
        registercallbacks(m)

        m.internalModelLoaded = true
    end

    if !all(isnan(m.colVal))
        if applicable(MathProgBase.setwarmstart!, m.internalModel, m.colVal)
            MathProgBase.setwarmstart!(m.internalModel, m.colVal)
        else
            !suppress_warnings && Base.warn_once("Solver does not appear to support providing initial feasible solutions.")
        end
    end

    if !load_model_only
        MathProgBase.optimize!(m.internalModel)
        stat = MathProgBase.status(m.internalModel)
    else
        stat = :NotSolved
    end

    if stat == :NotSolved
        # do nothing
    else
        if stat != :Optimal
            !suppress_warnings && warn("Not solved to optimality, status: ", string(stat))
        end
        # It's possible that we have a feasible solution if we're not optimal
        # TODO: Test this behavior on various solvers
        try
            # store solution values in model
            m.objVal = MathProgBase.getobjval(m.internalModel)
            m.objVal += m.obj.aff.constant
        catch
            m.objVal = NaN
        end
        try
            m.colVal = MathProgBase.getsolution(m.internalModel)
        catch
            m.colVal = fill(NaN, m.numCols)
        end
    end

    return stat
end

function buildInternalModel(m::Model)
    m.nlpdata == nothing || error("buildInternalModel not supported for nonlinear problems")

    anyInts = false
    for j = 1:m.numCols
        if m.colCat[j] != :Cont
            anyInts = true
            break
        end
    end

    if isa(m.solver,UnsetSolver) &&
      (length(m.obj.qvars1) > 0 || length(m.quadconstr) > 0)
        m.solver = MathProgBase.defaultQPsolver
    end
    if anyInts
        if isa(m.solver,UnsetSolver)
            m.solver = MathProgBase.defaultMIPsolver
        end
    else
        if isa(m.solver,UnsetSolver)
            m.solver = MathProgBase.defaultLPsolver
        end
    end
    m.internalModel = MathProgBase.model(m.solver)

    f, rowlb, rowub = prepProblemBounds(m)
    A = prepConstrMatrix(m)
    MathProgBase.loadproblem!(m.internalModel, A, m.colLower, m.colUpper, f, rowlb, rowub, m.objSense)
    addQuadratics(m)

    if anyInts # do MIP stuff
        MathProgBase.setvartype!(m.internalModel, m.colCat)
        addSOS(m)
        registercallbacks(m)
        if !all(isnan(m.colVal))
            if applicable(MathProgBase.setwarmstart!, m.internalModel, m.colVal)
                MathProgBase.setwarmstart!(m.internalModel, m.colVal)
            else
                Base.warn_once("Solver does not appear to support providing initial feasible solutions.")
            end
        end
    end
    m.internalModelLoaded = true
    nothing
end

# returns (unsorted) column indices and coefficient terms for merged vector
# assume that v is zero'd
function merge_duplicates{CoefType,IntType<:Integer}(::Type{IntType},aff::GenericAffExpr{CoefType,Variable}, v::IndexedVector{CoefType}, m::Model)
    resize!(v, m.numCols)
    for ind in 1:length(aff.coeffs)
        var = aff.vars[ind]
        is(var.m, m) || error("Variable does not belong to this model")
        addelt!(v, aff.vars[ind].col, aff.coeffs[ind])
    end
    indices = Array(IntType,v.nnz)
    coeffs = Array(CoefType,v.nnz)
    for i in 1:v.nnz
        idx = v.nzidx[i]
        indices[i] = idx
        coeffs[i] = v.elts[idx]
    end
    empty!(v)

    return indices, coeffs

end

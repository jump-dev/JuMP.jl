function solve(m::Model;IpoptOptions::Dict=Dict(),load_model_only=false, suppress_warnings=false)
    if m.nlpdata != nothing
        return solveIpopt(m, options=IpoptOptions, suppress_warnings=suppress_warnings)
    end
    # Analyze model to see if any integers
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
        setquadobjterms!(m.internalModel, Cint[v.col for v in m.obj.qvars1], Cint[v.col for v in m.obj.qvars2], m.obj.qcoeffs)
    end

    # Add quadratic constraint to solver
    for k in 1:length(m.quadconstr)
        qconstr = m.quadconstr[k]
        if !((s = string(qconstr.sense)[1]) in ['<', '>', '='])
            error("Invalid sense for quadratic constraint")
        end
        terms = qconstr.terms
        for ind in 1:length(terms.qvars1)
            if (terms.qvars1[ind].m != m) || (terms.qvars2[ind].m != m)
                error("Variable not owned by model present in constraints")
            end
        end
        addquadconstr!(m.internalModel, Cint[v.col for v in qconstr.terms.aff.vars], qconstr.terms.aff.coeffs, Cint[v.col for v in qconstr.terms.qvars1], Cint[v.col for v in qconstr.terms.qvars2], qconstr.terms.qcoeffs, s, -qconstr.terms.aff.constant)
    end
end

function addSOS(m::Model)
    try
        for i in 1:length(m.sosconstr)
            sos = m.sosconstr[i]
            indices = Int[v.col for v in sos.terms]
            if sos.sostype == :SOS1
                addsos1!(m.internalModel, indices, sos.weights)
            elseif sos.sostype == :SOS2
                addsos2!(m.internalModel, indices, sos.weights)
            end
        end
    catch
        for i in 1:length(m.sosconstr)
            sos = m.sosconstr[i]
            indices = Int[v.col for v in sos.terms]
            nvars = length(indices)
            if sos.sostype == :SOS1
                Base.warn_once("Current solver does not support SOS1 constraints, adding manually")
                addconstr!(m.internalModel, indices, ones(nvars), 0., 1.)
            elseif sos.sostype == :SOS2
                error("Current solver does not support SOS2 constraints")
            end
        end
    end
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
    A = rowmat'
end

function solveLP(m::Model; load_model_only=false, suppress_warnings=false)
    f, rowlb, rowub = prepProblemBounds(m)  

    # Ready to solve

    if m.internalModelLoaded
        try
            setvarLB!(m.internalModel, m.colLower)
            setvarUB!(m.internalModel, m.colUpper)
            setconstrLB!(m.internalModel, rowlb)
            setconstrUB!(m.internalModel, rowub)
            setobj!(m.internalModel, f)
            setsense!(m.internalModel, m.objSense)
        catch
            !suppress_warnings && Base.warn_once("Solver does not appear to support hot-starts. Problem will be solved from scratch.")
            m.internalModelLoaded = false
        end
        all_cont = true
        try # this fails for LPs for some unfathomable reason...but if it's an LP, we're good anyway
            all_cont = mapreduce(x->isequal([:Cont],x), &, MathProgBase.getvartype(m.internalModel))
        end
        if !all_cont
            setvartype!(m.internalModel, fill(:Cont,m.numCols))
        end
    end
    if !m.internalModelLoaded
        A = prepConstrMatrix(m)
        m.internalModel = model(m.solver)
        loadproblem!(m.internalModel, A, m.colLower, m.colUpper, f, rowlb, rowub, m.objSense)
        addQuadratics(m)
        m.internalModelLoaded = true
    end 

    if !load_model_only
        optimize!(m.internalModel)
        stat = status(m.internalModel)
    else
        stat = :NotSolved
    end

    if stat == :NotSolved
        # do nothing
    elseif stat != :Optimal
        !suppress_warnings && warn("Not solved to optimality, status: $stat")
        if stat == :Infeasible
            try
                m.linconstrDuals = getinfeasibilityray(m.internalModel)
            catch
                !suppress_warnings && warn("Infeasibility ray (Farkas proof) not available")
            end
        elseif stat == :Unbounded
            try
                m.colVal = getunboundedray(m.internalModel)
            catch
                !suppress_warnings && warn("Unbounded ray not available")
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
            !suppress_warnings && warn("Dual solutions not available")
        end
    end

    return stat
end

function solveMIP(m::Model; load_model_only=false, suppress_warnings=false)
    f, rowlb, rowub = prepProblemBounds(m)
    A = prepConstrMatrix(m)

    # Ready to solve
    if m.internalModelLoaded
        try
            setvarLB!(m.internalModel, m.colLower)
            setvarUB!(m.internalModel, m.colUpper)
            setconstrLB!(m.internalModel, rowlb)
            setconstrUB!(m.internalModel, rowub)
            setobj!(m.internalModel, f)
            setsense!(m.internalModel, m.objSense)
            setvartype!(m.internalModel, m.colCat)
        catch
            m.internalModelLoaded = false
        end
    end
    if !m.internalModelLoaded
        m.internalModel = model(m.solver)
        
        loadproblem!(m.internalModel, A, m.colLower, m.colUpper, f, rowlb, rowub, m.objSense)
        setvartype!(m.internalModel, m.colCat)

        addSOS(m)

        addQuadratics(m)
        registercallbacks(m)

        m.internalModelLoaded = true
    end

    if !all(isnan(m.colVal))
        try
            setwarmstart!(m.internalModel, m.colVal)
        catch
            !suppress_warnings && Base.warn_once("Solver does not appear to support providing initial feasible solutions.")
        end
    end

    if !load_model_only
        optimize!(m.internalModel)
        stat = status(m.internalModel)
    else
        stat = :NotSolved
    end

    if stat == :NotSolved
        # do nothing
    else
        # It's possible that we have a feasible solution if we're not optimal
        # TODO: Test this behavior on various solvers
        try
            # store solution values in model
            m.objVal = getobjval(m.internalModel)
            m.objVal += m.obj.aff.constant
            m.colVal = getsolution(m.internalModel)
        end
    end
    if stat != :Optimal
        !suppress_warnings && warn("Not solved to optimality, status: $stat")
    end

    return stat
end

# currently used only in callbacks
# returns (unsorted) column indices and coefficient terms for merged vector
# assume that v is zero'd and has the right size (total number of variables in the model)
function merge_duplicates{CoefType,IntType<:Integer}(::Type{IntType},aff::GenericAffExpr{CoefType,Variable}, v::IndexedVector{CoefType}, m::Model)
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

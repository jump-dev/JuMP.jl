#  Copyright 2015, Iain Dunning, Joey Huchette, Miles Lubin, and contributors
#  This Source Code Form is subject to the terms of the Mozilla Public
#  License, v. 2.0. If a copy of the MPL was not distributed with this
#  file, You can obtain one at http://mozilla.org/MPL/2.0/.

immutable ProblemTraits
    int::Bool
    lin::Bool
    qp ::Bool
    qc ::Bool
    nlp::Bool
    sdp::Bool
    sos::Bool
end

function problemclass(m::Model)
    int = any(c-> !(c == :Cont || c == :Fixed), m.colCat)
    qp = !isempty(m.obj.qvars1)
    qc = !isempty(m.quadconstr)
    nlp = m.nlpdata != nothing
    sdp = !isempty(m.sdpconstr)
    sos = !isempty(m.sosconstr)
    ProblemTraits(int, !(qp|qc|nlp|sdp|sos), qp, qc, nlp, sdp, sos)
end

function solve(m::Model; suppress_warnings=false, ignore_solve_hook=(m.solvehook==nothing), kwargs...)

    ignore_solve_hook || return m.solvehook(m; suppress_warnings=suppress_warnings, kwargs...)

    # Clear warning counters
    m.getvalue_counter = 0
    m.operator_counter = 0

    # Unfortunately we still have to take a different codepath for nlp
    if m.nlpdata != nothing
        if isa(m.solver,UnsetSolver)
            m.solver = MathProgBase.defaultNLPsolver
        end
        return solvenlp(m, suppress_warnings=suppress_warnings)
    end

    unset = m.solver == UnsetSolver()
    traits = problemclass(m)
    buildInternalModel(m, traits, suppress_warnings=suppress_warnings)

    MathProgBase.optimize!(m.internalModel)
    stat = MathProgBase.status(m.internalModel)

    if stat == :Optimal
        if !(traits.int | traits.sos | traits.sdp)
            m.redCosts = try
                MathProgBase.getreducedcosts(m.internalModel)
            catch
                fill(NaN, length(m.colVal))
            end

            m.linconstrDuals = try
                MathProgBase.getconstrduals(m.internalModel)
            catch
                fill(NaN, length(m.linconstr))
            end
        end
    else
        suppress_warnings || warn("Not solved to optimality, status: $stat")

        if stat == :Infeasible
            m.linconstrDuals = try
                MathProgBase.getinfeasibilityray(m.internalModel)
            catch
                suppress_warnings || warn("Infeasibility ray (Farkas proof) not available")
                fill(NaN, length(m.linconstr))
            end
        elseif stat == :Unbounded
            try
                m.colVal = MathProgBase.getunboundedray(m.internalModel)
            catch
                suppress_warnings || warn("Unbounded ray not available")
            end
        end
    end

    m.objVal = NaN
    m.colVal = fill(NaN, m.numCols)
    if !(stat == :Infeasible || stat == :Unbounded)
        try
            objVal = MathProgBase.getobjval(m.internalModel) + m.obj.aff.constant
            colVal = MathProgBase.getsolution(m.internalModel)
            m.objVal = objVal # Don't corrupt the answers if one of the above two calls fails
            m.colVal = colVal
        end
    end

    if traits.sdp && m.objSense == :Max
        m.objVal *= -1
    end

    if unset
        m.solver = UnsetSolver()
        if traits.int
            m.internalModelLoaded = false
        end
    end

    stat
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
    const sensemap = @compat Dict(:(<=) => '<', :(>=) => '>', :(==) => '=')
    for k in 1:length(m.quadconstr)
        qconstr = m.quadconstr[k]::QuadConstraint
        if !haskey(sensemap, qconstr.sense)
            error("Invalid sense for quadratic constraint")
        end
        s = sensemap[qconstr.sense]

        terms::QuadExpr = qconstr.terms
        assert_isfinite(terms)
        for ind in 1:length(terms.qvars1)
            if (terms.qvars1[ind].m != m) || (terms.qvars2[ind].m != m)
                error("Variable not owned by model present in constraints")
            end
        end
        affidx  = Cint[v.col for v in qconstr.terms.aff.vars]
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
    linconstr = m.linconstr::Vector{LinearConstraint}
    numRows = length(linconstr)
    rowlb = fill(-Inf, numRows)
    rowub = fill(+Inf, numRows)
    for c in 1:numRows
        rowlb[c] = linconstr[c].lb
        rowub[c] = linconstr[c].ub
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

function vartypes_without_fixed(m::Model)
    colCats = copy(m.colCat)
    for i in 1:length(colCats)
        if colCats[i] == :Fixed
            @assert m.colLower[i] == m.colUpper[i]
            colCats[i] = :Cont
        end
    end
    return colCats
end

function collect_expr!(m, tmprow, terms::AffExpr)
    empty!(tmprow)
    assert_isfinite(terms)
    coeffs = terms.coeffs
    vars = terms.vars
    # collect duplicates
    for ind in 1:length(coeffs)
        if !is(vars[ind].m, m)
            error("Variable not owned by model present in constraints")
        end
        addelt!(tmprow,vars[ind].col, coeffs[ind])
    end
    tmprow
end

function conicconstraintdata(m::Model)
    var_cones = Any[]
    con_cones = Any[]
    nnz = 0

    # find starting column indices for sdp matrices
    sdp_start, sdp_end = Int[], Int[]
    numSDPRows = 0
    numSymRows = 0
    for c in m.sdpconstr
        n = size(c.terms,1)
        @assert n == size(c.terms,2)
        @assert ndims(c.terms) == 2
        if isa(c.terms,OneIndexedArray)
            frst = c.terms[1,1].col
            last = c.terms[end,end].col
            push!(sdp_start, frst)
            push!(sdp_end, last)
            push!(var_cones, (:SDP,frst:last))
        else
            numSDPRows += convert(Int, n*(n+1)/2)
            for i in 1:n, j in i:n
                nnz += length(c.terms[i,j].coeffs)
            end
        end
        if !issym(c.terms)
            # symmetry constraints
            numSymRows += convert(Int, n*(n-1)/2)
        end
    end

    soc_cones  = Any[]
    rsoc_cones = Any[]
    numQuadRows = 0
    for qconstr in m.quadconstr
        q = copy(qconstr.terms)
        if qconstr.sense == :(>=)
            q *= -1
        end
        if !(isempty(q.aff.vars) && q.aff.constant == 0)
            error("Quadratic constraint $qconstr must be in second-order cone form")
        end
        n_pos_on_diag = 0
        off_diag_idx  = 0
        neg_diag_idx  = 0
        n = length(q.qvars1)
        for i in 1:n
            if q.qvars1[i].col == q.qvars2[i].col
                if q.qcoeffs[i] == 1
                    n_pos_on_diag += 1
                elseif q.qcoeffs[i] == -1
                    neg_diag_idx == off_diag_idx == 0 || error("Invalid SOC constraint $qconstr")
                    neg_diag_idx = i
                end
            else
                if q.qcoeffs[i] == -1
                    neg_diag_idx == off_diag_idx == 0 || error("Invalid rotated SOC constraint $qconstr")
                    off_diag_idx = i
                end
            end
        end
        cone = Array(Int, n)
        if n_pos_on_diag == n-1 && neg_diag_idx > 0
            cone[1] = q.qvars1[neg_diag_idx].col
            for i in 1:(neg_diag_idx-1); cone[i+1] = q.qvars1[i].col; end
            for i in (neg_diag_idx+1):n; cone[i]   = q.qvars1[i].col; end
            push!(soc_cones, cone)
        elseif n_pos_on_diag == n-1 && off_diag_idx > 0
            cone[1] = q.qvars1[off_diag_idx].col
            cone[2] = q.qvars2[off_diag_idx].col
            for i in 1:(off_diag_idx-1); cone[i+2] = q.qvars1[i].col; end
            for i in (off_diag_idx+1):n; cone[i+1] = q.qvars1[i].col; end
            push!(rsoc_cones, cone)
        else
            error("Quadratic constraint $qconstr is not conic representable")
        end
        numQuadRows += length(cone)
    end

    # Create sparse A matrix
    # First we build it row-wise, then use the efficient transpose
    # Theory is, this is faster than us trying to do it ourselves
    # Intialize storage
    linconstr = m.linconstr::Vector{LinearConstraint}
    numLinRows = length(linconstr)
    numBounds = 0
    nonNeg  = Int[]
    nonPos  = Int[]
    free    = Int[]
    zeroVar = Int[]
    in_sdp = false
    for i in 1:m.numCols
        lb, ub = m.colLower[i], m.colUpper[i]
        if i in sdp_start
            in_sdp = true
            @assert lb == -Inf && ub == Inf
        end

        if !(lb == 0 || lb == -Inf)
            numBounds += 1
        end
        if !(ub == 0 || ub == Inf)
            numBounds += 1
        end
        if lb == 0 && ub == 0
            push!(zeroVar, i)
        elseif lb == 0
            push!(nonNeg, i)
        elseif ub == 0
            push!(nonPos, i)
        elseif in_sdp
            # do nothing
        else
            push!(free, i)
        end
        if in_sdp
            @assert lb == -Inf && ub == Inf
            if i in sdp_end
                in_sdp = false
            end
        end
    end

    if !isempty(zeroVar)
        push!(var_cones, (:Zero,zeroVar))
    end
    if !isempty(nonNeg)
        push!(var_cones, (:NonNeg,nonNeg))
    end
    if !isempty(nonPos)
        push!(var_cones, (:NonPos,nonPos))
    end
    if !isempty(free)
        push!(var_cones, (:Free,free))
    end

    nnz += numBounds
    for c in 1:numLinRows
        nnz += length(linconstr[c].terms.coeffs)
    end

    numRows = numLinRows + numBounds + numQuadRows + numSDPRows + numSymRows

    b = Array(Float64, numRows)

    I = Int[]
    J = Int[]
    V = Float64[]
    @compat sizehint!(I, nnz)
    @compat sizehint!(J, nnz)
    @compat sizehint!(V, nnz)

    # Fill it up
    nnz = 0
    tmprow = IndexedVector(Float64,m.numCols)
    tmpelts = tmprow.elts
    tmpnzidx = tmprow.nzidx
    nonneg_rows = Int[]
    nonpos_rows = Int[]
    eq_rows     = Int[]
    for c in 1:numLinRows
        if linconstr[c].lb == -Inf
            b[c] = linconstr[c].ub
            push!(nonneg_rows, c)
        elseif linconstr[c].ub == Inf
            b[c] = linconstr[c].lb
            push!(nonpos_rows, c)
        elseif linconstr[c].lb == linconstr[c].ub
            b[c] = linconstr[c].lb
            push!(eq_rows, c)
        else
            error("We currently do not support ranged constraints with conic solvers")
        end

        assert_isfinite(linconstr[c].terms)
        coeffs = linconstr[c].terms.coeffs
        vars = linconstr[c].terms.vars
        # collect duplicates
        for ind in 1:length(coeffs)
            if !is(vars[ind].m, m)
                error("Variable not owned by model present in constraints")
            end
            addelt!(tmprow,vars[ind].col, coeffs[ind])
        end
        nnz = tmprow.nnz
        append!(I, fill(c, nnz))
        indices = tmpnzidx[1:nnz]
        append!(J, indices)
        append!(V, tmpelts[indices])
        empty!(tmprow)
    end

    c = numLinRows
    for idx in 1:m.numCols
        lb = m.colLower[idx]
        if !(lb == 0 || lb == -Inf)
            nnz += 1
            c   += 1
            push!(I, c)
            push!(J, idx)
            push!(V, 1.0)
            b[c] = lb
            push!(nonpos_rows, c)
        end
        ub = m.colUpper[idx]
        if !(ub == 0 || ub == Inf)
            c   += 1
            push!(I, c)
            push!(J, idx)
            push!(V, 1.0)
            b[c] = ub
            push!(nonneg_rows, c)
        end
    end

    if !isempty(nonneg_rows)
        push!(con_cones, (:NonNeg,nonneg_rows))
    end
    if !isempty(nonpos_rows)
        push!(con_cones, (:NonPos,nonpos_rows))
    end
    if !isempty(eq_rows)
        push!(con_cones, (:Zero,eq_rows))
    end
    @assert c == numLinRows + numBounds

    for cone in soc_cones
        n = length(cone)
        rng = (c+1):(c+n)
        append!(I, rng)
        append!(J, copy(cone))
        append!(V, [-1.0; ones(n-1)])
        push!(con_cones, (:SOC,rng))
        c += n
    end
    for cone in rsoc_cones
        n = length(cone)
        rng = (c+1):(c+n)
        append!(I, rng)
        append!(J, copy(cone))
        append!(V, [-1.0; -1.0; ones(n-2)])
        push!(con_cones, (:SOCRotated,rng))
        c += n
    end
    @assert c == numLinRows + numBounds + numQuadRows

    tmpelts = tmprow.elts
    tmpnzidx = tmprow.nzidx
    for con in m.sdpconstr
        if !isa(con.terms, OneIndexedArray)
            sdp_start = c + 1
            n = size(con.terms,1)
            for i in 1:n, j in i:n
                c += 1
                terms::AffExpr = con.terms[i,j]
                collect_expr!(m, tmprow, terms)
                nnz = tmprow.nnz
                indices = tmpnzidx[1:nnz]
                append!(I, fill(c, nnz))
                append!(J, indices)
                append!(V, -tmpelts[indices])
                b[c] = terms.constant
            end
            push!(con_cones, (:SDP, sdp_start:c))
        end
        if !issym(con.terms)
            sym_start = c + 1
            # add linear symmetry constraints
            for i in 1:n, j in 1:(i-1)
                c += 1
                collect_expr!(m, tmprow, con.terms[i,j] - con.terms[j,i])
                nnz = tmprow.nnz
                indices = tmpnzidx[1:nnz]
                append!(I, fill(c, nnz))
                append!(J, indices)
                append!(V, tmpelts[indices])
                b[c] = 0
            end
            push!(con_cones, (:Zero, sym_start:c))
        end
    end
    @assert c == numRows

    A = sparse(I, J, V, numRows, m.numCols)
    # @show full(A), b
    # @show var_cones, con_cones

    # TODO: uncomment these lines when they work with Mosek
    # supported = MathProgBase.supportedcones(m.internalModel)
    # @assert (:NonNeg in supported) && (:NonPos in supported) && (:Free in supported) && (:SDP in supported)
    A, b, var_cones, con_cones
end

function buildInternalModel(m::Model, traits=problemclass(m); suppress_warnings=false)
    traits.nlp && error("buildInternalModel not supported for nonlinear problems")

    # set default solver
    if isa(m.solver, UnsetSolver)
        m.solver = begin
            if traits.int | traits.sos
                MathProgBase.defaultMIPsolver
            elseif traits.sdp
                MathProgBase.defaultSDPsolver
            elseif traits.qp | traits.qc
                MathProgBase.defaultQPsolver
            else
                MathProgBase.defaultLPsolver
            end
        end
    end

    if traits.sdp # really should be conic here
        f,_,_ = prepProblemBounds(m)
        if m.objSense == :Max
            scale!(f, -1)
        end

        m.internalModel = MathProgBase.model(m.solver)

        A, b, var_cones, con_cones = conicconstraintdata(m)
        MathProgBase.loadconicproblem!(m.internalModel, f, A, b, con_cones, var_cones)
    else
        f, rowlb, rowub = prepProblemBounds(m)
        if m.internalModelLoaded
            if applicable(MathProgBase.setvarLB!, m.internalModel, m.colLower) &&
               applicable(MathProgBase.setvarUB!, m.internalModel, m.colUpper) &&
               applicable(MathProgBase.setconstrLB!, m.internalModel, rowlb) &&
               applicable(MathProgBase.setconstrUB!, m.internalModel, rowub) &&
               applicable(MathProgBase.setobj!, m.internalModel, f) &&
               applicable(MathProgBase.setsense!, m.internalModel, m.objSense)
                MathProgBase.setvarLB!(m.internalModel, m.colLower)
                MathProgBase.setvarUB!(m.internalModel, m.colUpper)
                MathProgBase.setconstrLB!(m.internalModel, rowlb)
                MathProgBase.setconstrUB!(m.internalModel, rowub)
                MathProgBase.setobj!(m.internalModel, f)
                MathProgBase.setsense!(m.internalModel, m.objSense)
            else
                suppress_warnings || Base.warn_once("Solver does not appear to support hot-starts. Model will be built from scratch.")
                m.internalModelLoaded = false
            end
        end
        if !m.internalModelLoaded
            m.internalModel = MathProgBase.model(m.solver)
            A = prepConstrMatrix(m)

            # if we have either:
            #   1) A solver that does not support the loadproblem! interface, or
            #   2) A QCP and a solver that does not support the addquadconstr! interface,
            # wrap everything in a ConicSolverWrapper
            if !applicable(MathProgBase.loadproblem!, m.internalModel, A, m.colLower, m.colUpper, f, rowlb, rowub, m.objSense) ||
                ( applicable(MathProgBase.supportedcones, m.solver) && # feel like this should have a && traits.qc as well...
                  !method_exists(MathProgBase.addquadconstr!, (typeof(m.internalModel), Vector{Int}, Vector{Float64}, Vector{Int}, Vector{Int}, Vector{Float64}, Char, Float64)) &&
                  :SOC in MathProgBase.supportedcones(m.solver) )

                m.internalModel = MathProgBase.model(MathProgBase.ConicSolverWrapper(m.solver))
            end

            MathProgBase.loadproblem!(m.internalModel, A, m.colLower, m.colUpper, f, rowlb, rowub, m.objSense)
            addQuadratics(m)
        end

        addSOS(m)
        registercallbacks(m)
    end

    if applicable(MathProgBase.setvartype!, m.internalModel, Symbol[])
        colCats = vartypes_without_fixed(m)
        MathProgBase.setvartype!(m.internalModel, colCats)
    elseif traits.int
        error("Solver does not support discrete variables")
    end

    # TODO: change this so you can warmstart continuous problems?
    if traits.int && !all(isnan(m.colVal))
        if applicable(MathProgBase.setwarmstart!, m.internalModel, m.colVal)
            MathProgBase.setwarmstart!(m.internalModel, m.colVal)
        else
            suppress_warnings || Base.warn_once("Solver does not appear to support providing initial feasible solutions.")
        end
    end

    if applicable(MathProgBase.updatemodel!, m.internalModel)
        MathProgBase.updatemodel!(m.internalModel)
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

#  Copyright 2015, Iain Dunning, Joey Huchette, Miles Lubin, and contributors
#  This Source Code Form is subject to the terms of the Mozilla Public
#  License, v. 2.0. If a copy of the MPL was not distributed with this
#  file, You can obtain one at http://mozilla.org/MPL/2.0/.
#############################################################################
# JuMP
# An algebraic modeling language for Julia
# See http://github.com/JuliaOpt/JuMP.jl
#############################################################################
# src/solvers.jl
# Handles conversion of the JuMP Model into a format that can be passed
# through the MathProgBase interface to solvers, and ongoing updating of
# that representation if supported by the solver.
#############################################################################

# Analyze a JuMP Model to determine its traits, and thus what solvers can
# be used to solve the problem
type ProblemTraits
    int::Bool  # has integer variables
    lin::Bool  # has only linear objectives and constraints
    qp ::Bool  # has a quadratic objective function
    qc ::Bool  # has a quadratic constraint
    nlp::Bool  # has general nonlinear objective or constraints
    soc::Bool  # has a second-order cone constraint
    sdp::Bool  # has an SDP constraint (or SDP variable bounds)
    sos::Bool  # has an SOS constraint
    conic::Bool  # has an SDP or SOC constraint
end
function ProblemTraits(m::Model)
    int = any(c-> !(c == :Cont || c == :Fixed), m.colCat)
    qp = !isempty(m.obj.qvars1)
    qc = !isempty(m.quadconstr)
    nlp = m.nlpdata !== nothing
    soc = !isempty(m.socconstr)
    # will need to change this when we add support for arbitrary variable cones
    sdp = !isempty(m.sdpconstr) || !isempty(m.varCones)
    sos = !isempty(m.sosconstr)
    ProblemTraits(int, !(qp|qc|nlp|soc|sdp|sos), qp, qc, nlp, soc, sdp, sos, soc|sdp)
end
function default_solver(traits::ProblemTraits)
    if traits.int || traits.sos
        MathProgBase.defaultMIPsolver
    elseif traits.sdp
        MathProgBase.defaultSDPsolver
    elseif traits.conic
        MathProgBase.defaultConicsolver
    elseif traits.qp || traits.qc
        MathProgBase.defaultQPsolver
    elseif traits.nlp
        MathProgBase.defaultNLPsolver
    else
        MathProgBase.defaultLPsolver
    end
end

function fillConicRedCosts(m::Model)
    bndidx = 0
    numlinconstr = length(m.linconstr)
    for i in 1:m.numCols
        lower = false
        upper = false
        lb, ub = m.colLower[i], m.colUpper[i]

        if lb != -Inf
            lower = true
            bndidx += 1
        end
        if ub != Inf
            upper = true
            bndidx += 1
        end

        if lower && !upper
            m.redCosts[i] = m.conicconstrDuals[numlinconstr + bndidx]
        elseif !lower && upper
            m.redCosts[i] = m.conicconstrDuals[numlinconstr + bndidx]
        elseif lower && upper
            m.redCosts[i] = m.conicconstrDuals[numlinconstr + bndidx]+m.conicconstrDuals[numlinconstr + bndidx-1]
        end
    end
end

function fillConicDuals(m::Model)

    numRows, numCols = length(m.linconstr), m.numCols

    numBndRows = getNumBndRows(m)
    numSOCRows = getNumSOCRows(m)
    m.conicconstrDuals = try
        MathProgBase.getdual(m.internalModel)
    catch
        fill(NaN, numRows+numBndRows+numSOCRows)
    end
    if m.conicconstrDuals[1] != NaN
        if m.objSense == :Min
            scale!(m.conicconstrDuals, -1)
        end
        m.linconstrDuals = m.conicconstrDuals[1:length(m.linconstr)]
        m.redCosts = zeros(numCols)
        if numBndRows > 0
            fillConicRedCosts(m)
        end
    end

end

function solve(m::Model; suppress_warnings=false,
                ignore_solve_hook=(m.solvehook===nothing),
                relaxation=false,
                kwargs...)
    # If the user or an extension has provided a solve hook, call
    # that instead of solving the model ourselves
    if !ignore_solve_hook
        return m.solvehook(m; suppress_warnings=suppress_warnings, kwargs...)::Symbol
    end

    isempty(kwargs) || error("Unrecognized keyword arguments: $(join([k[1] for k in kwargs], ", "))")

    # Clear warning counters
    m.getvalue_counter = 0
    m.operator_counter = 0

    # Remember if the solver was initially unset so we can restore
    # it to be unset later
    unset = m.solver == UnsetSolver()

    # Analyze the problems traits to determine what solvers we can use
    traits = ProblemTraits(m)

    # Build the MathProgBase model from the JuMP model
    buildInternalModel(m, traits, suppress_warnings=suppress_warnings, relaxation=relaxation)

    # If the model is a general nonlinear, use different logic in
    # nlp.jl to solve the problem
    traits.nlp && return solvenlp(m, traits, suppress_warnings=suppress_warnings)

    # Solve the problem
    MathProgBase.optimize!(m.internalModel)
    stat::Symbol = MathProgBase.status(m.internalModel)

    # Extract solution from the solver
    numRows, numCols = length(m.linconstr), m.numCols
    m.objVal = NaN
    m.colVal = fill(NaN, numCols)
    m.linconstrDuals = Array(Float64, 0)

    discrete = !relaxation && (traits.int || traits.sos)
    if stat == :Optimal
        # If we think dual information might be available, try to get it
        # If not, return an array of the correct length
        if !discrete && !traits.conic
            m.redCosts = try
                MathProgBase.getreducedcosts(m.internalModel)[1:numCols]
            catch
                fill(NaN, numCols)
            end

            m.linconstrDuals = try
                MathProgBase.getconstrduals(m.internalModel)[1:numRows]
            catch
                fill(NaN, numRows)
            end
        end
        # conic duals (currently, SOC only)
        if !discrete && traits.soc && !traits.qp && !traits.qc && !traits.sdp
            fillConicDuals(m)
        end
    else
        # Problem was not solved to optimality, attempt to extract useful
        # information anyway
        suppress_warnings || warn("Not solved to optimality, status: $stat")
        # Some solvers provide infeasibility rays (dual) or unbounded
        # rays (primal) for linear problems. Store these as the solution
        # if the exist.
        if traits.lin
            if stat == :Infeasible
                m.linconstrDuals = try
                    infray = MathProgBase.getinfeasibilityray(m.internalModel)
                    @assert length(infray) == numRows
                    infray
                catch
                    suppress_warnings || warn("Infeasibility ray (Farkas proof) not available")
                    fill(NaN, numRows)
                end
            elseif stat == :Unbounded
                m.colVal = try
                    unbdray = MathProgBase.getunboundedray(m.internalModel)
                    @assert length(unbdray) == numCols
                    unbdray
                catch
                    suppress_warnings || warn("Unbounded ray not available")
                    fill(NaN, numCols)
                end
            end
        end
        # conic duals (currently, SOC only)
        if !discrete && traits.soc && !traits.qp && !traits.qc && !traits.sdp
            if stat == :Infeasible
                fillConicDuals(m)
            end
        end
    end

    # If the problem was solved, or if it terminated prematurely, try
    # to extract a solution anyway. This commonly occurs when a time
    # limit or tolerance is set (:UserLimit)
    if !(stat == :Infeasible || stat == :Unbounded)
        try
            objVal = MathProgBase.getobjval(m.internalModel) + m.obj.aff.constant
            colVal = MathProgBase.getsolution(m.internalModel)[1:numCols]
            # Don't corrupt the answers if one of the above two calls fails
            m.objVal = objVal
            m.colVal = colVal
        end
    end

    # The MathProgBase interface defines a conic problem to always be
    # a minimization problem, so we need to flip the objective before
    # reporting it to the user
    if traits.conic && m.objSense == :Max
        m.objVal *= -1
    end

    # If the solver was initially not set, we will restore this status
    # and drop the internal MPB model. This is important for the case
    # where the solver used changes between solves because the user
    # has changed the problem class (e.g. LP to MILP)
    if unset
        m.solver = UnsetSolver()
        if traits.int
            m.internalModelLoaded = false
        end
    end

    # don't keep relaxed model in memory
    relaxation && (m.internalModelLoaded = false)

    # Return the solve status
    stat
end

function isquadsoc(m::Model)
    # check if all quadratic constraints are actually conic
    all_conic = true
    for qconstr in m.quadconstr
        q = copy(qconstr.terms)
        if qconstr.sense == :(>=)
            q *= -1
        end
        if !(all(t->t==0, q.aff.coeffs) && q.aff.constant == 0)
            all_conic = false
            break
        end
        n_pos_on_diag = 0
        off_diag_idx  = 0
        neg_diag_idx  = 0
        n = length(q.qvars1)
        nz = 0
        for i in 1:n
            q.qcoeffs[i] == 0 && continue
            nz += 1
            if q.qvars1[i].col == q.qvars2[i].col
                if q.qcoeffs[i] == 1
                    n_pos_on_diag += 1
                elseif q.qcoeffs[i] == -1
                    if !(neg_diag_idx == off_diag_idx == 0)
                        all_conic = false; break
                    end
                    neg_diag_idx = i
                else
                    all_conic = false; break
                end
            else
                if q.qcoeffs[i] == -1
                    if !(neg_diag_idx == off_diag_idx == 0)
                        all_conic = false; break
                    end
                    off_diag_idx = i
                else
                    all_conic = false; break
                end
            end
        end
        if n_pos_on_diag == nz-1 && neg_diag_idx > 0
            # Plain SOC
        elseif n_pos_on_diag == nz-1 && off_diag_idx > 0
            # Rotated SOC
        else
            all_conic = false
        end
    end
    return all_conic && length(m.quadconstr) > 0

end

# Converts the JuMP Model into a MathProgBase model based on the
# traits of the model
function buildInternalModel(m::Model, traits=ProblemTraits(m);
                            suppress_warnings=false, relaxation=false)
    # Set solver based on the model's traits if it hasn't provided
    if isa(m.solver, UnsetSolver)
        m.solver = default_solver(traits)
    end

    # If the model is nonlinear, use different logic in nlp.jl
    # to build the problem
    traits.nlp && return _buildInternalModel_nlp(m, traits)

    ## Temporary hack for Mosek, see https://github.com/JuliaOpt/Mosek.jl/pull/67
    if contains("$(typeof(m.solver))","MosekSolver") && isquadsoc(m)
        traits.conic = true
    end

    if traits.conic
        # If there are semicontinuous/semi-integer variables, we will have to
        # adjust the b vector below to construct a valid relaxation. This seems
        # like a pretty marginal case, so let's punt on it for now.
        if relaxation && any(x -> (x == :SemiCont || x == :SemiInt), m.colCat)
            error("Relaxations of conic problem with semi-integer/semicontinuous variables are not currently supported.")
        end

        # If the problem is conic then use only the objective
        # coefficients from prepProblemBounds
        f,_,_ = prepProblemBounds(m)

        # The conic MPB interface defines conic problems as
        # always being minimization problems, so flip if needed
        m.objSense == :Max && scale!(f, -1.0)

        # Obtain a fresh MPB model for the solver
        # If the problem is conic, we rebuild the problem from
        # scratch every time
        m.internalModel = MathProgBase.ConicModel(m.solver)

        # Build up the LHS, RHS and cones from the JuMP Model...
        A, b, var_cones, con_cones = conicconstraintdata(m)
        # ... and pass to the solver
        MathProgBase.loadproblem!(m.internalModel, f, A, b, con_cones, var_cones)
    else
        # Extract objective coefficients and linear constraint bounds
        f, rowlb, rowub = prepProblemBounds(m)
        # If we already have an MPB model for the solver...
        if m.internalModelLoaded
            # ... and if the solver supports updating bounds/objective
            if applicable(MathProgBase.setvarLB!, m.internalModel, m.colLower) &&
               applicable(MathProgBase.setvarUB!, m.internalModel, m.colUpper) &&
               applicable(MathProgBase.setconstrLB!, m.internalModel, rowlb) &&
               applicable(MathProgBase.setconstrUB!, m.internalModel, rowub) &&
               applicable(MathProgBase.setobj!, m.internalModel, f) &&
               applicable(MathProgBase.setsense!, m.internalModel, m.objSense)
                MathProgBase.setvarLB!(m.internalModel, copy(m.colLower))
                MathProgBase.setvarUB!(m.internalModel, copy(m.colUpper))
                MathProgBase.setconstrLB!(m.internalModel, rowlb)
                MathProgBase.setconstrUB!(m.internalModel, rowub)
                MathProgBase.setobj!(m.internalModel, f)
                MathProgBase.setsense!(m.internalModel, m.objSense)
            else
                # The solver doesn't support changing bounds/objective
                # We need to build the model from scratch
                if !suppress_warnings
                    Base.warn_once("Solver does not appear to support hot-starts. Model will be built from scratch.")
                end
                m.internalModelLoaded = false
            end
        end
        # If we don't already have a MPB model
        if !m.internalModelLoaded
            # Obtain a fresh MPB model for the solver
            m.internalModel = MathProgBase.LinearQuadraticModel(m.solver)
            # Construct a LHS matrix from the linear constraints
            A = prepConstrMatrix(m)

            # Load the problem data into the model...
            collb = copy(m.colLower)
            colub = copy(m.colUpper)
            if relaxation
                for i in 1:m.numCols
                    if m.colCat[i] in (:SemiCont,:SemiInt)
                        collb[i] = min(0.0, collb[i])
                        colub[i] = max(0.0, colub[i])
                    end
                end
            end
            MathProgBase.loadproblem!(m.internalModel, A, collb, colub, f, rowlb, rowub, m.objSense)
            # ... and add quadratic and SOS constraints separately
            addQuadratics(m)
            if !relaxation
                addSOS(m)
            end
        end
        # Update solver callbacks, if any
        if !relaxation
            registercallbacks(m)
        end
    end

    # Update the type of each variable
    if applicable(MathProgBase.setvartype!, m.internalModel, Symbol[])
        if relaxation
            MathProgBase.setvartype!(m.internalModel, fill(:Cont, m.numCols))
        else
            colCats = vartypes_without_fixed(m)
            MathProgBase.setvartype!(m.internalModel, colCats)
        end
    elseif traits.int
        # Solver that do not implement anything other than continuous
        # variables do not need to implement this method, so throw an
        # error if the model has anything but continuous
        error("Solver does not support discrete variables")
    end

    # Provide a primal solution to the solve, if the problem is integer
    # and the user has provided one.
    # TODO: change this so you can warm start continuous problems?
    if !relaxation && traits.int && !all(isnan(m.colVal))
        if applicable(MathProgBase.setwarmstart!, m.internalModel, m.colVal)
            MathProgBase.setwarmstart!(m.internalModel, m.colVal)
        else
            suppress_warnings || Base.warn_once("Solver does not appear to support providing initial feasible solutions.")
        end
    end

    # Some solvers need to have an explicit "update" phase, e.g. Gurobi
    if applicable(MathProgBase.updatemodel!, m.internalModel)
        MathProgBase.updatemodel!(m.internalModel)
    end

    # Record that we have a MPB model constructed
    m.internalModelLoaded = true
    nothing
end

# Add the quadratic part of the objective and all quadratic constraints
# to the internal MPB model
function addQuadratics(m::Model)
    # The objective function is always a quadratic expression, but
    # may have no quadratic terms (i.e. be just affine)
    if length(m.obj.qvars1) != 0
        # Check that no coefficients are NaN/Inf
        assert_isfinite(m.obj)
        # Check that quadratic term variables belong to this model
        # Affine portion is checked in prepProblemBounds
        if !(verify_ownership(m, m.obj.qvars1) &&
                verify_ownership(m, m.obj.qvars2))
            error("Variable not owned by model present in objective")
        end
        # Check for solver support for quadratic objectives happens in MPB
        MathProgBase.setquadobjterms!(m.internalModel,
            Cint[v.col for v in m.obj.qvars1],
            Cint[v.col for v in m.obj.qvars2], m.obj.qcoeffs)
    end

    # Add quadratic constraint to solver
    const sensemap = Dict(:(<=) => '<', :(>=) => '>', :(==) => '=')
    for k in 1:length(m.quadconstr)
        qconstr = m.quadconstr[k]::QuadConstraint
        if !haskey(sensemap, qconstr.sense)
            error("Invalid sense for quadratic constraint")
        end
        s = sensemap[qconstr.sense]

        terms::QuadExpr = qconstr.terms
        # Check that no coefficients are NaN/Inf
        assert_isfinite(terms)
        # Check that quadratic and affine term variables belong to this model
        if !(verify_ownership(m, terms.qvars1) &&
                verify_ownership(m, terms.qvars2) &&
                verify_ownership(m, terms.aff.vars))
            error("Variable not owned by model present in quadratic constraint")
        end
        # Extract indices for MPB, and add the constraint (if we can)
        affidx  = Cint[v.col for v in terms.aff.vars]
        var1idx = Cint[v.col for v in terms.qvars1]
        var2idx = Cint[v.col for v in terms.qvars2]
        if applicable(MathProgBase.addquadconstr!, m.internalModel, affidx, terms.aff.coeffs, var1idx, var2idx, terms.qcoeffs, s, -terms.aff.constant)
            MathProgBase.addquadconstr!(m.internalModel,
                affidx, terms.aff.coeffs,           # aᵀx +
                var1idx, var2idx, terms.qcoeffs,    # xᵀQx
                s, -terms.aff.constant)             # ≤/≥ b
        else
            error("Solver does not support quadratic constraints")
        end
    end
    nothing
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

# Returns coefficients for the affine part of the objective and the
# affine constraint lower and upper bounds, all as dense vectors
function prepProblemBounds(m::Model)

    # Create dense objective vector
    objaff::AffExpr = m.obj.aff
    # Check that no coefficients are NaN/Inf
    assert_isfinite(objaff)
    if !verify_ownership(m, objaff.vars)
        error("Variable not owned by model present in objective")
    end
    f = zeros(m.numCols)
    @inbounds for ind in 1:length(objaff.vars)
        f[objaff.vars[ind].col] += objaff.coeffs[ind]
    end

    # Create dense affine constraint bound vectors
    linconstr = m.linconstr::Vector{LinearConstraint}
    numRows = length(linconstr)
    # -Inf means no lower bound, +Inf means no upper bound
    rowlb = fill(-Inf, numRows)
    rowub = fill(+Inf, numRows)
    @inbounds for ind in 1:numRows
        rowlb[ind] = linconstr[ind].lb
        rowub[ind] = linconstr[ind].ub
    end

    return f, rowlb, rowub
end

# Convert all the affine constraints into a sparse column-wise
# matrix of coefficients.
function prepConstrMatrix(m::Model)

    linconstr = m.linconstr::Vector{LinearConstraint}
    numRows = length(linconstr)
    # Calculate the maximum number of nonzeros
    # The actual number may be less because of cancelling or
    # zero-coefficient terms
    nnz = 0
    for c in 1:numRows
        nnz += length(linconstr[c].terms.coeffs)
    end
    # Non-zero row indices
    I = Array(Int,nnz)
    # Non-zero column indices
    J = Array(Int,nnz)
    # Non-zero values
    V = Array(Float64,nnz)

    # Fill it up!
    # Number of nonzeros seen so far
    nnz = 0
    for c in 1:numRows
        # Check that no coefficients are NaN/Inf
        assert_isfinite(linconstr[c].terms)
        coeffs = linconstr[c].terms.coeffs
        vars   = linconstr[c].terms.vars
        # Check that variables belong to this model
        if !verify_ownership(m, vars)
            error("Variable not owned by model present in a constraint")
        end
        # Record all (i,j,v) triplets
        @inbounds for ind in 1:length(coeffs)
            nnz += 1
            I[nnz] = c
            J[nnz] = vars[ind].col
            V[nnz] = coeffs[ind]
        end
    end

    # sparse() handles merging duplicate terms and removing zeros
    A = sparse(I,J,V,numRows,m.numCols)
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
    var_cones = Any[cone for cone in m.varCones]
    con_cones = Any[]
    nnz = 0

    # find starting column indices for sdp matrices
    numSDPRows = 0
    numSymRows = 0
    for c in m.sdpconstr
        n = size(c.terms,1)
        @assert n == size(c.terms,2)
        @assert ndims(c.terms) == 2
        numSDPRows += convert(Int, n*(n+1)/2)
        for i in 1:n, j in i:n
            nnz += length(c.terms[i,j].coeffs)
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
        if !(all(t->t==0, q.aff.coeffs) && q.aff.constant == 0)
            error("Quadratic constraint $qconstr must be in second-order cone form")
        end
        n_pos_on_diag = 0
        off_diag_idx  = 0
        neg_diag_idx  = 0
        n = length(q.qvars1)
        nz = 0
        for i in 1:n
            q.qcoeffs[i] == 0 && continue
            nz += 1
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
        cone = Array(Int, nz)
        if n_pos_on_diag == nz-1 && neg_diag_idx > 0
            cone[1] = q.qvars1[neg_diag_idx].col
            r = 1
            for i in 1:n
                (q.qcoeffs[i] == 0 || i == neg_diag_idx) && continue
                r += 1
                cone[r] = q.qvars1[i].col
            end
            push!(soc_cones, cone)
        elseif n_pos_on_diag == nz-1 && off_diag_idx > 0
            cone[1] = q.qvars1[off_diag_idx].col
            cone[2] = q.qvars2[off_diag_idx].col
            r = 2
            for i in 1:n
                (q.qcoeffs[i] == 0 || i == off_diag_idx) && continue
                r += 1
                cone[r] = q.qvars1[i].col
            end
            push!(rsoc_cones, cone)
        else
            error("Quadratic constraint $qconstr is not conic representable")
        end
        numQuadRows += length(cone)
    end

    linconstr = m.linconstr::Vector{LinearConstraint}
    numLinRows = length(linconstr)
    numBounds = 0
    nonNeg  = Int[]
    nonPos  = Int[]
    free    = Int[]
    zeroVar = Int[]
    for i in 1:m.numCols
        seen = false
        lb, ub = m.colLower[i], m.colUpper[i]
        for (_,cone) in m.varCones
            if i in cone
                seen = true
                @assert lb == -Inf && ub == Inf
                break
            end
        end

        if !seen
            if lb != -Inf
                numBounds += 1
            end
            if ub != Inf
                numBounds += 1
            end
            if lb == 0 && ub == 0
                push!(zeroVar, i)
            elseif lb == 0
                push!(nonNeg, i)
            elseif ub == 0
                push!(nonPos, i)
            else
                push!(free, i)
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

    numSOCRows = 0
    numNormRows = 0
    for con in m.socconstr
        numNormRows += 1
        numSOCRows += length(con.normexpr.norm.terms) + 1
    end
    numRows = numLinRows + numBounds + numQuadRows + numSOCRows + numSDPRows + numSymRows

    # should maintain the order of constraints in the above form
    # throughout the code c is the conic constraint index
    # TODO: only added linear+bound+soc support, extend to all
    constr_dual_map = Array(Vector{Int}, numLinRows + numBounds + numNormRows)

    b = Array(Float64, numRows)

    I = Int[]
    J = Int[]
    V = Float64[]
    sizehint!(I, nnz)
    sizehint!(J, nnz)
    sizehint!(V, nnz)

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
        constr_dual_map[c] = collect(c)
    end

    c = numLinRows
    bndidx = 0
    for idx in 1:m.numCols
        lb = m.colLower[idx]
        if lb != -Inf
            bndidx += 1
            nnz += 1
            c   += 1
            push!(I, c)
            push!(J, idx)
            push!(V, 1.0)
            b[c] = lb
            push!(nonpos_rows, c)
            constr_dual_map[numLinRows + bndidx] = collect(c)
        end
        ub = m.colUpper[idx]
        if ub != Inf
            bndidx += 1
            c   += 1
            push!(I, c)
            push!(J, idx)
            push!(V, 1.0)
            b[c] = ub
            push!(nonneg_rows, c)
            constr_dual_map[numLinRows + bndidx] = collect(c)
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
        b[rng] = 0
        c += n
    end
    for cone in rsoc_cones
        n = length(cone)
        rng = (c+1):(c+n)
        append!(I, rng)
        append!(J, copy(cone))
        append!(V, [-1.0; -1.0; ones(n-2)])
        push!(con_cones, (:SOCRotated,rng))
        b[rng] = 0
        c += n
    end
    @assert c == numLinRows + numBounds + numQuadRows

    tmpelts = tmprow.elts
    tmpnzidx = tmprow.nzidx
    socidx = 0
    for con in m.socconstr
        socidx += 1
        expr = con.normexpr
        c += 1
        soc_start = c
        collect_expr!(m, tmprow, expr.aff)
        nnz = tmprow.nnz
        indices = tmpnzidx[1:nnz]
        append!(I, fill(c, nnz))
        append!(J, indices)
        append!(V, tmpelts[indices])
        b[c] = -expr.aff.constant
        for term in expr.norm.terms
            c += 1
            collect_expr!(m, tmprow, term)
            nnz = tmprow.nnz
            indices = tmpnzidx[1:nnz]
            append!(I, fill(c, nnz))
            append!(J, indices)
            append!(V, -expr.coeff*tmpelts[indices])
            b[c] = expr.coeff*term.constant
        end
        push!(con_cones, (:SOC, soc_start:c))
        constr_dual_map[numLinRows + numBounds + socidx] = collect(soc_start:c)
    end
    @assert c == numLinRows + numBounds + numQuadRows + numSOCRows

    for con in m.sdpconstr
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

    m.constrDualMap = constr_dual_map

    A = sparse(I, J, V, numRows, m.numCols)
    # @show full(A), b
    # @show var_cones, con_cones

    # TODO: uncomment these lines when they work with Mosek
    # supported = MathProgBase.supportedcones(m.internalModel)
    # @assert (:NonNeg in supported) && (:NonPos in supported) && (:Free in supported) && (:SDP in supported)
    A, b, var_cones, con_cones
end

getConstraintBounds(m::Model) = getConstraintBounds(m, ProblemTraits(m))

function getConstraintBounds(m::Model,traits::ProblemTraits)

    if traits.conic
        error("Not implemented for conic problems")
    elseif traits.sos
        error("Not implemented for SOS constraints")
    end

    linobj, linrowlb, linrowub = prepProblemBounds(m)

    quadrowlb = Float64[]
    quadrowub = Float64[]
    for c::QuadConstraint in m.quadconstr
        if c.sense == :(<=)
            push!(quadrowlb, -Inf)
            push!(quadrowub, 0.0)
        elseif c.sense == :(>=)
            push!(quadrowlb, 0.0)
            push!(quadrowub, Inf)
        else
            error("Unrecognized quadratic constraint sense $(c.sense)")
        end
    end

    nlrowlb = Float64[]
    nlrowub = Float64[]

    if traits.nlp
        nldata::NLPData = m.nlpdata
        for c in nldata.nlconstr
            push!(nlrowlb, c.lb)
            push!(nlrowub, c.ub)
        end
    end

    lb = [linrowlb;quadrowlb;nlrowlb]
    ub = [linrowub;quadrowub;nlrowub]

    return lb, ub

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

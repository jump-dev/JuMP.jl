typealias NonlinearConstraint GenericRangeConstraint{ReverseDiffSparse.SymbolicOutput}

type NLPData
    nlobj
    nlconstr::Vector{NonlinearConstraint}
    nlconstrlist::ReverseDiffSparse.ExprList
    evaluator
end

NLPData() = NLPData(nothing, NonlinearConstraint[], ExprList(), nothing)

Base.copy(::NLPData) = error("Copying nonlinear problems not yet implemented")

function initNLP(m::Model)
    if m.nlpdata === nothing
        m.nlpdata = NLPData()
    end
end

type JuMPNLPEvaluator <: MathProgBase.AbstractNLPEvaluator
    m::Model
    A::SparseMatrixCSC{Float64,Int} # linear constraint matrix
    has_nlobj::Bool
    linobj::Vector{Float64}
    eval_f_nl::Function
    eval_fgrad_nl::Function
    eval_h_nl::Function
    nnz_hess_obj::Int
    jac_I::Vector{Int}
    jac_J::Vector{Int}
    hess_I::Vector{Int}
    hess_J::Vector{Int}
    # timers
    eval_f_timer::Float64
    eval_g_timer::Float64
    eval_grad_f_timer::Float64
    eval_jac_g_timer::Float64
    eval_hesslag_timer::Float64
    function JuMPNLPEvaluator(m::Model, A)
        d = new(m,A)
        d.eval_f_timer = 0
        d.eval_g_timer = 0
        d.eval_grad_f_timer = 0
        d.eval_jac_g_timer = 0
        d.eval_hesslag_timer = 0
        return d
    end
end

function MathProgBase.initialize(d::JuMPNLPEvaluator, requested_features::Vector{Symbol})
    for feat in requested_features
        if !(feat in [:Grad, :Jac, :Hess, :ExprGraph])
            error("Unsupported feature $feat")
            # TODO: implement Jac-vec and Hess-vec products
            # for solvers that need them
        end
    end
    if d.eval_f_timer != 0
        # we've already been initialized
        # assume no new features are being requested.
        return
    end

    if :ExprGraph in requested_features
        need_expr = true
    else
        need_expr = false
    end

    initNLP(d.m) #in case the problem is purely linear/quadratic thus far
    nldata::NLPData = d.m.nlpdata
    d.has_nlobj = isa(nldata.nlobj, ReverseDiffSparse.SymbolicOutput)
    if d.has_nlobj
        @assert length(d.m.obj.qvars1) == 0 && length(d.m.obj.aff.vars) == 0
    end
    
    tic()

    d.linobj, linrowlb, linrowub = prepProblemBounds(d.m)

    A = d.A

    n_nlconstr = length(nldata.nlconstr)

    constrhessI, constrhessJ = prep_sparse_hessians(nldata.nlconstrlist, d.m.numCols, need_expr=need_expr)
    nljacI, nljacJ = jac_nz(nldata.nlconstrlist) # nonlinear jacobian components
    nnz_jac::Int = nnz(A) + length(nljacI)
    nnz_hess = length(constrhessI)

    for c in d.m.quadconstr
        nnz_jac += 2*length(c.terms.qvars1)+length(c.terms.aff.vars)
        nnz_hess += length(c.terms.qvars1)
    end
    
    if d.has_nlobj
        d.eval_fgrad_nl = genfgrad_simple(nldata.nlobj)
        d.eval_f_nl = genfval_simple(nldata.nlobj)
    end

    # Jacobian structure
    jac_I = zeros(nnz_jac)
    jac_J = zeros(nnz_jac)
    let
        idx = 1
        for col = 1:size(A,2)
            for pos = A.colptr[col]:(A.colptr[col+1]-1)
                jac_I[idx] = A.rowval[pos]
                jac_J[idx] = col
                idx += 1
            end
        end
        rowoffset = size(A,1)+1
        for c::QuadConstraint in d.m.quadconstr
            aff = c.terms.aff
            for k in 1:length(aff.vars)
                jac_I[idx] = rowoffset
                jac_J[idx] = aff.vars[k].col
                idx += 1
            end
            for k in 1:length(c.terms.qvars1)
                jac_I[idx] = rowoffset
                jac_J[idx] = c.terms.qvars1[k].col
                jac_I[idx+1] = rowoffset
                jac_J[idx+1] = c.terms.qvars2[k].col
                idx += 2
            end
            rowoffset += 1
        end
        for k in 1:length(nljacI)
            jac_I[idx] = nljacI[k]+rowoffset-1
            jac_J[idx] = nljacJ[k]
            idx += 1
        end
        @assert idx-1 == nnz_jac
    end
    d.jac_I = jac_I
    d.jac_J = jac_J



    if d.has_nlobj
        hI, hJ, d.eval_h_nl = gen_hessian_sparse_color_parametric(nldata.nlobj, d.m.numCols)
        nnz_hess += length(hI)
    else
        hI = []
        hJ = []
        d.eval_h_nl = (x,y,z) -> nothing
        nnz_hess += length(d.m.obj.qvars1)
    end
    d.nnz_hess_obj = length(hI)

    hess_I = zeros(nnz_hess)
    hess_J = zeros(nnz_hess)

    let
        qobj::QuadExpr = d.m.obj
        for i in 1:length(hI)
            # not guaranteed to be lower-triangular
            ix1 = nldata.nlobj.mapfromcanonical[hI[i]]
            ix2 = nldata.nlobj.mapfromcanonical[hJ[i]]
            if ix2 > ix1
                ix1, ix2 = ix2, ix1
            end
            hess_I[i] = ix1
            hess_J[i] = ix2
        end
        idx = length(hI)+1
        for k in 1:length(qobj.qvars1)
            qidx1 = qobj.qvars1[k].col
            qidx2 = qobj.qvars2[k].col
            if qidx2 > qidx1
                qidx1, qidx2 = qidx2, qidx1
            end
            hess_I[idx] = qidx1
            hess_J[idx] = qidx2
            idx += 1
        end
        # quadratic constraints
        for c::QuadConstraint in d.m.quadconstr
            for k in 1:length(c.terms.qvars1)
                qidx1 = c.terms.qvars1[k].col
                qidx2 = c.terms.qvars2[k].col
                if qidx2 > qidx1
                    qidx1, qidx2 = qidx2, qidx1
                end
                hess_I[idx] = qidx1
                hess_J[idx] = qidx2
                idx += 1
            end
        end

        hess_I[idx:end] = constrhessI
        hess_J[idx:end] = constrhessJ
    end

    d.hess_I = hess_I
    d.hess_J = hess_J

    numconstr = length(d.m.linconstr)+length(d.m.quadconstr)+n_nlconstr
    # call functions once to pre-compile
    MathProgBase.eval_f(d, d.m.colVal)
    MathProgBase.eval_grad_f(d, Array(Float64,d.m.numCols), d.m.colVal)
    MathProgBase.eval_g(d, Array(Float64,numconstr), d.m.colVal)
    MathProgBase.eval_jac_g(d, Array(Float64,nnz_jac), d.m.colVal)
    MathProgBase.eval_hesslag(d, Array(Float64,nnz_hess), d.m.colVal, 1.0, ones(numconstr))

    tprep = toq()
    #println("Prep time: $tprep")
    
    # reset timers
    d.eval_f_timer = 0
    d.eval_grad_f_timer = 0
    d.eval_g_timer = 0
    d.eval_jac_g_timer = 0
    d.eval_hesslag_timer = 0

    nothing
end

MathProgBase.features_available(d::JuMPNLPEvaluator) = [:Grad, :Jac, :Hess, :ExprGraph]

function MathProgBase.eval_f(d::JuMPNLPEvaluator, x)
    tic()
    if d.has_nlobj
        v = d.eval_f_nl(x)::Float64
    else
        qobj = d.m.obj::QuadExpr
        v = dot(x,d.linobj) + qobj.aff.constant
        for k in 1:length(qobj.qvars1)
            v += qobj.qcoeffs[k]*x[qobj.qvars1[k].col]*x[qobj.qvars2[k].col]
        end
    end
    d.eval_f_timer += toq()
    return v
end

function MathProgBase.eval_grad_f(d::JuMPNLPEvaluator, g, x)
    tic()
    if d.has_nlobj
        d.eval_fgrad_nl(x,g)
    else
        copy!(g,d.linobj)
        qobj::QuadExpr = d.m.obj
        for k in 1:length(qobj.qvars1)
            coef = qobj.qcoeffs[k]
            g[qobj.qvars1[k].col] += coef*x[qobj.qvars2[k].col]
            g[qobj.qvars2[k].col] += coef*x[qobj.qvars1[k].col]
        end
    end
    d.eval_grad_f_timer += toq()
    return
end

function MathProgBase.eval_g(d::JuMPNLPEvaluator, g, x)
    tic()
    A = d.A
    for i in 1:size(A,1); g[i] = 0.0; end
    #fill!(subarr(g,1:size(A,1)), 0.0)
    A_mul_B!(subarr(g,1:size(A,1)),A,x)
    idx = size(A,1)+1
    quadconstr = d.m.quadconstr::Vector{QuadConstraint}
    for c::QuadConstraint in quadconstr
        aff = c.terms.aff
        v = aff.constant
        for k in 1:length(aff.vars)
            v += aff.coeffs[k]*x[aff.vars[k].col]
        end
        for k in 1:length(c.terms.qvars1)
            v += c.terms.qcoeffs[k]*x[c.terms.qvars1[k].col]*x[c.terms.qvars2[k].col]
        end
        g[idx] = v
        idx += 1
    end
    eval_g!(subarr(g,idx:length(g)), (d.m.nlpdata::NLPData).nlconstrlist, x)
    
    d.eval_g_timer += toq()
    #print("x = ");show(x);println()
    #println(size(A,1), " g(x) = ");show(g);println()
    return
end

function MathProgBase.eval_jac_g(d::JuMPNLPEvaluator, J, x)
    tic()
    fill!(J,0.0)
    idx = 1
    A = d.A
    for col = 1:size(A,2)
        for pos = A.colptr[col]:(A.colptr[col+1]-1)
            J[idx] = A.nzval[pos]
            idx += 1
        end
    end
    quadconstr = d.m.quadconstr::Vector{QuadConstraint}
    for c::QuadConstraint in quadconstr
        aff = c.terms.aff
        for k in 1:length(aff.vars)
            J[idx] = aff.coeffs[k]
            idx += 1
        end
        for k in 1:length(c.terms.qvars1)
            coef = c.terms.qcoeffs[k]
            qidx1 = c.terms.qvars1[k].col
            qidx2 = c.terms.qvars2[k].col

            J[idx] = coef*x[qidx2]
            J[idx+1] = coef*x[qidx1]
            idx += 2
        end
    end
    eval_jac_g!(subarr(J,idx:length(J)), (d.m.nlpdata::NLPData).nlconstrlist, x)
    
    d.eval_jac_g_timer += toq()
    #print("x = ");show(x);println()
    #print("V ");show(J);println()
    return
end

function MathProgBase.eval_hesslag(
    d::JuMPNLPEvaluator,
    H::Vector{Float64},         # Sparse hessian entry vector
    x::Vector{Float64},         # Current solution
    obj_factor::Float64,        # Lagrangian multiplier for objective
    lambda::Vector{Float64})    # Multipliers for each constraint

    qobj = d.m.obj::QuadExpr
    nldata = d.m.nlpdata::NLPData
    
    tic()
    d.eval_h_nl(x, subarr(H, 1:d.nnz_hess_obj), nldata.nlobj)
    for i in 1:d.nnz_hess_obj; H[i] *= obj_factor; end
    #scale!(subarr(H, 1:length(hI)), obj_factor)
    # quadratic objective
    idx = 1+d.nnz_hess_obj
    for k in 1:length(qobj.qvars1)
        if qobj.qvars1[k].col == qobj.qvars2[k].col
            H[idx] = obj_factor*2*qobj.qcoeffs[k]
        else
            H[idx] = obj_factor*qobj.qcoeffs[k]
        end
        idx += 1
    end
    # quadratic constraints
    quadconstr = d.m.quadconstr::Vector{QuadConstraint}
    for i in 1:length(quadconstr)
        c = quadconstr[i]
        l = lambda[length(d.m.linconstr)+i]
        for k in 1:length(c.terms.qvars1)
            if c.terms.qvars1[k].col == c.terms.qvars2[k].col
                H[idx] = l*2*c.terms.qcoeffs[k]
            else
                H[idx] = l*c.terms.qcoeffs[k]
            end
            idx += 1
        end
    end

    eval_hess!(subarr(H, idx:length(H)), nldata.nlconstrlist, x, subarr(lambda, (length(d.m.linconstr::Vector{LinearConstraint})+length(quadconstr)+1):length(lambda)))
    d.eval_hesslag_timer += toq()
    return

end

MathProgBase.isobjlinear(d::JuMPNLPEvaluator) = !(isa(d.m.nlpdata.nlobj, ReverseDiffSparse.SymbolicOutput)) && (length(d.m.obj.qvars1) == 0)
# interpret quadratic to include purely linear
MathProgBase.isobjquadratic(d::JuMPNLPEvaluator) = !(isa(d.m.nlpdata.nlobj, ReverseDiffSparse.SymbolicOutput)) 

MathProgBase.isconstrlinear(d::JuMPNLPEvaluator, i::Integer) = (i <= length(d.m.linconstr))

MathProgBase.jac_structure(d::JuMPNLPEvaluator) = d.jac_I, d.jac_J
MathProgBase.hesslag_structure(d::JuMPNLPEvaluator) = d.hess_I, d.hess_J

# currently don't merge duplicates (this isn't required by MPB standard)
function affToExpr(aff::AffExpr, constant::Bool)
    ex = Expr(:call,:+)
    for k in 1:length(aff.vars)
        push!(ex.args, Expr(:call,:*,aff.coeffs[k],:(x[$(aff.vars[k].col)])))
    end
    if constant && aff.constant != 0
        push!(ex.args, aff.constant)
    end
    return ex
end

function quadToExpr(q::QuadExpr,constant::Bool)
    ex = Expr(:call,:+)
    for k in 1:length(q.qvars1)
        push!(ex.args, Expr(:call,:*,q.qcoeffs[k],:(x[$(q.qvars1[k].col)]), :(x[$(q.qvars2[k].col)])))
    end
    append!(ex.args, affToExpr(q.aff,constant).args[2:end])
    return ex
end

function MathProgBase.obj_expr(d::JuMPNLPEvaluator)
    if isa(d.m.nlpdata.nlobj, ReverseDiffSparse.SymbolicOutput)
        return ReverseDiffSparse.to_flat_expr(d.m.nlpdata.nlobj)
    else
        return quadToExpr(d.m.obj, true)
    end
end

function MathProgBase.constr_expr(d::JuMPNLPEvaluator,i::Integer)
    nlin = length(d.m.linconstr)
    nquad = length(d.m.quadconstr)
    if i <= nlin
        constr = d.m.linconstr[i]
        ex = affToExpr(constr.terms, false)
        if sense(constr) == :range
            return Expr(:comparison, constr.lb, ex, constr.ub)
        else
            return Expr(:comparison, ex, sense(constr), rhs(constr))
        end
    elseif i > nlin && i <= nlin + nquad
        i -= nlin
        qconstr = d.m.quadconstr[i]
        return Expr(:comparison, quadToExpr(qconstr.terms, true), qconstr.sense, 0)
    else
        i -= nlin + nquad
        ex = ReverseDiffSparse.to_flat_expr(d.m.nlpdata.nlconstrlist, i)
        return Expr(:comparison, ex, sense(d.m.nlpdata.nlconstr[i]), rhs(d.m.nlpdata.nlconstr[i]))
    end
end




function solvenlp(m::Model; suppress_warnings=false)
    
    linobj, linrowlb, linrowub = prepProblemBounds(m)

    nldata::NLPData = m.nlpdata
    if m.internalModelLoaded
        @assert isa(nldata.evaluator, JuMPNLPEvaluator)
        d = nldata.evaluator
    else
        A = prepConstrMatrix(m)
        d = JuMPNLPEvaluator(m,A)
        nldata.evaluator = d
    end

    numConstr = length(m.linconstr) + length(m.quadconstr) + length(nldata.nlconstr)

    nlrowlb = Float64[]
    nlrowub = Float64[]
    for c in nldata.nlconstr
        push!(nlrowlb, c.lb)
        push!(nlrowub, c.ub)
    end
    quadrowlb = Float64[]
    quadrowub = Float64[]
    for c::QuadConstraint in d.m.quadconstr
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
    
    #print("LB: ");show([linrowlb,nlrowlb]);println()
    #print("UB: ");show([linrowub,nlrowub]);println()

    m.internalModel = MathProgBase.model(m.solver)

    MathProgBase.loadnonlinearproblem!(m.internalModel, m.numCols, numConstr, m.colLower, m.colUpper, [linrowlb,quadrowlb,nlrowlb], [linrowub,quadrowub,nlrowub], m.objSense, d)
    if any(x->x!=:Cont, m.colCat)
        if applicable(MathProgBase.setvartype!, m.internalModel, m.colCat)
            MathProgBase.setvartype!(m.internalModel, m.colCat)
        else
            error("Solver does not support discrete variables")
        end
    end

    if !any(isnan(m.colVal))
        MathProgBase.setwarmstart!(m.internalModel, m.colVal)
    else
        # solve LP to find feasible point
        # do we need an iterior point?
        lpsol = MathProgBase.linprog(zeros(m.numCols), d.A, linrowlb, linrowub, m.colLower, m.colUpper)
        @assert lpsol.status == :Optimal
        MathProgBase.setwarmstart!(m.internalModel, lpsol.sol)
    end

    MathProgBase.optimize!(m.internalModel)
    stat = MathProgBase.status(m.internalModel)
    
    m.objVal = MathProgBase.getobjval(m.internalModel)
    m.colVal = MathProgBase.getsolution(m.internalModel)
    
    if stat != :Optimal && !suppress_warnings
        warn("Not solved to optimality, status: $stat")
    end

    m.internalModelLoaded = true

    #println("feval $(d.eval_f_timer)\nfgrad $(d.eval_grad_f_timer)\ngeval $(d.eval_g_timer)\njaceval $(d.eval_jac_g_timer)\nhess $(d.eval_hesslag_timer)")
    
    return stat

end





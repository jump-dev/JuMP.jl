#  Copyright 2015, Iain Dunning, Joey Huchette, Miles Lubin, and contributors
#  This Source Code Form is subject to the terms of the Mozilla Public
#  License, v. 2.0. If a copy of the MPL was not distributed with this
#  file, You can obtain one at http://mozilla.org/MPL/2.0/.

typealias NonlinearConstraint GenericRangeConstraint{ReverseDiffSparse.SymbolicOutput}

type NLPData
    nlobj
    nlconstr::Vector{NonlinearConstraint}
    nlconstrlist::ReverseDiffSparse.ExprList
    nlconstrDuals::Vector{Float64}
    evaluator
end

NLPData() = NLPData(nothing, NonlinearConstraint[], ExprList(), Float64[], nothing)

Base.copy(::NLPData) = error("Copying nonlinear problems not yet implemented")

function initNLP(m::Model)
    if m.nlpdata === nothing
        m.nlpdata = NLPData()
    end
end

function getDual(c::ConstraintRef{NonlinearConstraint})
    initNLP(c.m)
    nldata::NLPData = c.m.nlpdata
    if length(nldata.nlconstrDuals) != length(nldata.nlconstr)
        error("Dual solution not available. Check that the model was properly solved.")
    end
    return nldata.nlconstrDuals[c.idx]
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
    requested_hessian::Bool
    # timers
    eval_f_timer::Float64
    eval_g_timer::Float64
    eval_grad_f_timer::Float64
    eval_jac_g_timer::Float64
    eval_hesslag_timer::Float64
    function JuMPNLPEvaluator(m::Model)
        d = new(m)
        d.A = prepConstrMatrix(m)
        d.eval_f_timer = 0
        d.eval_g_timer = 0
        d.eval_grad_f_timer = 0
        d.eval_jac_g_timer = 0
        d.eval_hesslag_timer = 0
        return d
    end
end

@Base.deprecate JuMPNLPEvaluator(m::Model,A) JuMPNLPEvaluator(m)

function MathProgBase.initialize(d::JuMPNLPEvaluator, requested_features::Vector{Symbol})
    for feat in requested_features
        if !(feat in [:Grad, :Jac, :Hess, :HessVec, :ExprGraph])
            error("Unsupported feature $feat")
            # TODO: implement Jac-vec products
            # for solvers that need them
        end
    end
    if d.eval_f_timer != 0
        # we've already been initialized
        # assume no new features are being requested.
        return
    end

    initNLP(d.m) #in case the problem is purely linear/quadratic thus far
    nldata::NLPData = d.m.nlpdata

    if :ExprGraph in requested_features
        prep_expression_output(nldata.nlconstrlist)
        if length(requested_features) == 1 # don't need to do anything else
            return
        end
    end

    need_hessian = (:Hess in requested_features)
    d.requested_hessian = need_hessian

    d.has_nlobj = isa(nldata.nlobj, ReverseDiffSparse.SymbolicOutput)
    if d.has_nlobj
        @assert length(d.m.obj.qvars1) == 0 && length(d.m.obj.aff.vars) == 0
    end

    tic()

    d.linobj, linrowlb, linrowub = prepProblemBounds(d.m)

    A = d.A

    n_nlconstr = length(nldata.nlconstr)

    constrhessI, constrhessJ = prep_sparse_hessians(nldata.nlconstrlist, d.m.numCols, hessvec_only = !need_hessian)
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
    jac_I = zeros(Int,nnz_jac)
    jac_J = zeros(Int,nnz_jac)
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

    if need_hessian
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

        hess_I = zeros(Int,nnz_hess)
        hess_J = zeros(Int,nnz_hess)

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
        d.hess_I = hess_I
        d.hess_J = hess_J
    end

    numconstr = length(d.m.linconstr)+length(d.m.quadconstr)+n_nlconstr
    # call functions once to pre-compile
    MathProgBase.eval_f(d, d.m.colVal)
    MathProgBase.eval_grad_f(d, Array(Float64,d.m.numCols), d.m.colVal)
    MathProgBase.eval_g(d, Array(Float64,numconstr), d.m.colVal)
    MathProgBase.eval_jac_g(d, Array(Float64,nnz_jac), d.m.colVal)
    need_hessian && MathProgBase.eval_hesslag(d, Array(Float64,nnz_hess), d.m.colVal, 1.0, ones(numconstr))

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

MathProgBase.features_available(d::JuMPNLPEvaluator) = [:Grad, :Jac, :Hess, :HessVec, :ExprGraph]

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

import DualNumbers: Dual, epsilon

function MathProgBase.eval_hesslag_prod(
    d::JuMPNLPEvaluator,
    h::Vector{Float64}, # output vector
    x::Vector{Float64}, # current solution
    v::Vector{Float64}, # rhs vector
    σ::Float64,         # multiplier for objective
    μ::Vector{Float64}) # multipliers for each constraint

    nldata = d.m.nlpdata::NLPData

    # evaluate directional derivative of the gradient
    dualvec = reinterpret(Dual{Float64}, nldata.nlconstrlist.dualvec)
    dualout = reinterpret(Dual{Float64}, nldata.nlconstrlist.dualout)
    @assert length(dualvec) >= length(x)
    for i in 1:length(x)
        dualvec[i] = Dual(x[i], v[i])
        dualout[i] = zero(Dual{Float64})
    end
    MathProgBase.eval_grad_f(d, dualout, dualvec)
    for i in 1:length(x)
        h[i] = σ*epsilon(dualout[i])
    end

    row = size(d.A,1)+1
    quadconstr = d.m.quadconstr::Vector{QuadConstraint}
    for c in quadconstr
        l = μ[row]
        for k in 1:length(c.terms.qvars1)
            col1 = c.terms.qvars1[k].col
            col2 = c.terms.qvars2[k].col
            coef = c.terms.qcoeffs[k]
            if col1 == col2
                h[col1] += l*2*coef*v[col1]
            else
                h[col1] += l*coef*v[col2]
                h[col2] += l*coef*v[col1]
            end
        end
        row += 1
    end

    ReverseDiffSparse.eval_hessvec!(h, v, nldata.nlconstrlist, x, subarr(μ,row:length(μ)))

end

function MathProgBase.eval_hesslag(
    d::JuMPNLPEvaluator,
    H::Vector{Float64},         # Sparse hessian entry vector
    x::Vector{Float64},         # Current solution
    obj_factor::Float64,        # Lagrangian multiplier for objective
    lambda::Vector{Float64})    # Multipliers for each constraint

    qobj = d.m.obj::QuadExpr
    nldata = d.m.nlpdata::NLPData

    d.requested_hessian || error("Hessian computations were not requested on the call to MathProgBase.initialize.")

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
function MathProgBase.hesslag_structure(d::JuMPNLPEvaluator)
    d.requested_hessian || error("Hessian computations were not requested on the call to MathProgBase.initialize.")
    return d.hess_I, d.hess_J
end

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
            return Expr(:comparison, constr.lb, :(<=), ex, :(<=), constr.ub)
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
        constr = d.m.nlpdata.nlconstr[i]
        if sense(constr) == :range
            return Expr(:comparison, constr.lb, :(<=), ex, :(<=), constr.ub)
        else
            return Expr(:comparison, ex, sense(constr), rhs(constr))
        end
    end
end

function _buildInternalModel_nlp(m::Model, traits)

    @assert isempty(m.sdpconstr) && isempty(m.socconstr)

    linobj, linrowlb, linrowub = prepProblemBounds(m)

    nldata::NLPData = m.nlpdata
    if m.internalModelLoaded
        @assert isa(nldata.evaluator, JuMPNLPEvaluator)
        d = nldata.evaluator
    else
        d = JuMPNLPEvaluator(m)
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

    m.internalModel = MathProgBase.model(m.solver)

    MathProgBase.loadnonlinearproblem!(m.internalModel, m.numCols, numConstr, m.colLower, m.colUpper, [linrowlb;quadrowlb;nlrowlb], [linrowub;quadrowub;nlrowub], m.objSense, d)
    if traits.int
        if applicable(MathProgBase.setvartype!, m.internalModel, m.colCat)
            MathProgBase.setvartype!(m.internalModel, vartypes_without_fixed(m))
        else
            error("Solver does not support discrete variables")
        end
    end

    if !any(isnan,m.colVal)
        MathProgBase.setwarmstart!(m.internalModel, m.colVal)
    else
        initval = copy(m.colVal)
        initval[isnan(m.colVal)] = 0
        MathProgBase.setwarmstart!(m.internalModel, min(max(m.colLower,initval),m.colUpper))
    end

    m.internalModelLoaded = true

end


function solvenlp(m::Model, traits; suppress_warnings=false)

    @assert m.internalModelLoaded

    MathProgBase.optimize!(m.internalModel)
    stat = MathProgBase.status(m.internalModel)

    if stat != :Infeasible && stat != :Unbounded
        m.objVal = MathProgBase.getobjval(m.internalModel)
        m.colVal = MathProgBase.getsolution(m.internalModel)
    end

    if stat != :Optimal
        suppress_warnings || warn("Not solved to optimality, status: $stat")
    end
    if stat == :Optimal && !traits.int
        if applicable(MathProgBase.getconstrduals, m.internalModel) && applicable(MathProgBase.getreducedcosts, m.internalModel)
            nlduals = MathProgBase.getconstrduals(m.internalModel)
            m.linconstrDuals = nlduals[1:length(m.linconstr)]
            # quadratic duals currently not available, formulate as nonlinear constraint if needed
            m.nlpdata.nlconstrDuals = nlduals[length(m.linconstr)+length(m.quadconstr)+1:end]
            m.redCosts = MathProgBase.getreducedcosts(m.internalModel)
        else
            suppress_warnings || Base.warn_once("Nonlinear solver does not provide dual solutions")
        end
    end

    #println("feval $(d.eval_f_timer)\nfgrad $(d.eval_grad_f_timer)\ngeval $(d.eval_g_timer)\njaceval $(d.eval_jac_g_timer)\nhess $(d.eval_hesslag_timer)")

    return stat

end

# getValue for nonlinear subexpressions
@compat function getValue(x::Union{ReverseDiffSparse.ParametricExpressionWithParams,ReverseDiffSparse.ParametricExpression{0}})
    # messy check to extract model object
    found = false
    m = nothing
    for item in ReverseDiffSparse.expression_data(x)
        if isa(item, JuMPContainer)
            found = true
            m = getmeta(item, :model)
            break
        elseif isa(item, Array{Variable}) && !isempty(item)
            found = true
            m = first(item).m
        elseif isa(item, Variable)
            found = true
            m = item.m
            break
        end
    end
    found || error("Unable to determine which model this expression belongs to. Are there any variables present?")
    return ReverseDiffSparse.getvalue(x, m.colVal)
end

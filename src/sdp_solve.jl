parseScalarExpr(m::Model,d::MatrixFuncVar;debug=false) = parseScalarExpr(m,convert(ScalarExpr,d);debug=debug)
function parseScalarExpr(m::Model, d::ScalarExpr; debug=false)
    scalvaridx  = Int[x.col for x in d.aff.vars]
    scalcoefidx = d.aff.coeffs
    matdict = Dict{Int64,AbstractArray}()
    bnd_offset = d.aff.constant
    for (it,el) in enumerate(d.matvars)
        expr = el.expr
        sinfo = m.sdpdata.solverinfo[expr.elem[1].index]
        sgn = (sinfo.psd ? +1.0 : -1.0)
        if el.func == :trace
            all(x->isa(x,SDPVar),expr.elem) || error("Cannot have nested structure inside trace operator")
            mat = sgn*expr.post[1] * expr.pre[1] # exploit cyclic property of trace
            (!isa(mat,UniformScaling) && countnz(mat) == 0) && continue
            isa(mat,UniformScaling) && (mat *= speye(size(el.expr)...))
            debug && println("trace matrix: $mat")
        elseif el.func == :ref # post matrix and constant will be empty
            @assert expr.post[1] == 𝕀
            mat = sgn*expr.pre[1]
            countnz(mat) == 0 && continue
            debug && println("ref matrix: $mat")
        else
            error("Only trace operator or reference is currently supported")
        end
        if haskey(matdict, sinfo.id)
            matdict[sinfo.id] += mat
        else
            matdict[sinfo.id]  = mat
        end
        bnd_offset += sgn*trace(mat*sinfo.offset) + trace(el.expr.constant)
    end
    matvaridx  = Array(Int, length(keys(matdict)))
    matcoefidx = Array(Int, length(keys(matdict)))
    for (it,key) in enumerate(keys(matdict))
        matvaridx[it] = key
        idx = addsdpmatrix!(m.internalModel,matdict[key])
        matcoefidx[it] = idx
    end
    return scalvaridx, scalcoefidx, matvaridx, matcoefidx, bnd_offset
end

function addPrimalConstraint(m::Model,c::PrimalConstraint; debug=true)
    scalvaridx, scalcoefidx, matvaridx, matcoefidx, bnd_offset = 
        parseScalarExpr(m, c.terms;debug=debug)

    nexpr = c.terms.normexpr
    if length(nexpr.vars) > 0
        if c.lb == -Inf
            all(x->(x>=0), nexpr.coeffs) || error("No support for constraints of the form ||x|| >= rhs")
        elseif c.ub ==  Inf
            all(x->(x<=0), nexpr.coeffs) || error("No support for constraints of the form ||x|| >= rhs")
        else
            error("Cannot have equality or range constraints with normed terms")
        end
        for it in 1:length(nexpr.vars)
            if nexpr.form[it] == :norm2
                idx = addnorm2(m, nexpr.vars[it])
                push!(scalvaridx, idx)
                push!(scalcoefidx, nexpr.coeffs[it])
            elseif nexpr.form[it] == :normfrob
                idx = addnormfrob(m, nexpr.vars[it])
                push!(scalvaridx, idx)
                push!(scalcoefidx, nexpr.coeffs[it])
            elseif nexpr.form[it] == :abs
                idx = addabs(m, nexpr.vars[it])
                push!(scalvaridx, idx)
                push!(scalcoefidx, nexpr.coeffs[it])
            end   
        end
    end

    if debug
        println("matvaridx = $matvaridx")
        println("scalvaridx = $scalvaridx")
        println("scalcoefidx = $scalcoefidx")
        println("lb = $(c.lb-bnd_offset)")
        println("ub = $(c.ub-bnd_offset)")
    end
    addsdpconstr!(m.internalModel,
                 matvaridx,
                 matcoefidx,
                 scalvaridx,
                 scalcoefidx,
                 c.lb-bnd_offset,
                 c.ub-bnd_offset)
    return nothing
end

function addInternalVar(m::Model, dim::Int64)
    idx = addsdpvar!(m.internalModel, dim)
    var = SDPVar(m,length(m.sdpdata.sdpvar)+1,dim)
    push!(m.sdpdata.sdpvar, var)
    push!(m.sdpdata.lb, 0.0)
    push!(m.sdpdata.ub, Inf)
    push!(m.sdpdata.varname, "_internalvar")
    push!(m.sdpdata.solverinfo, SolverInfo(idx,true,spzeros(dim,dim)))
    return var
end

function addMatrixConstraint(model::Model,d::MatrixConstraint)
    # issym(d.terms) || error("Matrix expression must be symmetric")
    m, n = size(d.terms,1), size(d.terms,2)
    if d.sense == :(==)
        for i in 1:m, j in 1:n
            addPrimalConstraint(model, d.terms[i,j] == 0.0)
        end
    elseif d.sense == :(>=)
        _internalvar = addInternalVar(model,n)
        for i in 1:m, j in i:n
            addPrimalConstraint(model, d.terms[i,j] == _internalvar[i,j])
        end
    elseif d.sense == :(<=)
        _internalvar = addInternalVar(model,n)
        for i in 1:m, j in i:n
            addPrimalConstraint(model, -d.terms[i,j] == _internalvar[i,j])
        end
    elseif d.sense == :(.>=)
        for i in 1:m, j in 1:n
            addPrimalConstraint(model, d.terms[i,j] >= 0.0)
        end    
    elseif d.sense == :(.<=)
        for i in 1:m, j in 1:n
            addPrimalConstraint(model, d.terms[i,j] <= 0.0)
        end
    end
end

function addDualConstraint(m::Model, d::DualConstraint)
    issym(d.terms.constant) || error("Dual constant matrix must be symmetric")
    n  = size(d.terms.constant, 1)
    for c in d.terms.coeffs
        issym(c) || error("Dual constraint must be symmetric")
        size(c,1) == n || error("Coefficient matrices must be of compatible sizes")
    end
    if d.sense == :(==)  
        for i in 1:n, j in i:n
            con = d.terms.constant[i,j]
            coef = map(x->x[i,j], d.terms.coeffs)
            addconstr!(m.internalModel,[v.col for v in d.terms.vars],coef,-con,-con)
        end
    elseif d.sense == :(>=)
        _internalvar = addInternalVar(m,n)
        for i in 1:n, j in i:n
            terms = ScalarExpr(d.terms.constant[i,j])
            for it in 1:length(d.terms.vars)
                terms += d.terms.vars[it] * d.terms.coeffs[it][i,j]
            end
            addPrimalConstraint(m, terms ==   _internalvar[i,j] )
        end
    elseif d.sense == :(<=)
        _internalvar = addInternalVar(m,n)
        for i in 1:n, j in i:n
            for it in 1:length(d.terms.vars)
                terms += d.terms.vars[it] * d.terms.coeffs[it][i,j]
            end
            addPrimalConstraint(m, terms == -(_internalvar[i,j]))
        end
    elseif d.sense == :(.>=)
        _internalvar = addInternalVar(m,n)  
        for i in 1:n, j in i:n
            con = d.terms.constant[i,j]
            coef = map(x->x[i,j], d.terms.coeffs)
            addconstr!(m.internalModel,[v.col for v in d.terms.vars],coef,-con,Inf)
        end
    elseif d.sense == :(.<=)
        _internalvar = addInternalVar(m,n)  
        for i in 1:n, j in i:n
            con = d.terms.constant[i,j]
            coef = map(x->x[i,j], d.terms.coeffs)
            addconstr!(m.internalModel,[v.col for v in d.terms.vars],coef,-Inf,-con)
        end
    end
end

function addInternalScalarVar(m::Model, lb::Float64, ub::Float64)
    addvar!(m.internalModel, lb, ub, 0.0)
    # return numvar(m.internalModel) #Mosek.getnumvar(m.internalModel.task)
    return Mosek.getnumvar(m.internalModel.task)
end

addabs(m::Model, v::Variable) = addabs(m, convert(AffExpr,v))
function addabs(m::Model, aff::AffExpr)
    idx = addInternalScalarVar(m, 0.0, Inf)
    addconstr!(m.internalModel, vcat([x.col for x in aff.vars], idx),
                                vcat(aff.coeffs, -1.0),
                                -Inf,
                                0.0)
    addconstr!(m.internalModel, vcat([x.col for x in aff.vars], idx),
                                vcat(aff.coeffs, 1.0),
                                0.0,
                                Inf)
    return idx
end

function addnorm2(m::Model, ex::MatrixExpr)
    sx,sy = size(ex)
    (sx == 1 || sy == 1) || error("2-norm does not work on matrices")
    varlist = Int[]
    for i in 1:max(sx,sy)
        elem = ex[i]
        if isa(elem, Variable)
            push!(varlist, elem.col)
        elseif isa(elem, AffExpr)
            if length(elem.vars) == 1 && elem.coeffs[1] == 1.0 && elem.constant == 0.0
                # a Variable disguised as AffExpr
                push!(varlist, elem.vars[1].col)
                continue
            end
            _internalvar = addInternalScalarVar(m, -Inf, Inf)
            push!(varlist, _internalvar)
            addconstr!(m.internalModel,
                       vcat(Int[x.col for x in elem.vars], _internalvar), 
                       vcat(elem.coeffs, -1.0), 
                       -elem.constant, 
                       -elem.constant)
        elseif isa(elem, MatrixFuncVar) || isa(elem, ScalarExpr)
            _internalvar = addInternalScalarVar(m, -Inf, Inf)
            push!(varlist, _internalvar)
            scalvaridx, scalcoefidx, matvaridx, matcoefidx, bnd_offset = 
                parseScalarExpr(m,elem)
            addsdpconstr!(m.internalModel,
                         matvaridx,
                         matcoefidx,
                         vcat(scalvaridx, _internalvar),
                         vcat(scalcoefidx, -1.0),
                         -bnd_offset,
                         -bnd_offset)

        else 
            error("Unrecognized element of type $(typeof(elem))")
        end
    end
    _internalvar = addInternalScalarVar(m, 0.0, Inf)
    addquadconstr!(m.internalModel, Int[],
                                    Float64[],
                                    vcat(varlist, _internalvar),
                                    vcat(varlist, _internalvar),
                                    vcat(fill(1.0, length(varlist)), -1.0),
                                    '<',
                                    0.0)
    return _internalvar
end

function addnormfrob(m::Model, ex::MatrixExpr)
    sx,sy = size(ex)
    varlist = Int[]
    for i in 1:sx, j in 1:sy
        elem = ex[i,j]
        if isa(elem, Variable)
            push!(varlist, elem.col)
        elseif isa(elem, AffExpr)
            if length(elem.vars) == 1 && elem.coeffs[1] == 1.0 && elem.constant == 0.0
                # a Variable disguised as AffExpr
                push!(varlist, elem.vars[1].col)
                continue
            end
            _internalvar = addInternalScalarVar(m, -Inf, Inf)
            push!(varlist, _internalvar)
            addconstr!(m.internalModel,
                       vcat(Int[x.col for x in elem.vars], _internalvar), 
                       vcat(elem.coeffs, -1.0), 
                       -elem.constant, 
                       -elem.constant)
        elseif isa(elem, MatrixFuncVar) || isa(elem, ScalarExpr)
            _internalvar = addInternalScalarVar(m, -Inf, Inf)
            push!(varlist, _internalvar)
            scalvaridx, scalcoefidx, matvaridx, matcoefidx, bnd_offset = 
                parseScalarExpr(m,elem)
            addsdpconstr!(m.internalModel,
                         matvaridx,
                         matcoefidx,
                         vcat(scalvaridx, _internalvar),
                         vcat(scalcoefidx, -1.0),
                         -bnd_offset,
                         -bnd_offset)

        else 
            error("Unrecognized element of type $(typeof(elem))")
        end
    end
    _internalvar = addInternalScalarVar(m, 0.0, Inf)
    addquadconstr!(m.internalModel, Int[],
                                    Float64[],
                                    vcat(varlist, _internalvar),
                                    vcat(varlist, _internalvar),
                                    vcat(fill(1.0, length(varlist)), -1.0),
                                    '<',
                                    0.0)
    return _internalvar
end

function setupSDPVar(m::Model, it::Int64)
    var   = m.sdpdata.sdpvar[it]
    if isa(var, SDPVar)
        lb    = m.sdpdata.lb[it]
        ub    = m.sdpdata.ub[it]
        sinfo = m.sdpdata.solverinfo[it]
        sinfo.id = addsdpvar!(m.internalModel, var.dim)
        if lb == 0.0 || all(x->(x==0),lb)
            sinfo.psd = true
            sinfo.offset = spzeros(size(var)...)
            if ub == Inf || all(x->(x==Inf),ub) # X >= 0
                # do nothing
            else
                if ub == 0.0 || all(x->(x==0),ub) # X == 0
                    addMatrixConstraint(m, var == spzeros(size(var)...))
                else # 0 <= X <= C
                    addMatrixConstraint(m, var <= ub)
                end
            end
        elseif ub == 0.0 || all(x->(x==0),ub)
            sinfo.psd = false
            sinfo.offset = spzeros(size(var)...)
            if lb == -Inf || all(x->(x==-Inf),lb) # X <= 0
                # do nothing (here, at least)
            else # C <= X <= 0
                addMatrixConstraint(m, var >= lb)
            end
        else
            if lb == -Inf || all(x->(x==-Inf),lb) # X <= D
                sinfo.psd    = false
                sinfo.offset = ub
            elseif ub == Inf || all(x->(x==Inf),ub) # X >= C
                sinfo.psd    = true
                sinfo.offset = lb
            else # C <= X <= D
                sinfo.psd    = true
                sinfo.offset = lb
                addMatrixConstraint(m, var <= ub)
            end
        end
    end
end

function solveSDP(m::Model)
    for j = 1:m.numCols
        m.colCat[j] == INTEGER && error("Integer variables present in SDP problem")
    end

    sdp = m.sdpdata
    # make this solver-independent when CSDP is working
    m.solver = Mosek.MosekSolver()
    m.internalModel = model(m.solver)

    # add linear (scalar) constraints
    f, rowlb, rowub = prepProblemBounds(m)
    A = prepConstrMatrix(m)
    loadproblem!(m.internalModel, A, m.colLower, m.colUpper, f, rowlb, rowub, m.objSense)

    for it in 1:length(sdp.sdpvar)
        setupSDPVar(m, it)
    end

    # TODO: make this work for nested structure
    # TODO: deal with bounds changing in objective
    scalcost = zeros(Float64, m.numCols)
    for (it,vari) in enumerate(sdp.sdpobj.aff.vars)
        scalcost[vari.col] = sdp.sdpobj.aff.coeffs[it]
    end
    nexpr = sdp.sdpobj.normexpr
    auxcoef = Float64[]
    if length(nexpr.vars) > 0
        if m.objSense == :Min
            mapreduce(x->(x>=0), &, nexpr.coeffs) || error("Cannot minimize -||x||")
        elseif m.objSense == :Max
            mapreduce(x->(x<=0), &, nexpr.coeffs) || error("Cannot maximize ||x||")
        else
            error("Cannot have equality or range constraints with normed terms")
        end
        for it in 1:length(nexpr.vars)
            if nexpr.form[it] == :norm2
                idx = addnorm2(m, nexpr.vars[it])
                push!(auxcoef, nexpr.coeffs[it])
            elseif nexpr.form[it] == :normfrob
                idx = addnormfrob(m, nexpr.vars[it])
                push!(auxcoef, nexpr.coeffs[it])
            elseif nexpr.form[it] == :abs
                idx = addabs(m, nexpr.vars[it])
                push!(auxcoef, nexpr.coeffs[it])
            end   
        end
    end
    matvaridx  = Int[]
    matcoefidx = Int[]
    for (it,var) in enumerate(sdp.sdpobj.matvars)
        @assert length(var.expr.elem) == 1 # TODO: deal with nested structure
        sgn = sdp.solverinfo[it].psd ? +1.0 : -1.0
        mat = var.expr.pre[1]
        if isa(var.expr.pre[1], UniformScaling)
            idx = addsdpmatrix!(m.internalModel, sgn*mat.λ*speye(size(var.expr.elem[1])...))
        else
            idx = addsdpmatrix!(m.internalModel, sgn*mat) # TODO: deal with post case as well
        end
        push!(matcoefidx, idx)
        push!(matvaridx, sdp.solverinfo[it].id)
    end
    setsdpobj!(m.internalModel, matvaridx, matcoefidx)

    setobj!(m.internalModel, vcat(scalcost+f, auxcoef))
    setsense!(m.internalModel,m.objSense)

    # add quadratic terms
    addQuadratics(m)

    # add primal constraints
    for c in sdp.primalconstr
        addPrimalConstraint(m,c)
    end

    # add matrix constraints
    for d in sdp.matrixconstr
        addMatrixConstraint(m,d)
    end

    # add dual constraints
    for d in sdp.dualconstr
        addDualConstraint(m,d)
    end

    optimize!(m.internalModel)
    stat = status(m.internalModel)

    if stat == :NotSolved
        # do nothing
    elseif stat != :Optimal
        warn("SDP not solved to optimality, status: $stat")
    else
        m.colVal = MathProgBase.getsolution(m.internalModel)
        sdp.sdpval = Array(AbstractArray, length(sdp.sdpvar))
        for it in 1:length(sdp.sdpvar)
            if isa(sdp.sdpvar[it],SDPVar)
                idx    = sdp.solverinfo[it].id
                sgn    = sdp.solverinfo[it].psd ? +1.0 : -1.0
                offset = sdp.solverinfo[it].offset
                sdp.sdpval[it] = sgn*MathProgBase.getsdpsolution(m.internalModel, idx) + offset
            elseif isa(sdp.sdpvar[it],MatrixVar)
                sdp.sdpval[it] = m.colVal[sdp.solverinfo[it].id]
            end
        end
        m.objVal = MathProgBase.getobjval(m.internalModel) + m.obj.aff.constant + sdp.sdpobj.aff.constant
    end
    return stat
end

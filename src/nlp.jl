typealias NonlinearConstraint GenericRangeConstraint{ReverseDiffSparse.SymbolicOutput}

type NLPData
    nlobj
    nlconstr::Vector{NonlinearConstraint}
    nlconstrlist::ReverseDiffSparse.ExprList
end

NLPData() = NLPData(nothing, NonlinearConstraint[], ExprList())

copy(::NLPData) = error("Copying nonlinear problems not yet implemented")

function initNLP(m::Model)
    if m.nlpdata === nothing
        m.nlpdata = NLPData()
    end
end

if Pkg.installed("Ipopt") != nothing
    eval(Expr(:using,:Ipopt))
end

function solveIpopt(m::Model; options::Dict=Dict(), suppress_warnings=false)
    if Pkg.installed("Ipopt") === nothing
        error("Cannot solve nonlinear instances without Ipopt solver. Please run Pkg.add(\"Ipopt\").")
    end
    # check that there are no integer variables
    for j = 1:m.numCols
        if m.colCat[j] == INTEGER
            error("Integer variables present in nonlinear problem")
        end
    end
    tic()
    nldata::NLPData = m.nlpdata
    has_nlobj = isa(nldata.nlobj, ReverseDiffSparse.SymbolicOutput)
    if has_nlobj
        @assert length(m.obj.qvars1) == 0 && length(m.obj.aff.vars) == 0
    end
    objscale = 1.0
    if m.objSense == :Max
        objscale = -1.0
    end

    linobj, linrowlb, linrowub = prepProblemBounds(m)
    A = prepConstrMatrix(m)

    nlrowlb = Float64[]
    nlrowub = Float64[]
    n_nlconstr = length(nldata.nlconstr)

    constrhessI, constrhessJ = prep_sparse_hessians(nldata.nlconstrlist)
    jacI, jacJ = jac_nz(nldata.nlconstrlist)
    nnz_jac = length(A.nzval) + length(jacI)
    nnz_hess = length(constrhessI)
    for c in nldata.nlconstr
        push!(nlrowlb, c.lb)
        push!(nlrowub, c.ub)
    end
    quadrowlb = Float64[]
    quadrowub = Float64[]
    for c in m.quadconstr
        nnz_jac += 2*length(c.terms.qvars1)+length(c.terms.aff.vars)
        nnz_hess += length(c.terms.qvars1)
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
    #println("Prep time: $tprep")
    tf, tgf, tg, tjg, th = zeros(5)
    
    if has_nlobj
        fg = genfgrad_simple(nldata.nlobj)
        f = genfval_simple(nldata.nlobj)
        function eval_f(x)
            #print("x = ");show(x);println()
            #println("f(x) = ", f(x))
            tic()
            v = f(x)
            tf += toq()
            return objscale*v
        end

        function eval_grad_f(x, g)
            tic()
            fg(x,g)
            scale!(g,objscale)
            tgf += toq()
            #print("x = ");show(x);println()
            #println("gradf(x) = ");show(g);println()
        end
    else
        # linear and quadratic
        function eval_f(x)
            tic()
            v = dot(linobj,x) + m.obj.aff.constant 
            qobj::QuadExpr = m.obj
            for k in 1:length(qobj.qvars1)
                v += qobj.qcoeffs[k]*x[qobj.qvars1[k].col]*x[qobj.qvars2[k].col]
            end
            tf += toq()
            return objscale*v
        end

        function eval_grad_f(x,g)
            tic()
            copy!(g,linobj)
            qobj::QuadExpr = m.obj
            for k in 1:length(qobj.qvars1)
                coef = qobj.qcoeffs[k]
                g[qobj.qvars1[k].col] += coef*x[qobj.qvars2[k].col]
                g[qobj.qvars2[k].col] += coef*x[qobj.qvars1[k].col]
            end
            scale!(g,objscale)
            tf += toq()
        end
    end




    function eval_g(x, g)
        tic()
        fill!(sub(g,1:size(A,1)), 0.0)
        if VERSION < v"0.3-"
            g[1:size(A,1)] = A*x
        else
            A_mul_B!(sub(g,1:size(A,1)),A,x)
        end
        idx = size(A,1)+1
        for c::QuadConstraint in m.quadconstr
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
        eval_g!(sub(g,idx:length(g)), nldata.nlconstrlist, x)
        
        tg += toq()
        #print("x = ");show(x);println()
        #println(size(A,1), " g(x) = ");show(g);println()
    end

    function eval_jac_g(x, mode, rows, cols, values)
        if mode == :Structure
            # Convert column wise sparse to triple format
            idx = 1
            for col = 1:size(A,2)
                for pos = A.colptr[col]:(A.colptr[col+1]-1)
                    rows[idx] = A.rowval[pos]
                    cols[idx] = col
                    idx += 1
                end
            end
            rowoffset = size(A,1)+1
            for c::QuadConstraint in m.quadconstr
                aff = c.terms.aff
                for k in 1:length(aff.vars)
                    rows[idx] = rowoffset
                    cols[idx] = aff.vars[k].col
                    idx += 1
                end
                for k in 1:length(c.terms.qvars1)
                    rows[idx] = rowoffset
                    cols[idx] = c.terms.qvars1[k].col
                    rows[idx+1] = rowoffset
                    cols[idx+1] = c.terms.qvars2[k].col
                    idx += 2
                end
                rowoffset += 1
            end
            for k in 1:length(jacI)
                rows[idx] = jacI[k]+rowoffset-1
                cols[idx] = jacJ[k]
                idx += 1
            end
            @assert idx-1 == nnz_jac
            #print("I ");show(rows);println()
            #print("J ");show(cols);println()


        else
            # Values
            tic()
            fill!(values,0.0)
            idx = 1
            for col = 1:size(A,2)
                for pos = A.colptr[col]:(A.colptr[col+1]-1)
                    values[idx] = A.nzval[pos]
                    idx += 1
                end
            end
            for c::QuadConstraint in m.quadconstr
                aff = c.terms.aff
                for k in 1:length(aff.vars)
                    values[idx] = aff.coeffs[k]
                    idx += 1
                end
                for k in 1:length(c.terms.qvars1)
                    coef = c.terms.qcoeffs[k]
                    qidx1 = c.terms.qvars1[k].col
                    qidx2 = c.terms.qvars2[k].col

                    values[idx] = coef*x[qidx2]
                    values[idx+1] = coef*x[qidx1]
                    idx += 2
                end
            end
            eval_jac_g!(sub(values,idx:length(values)), nldata.nlconstrlist, x)
            
            tjg += toq()
            #print("x = ");show(x);println()
            #print("V ");show(values);println()
        end
    end

    if has_nlobj
        hI, hJ, hfunc = gen_hessian_sparse_color_parametric(nldata.nlobj)
        nnz_hess += length(hI)
    else
        hI = []
        hJ = []
        hfunc = (x,y,z) -> nothing
        nnz_hess += length(m.obj.qvars1)
    end

    function eval_h(
        x::Vector{Float64},         # Current solution
        mode,                       # Either :Structure or :Values
        rows::Vector{Int32},        # Sparsity structure - row indices
        cols::Vector{Int32},        # Sparsity structure - column indices
        obj_factor::Float64,        # Lagrangian multiplier for objective
        lambda::Vector{Float64},    # Multipliers for each constraint
        values::Vector{Float64})    # The values of the Hessian

        qobj::QuadExpr = m.obj
        if mode == :Structure
            for i in 1:length(hI)
                rows[i] = nldata.nlobj.mapfromcanonical[hI[i]]
                cols[i] = nldata.nlobj.mapfromcanonical[hJ[i]]
            end
            idx = length(hI)+1
            for k in 1:length(qobj.qvars1)
                qidx1 = qobj.qvars1[k].col
                qidx2 = qobj.qvars2[k].col
                if qidx2 > qidx1
                    qidx1, qidx2 = qidx2, qidx1
                end
                rows[idx] = qidx1
                cols[idx] = qidx2
                idx += 1
            end
            # quadratic constraints
            for c::QuadConstraint in m.quadconstr
                for k in 1:length(c.terms.qvars1)
                    qidx1 = c.terms.qvars1[k].col
                    qidx2 = c.terms.qvars2[k].col
                    if qidx2 > qidx1
                        qidx1, qidx2 = qidx2, qidx1
                    end
                    rows[idx] = qidx1
                    cols[idx] = qidx2
                    idx += 1
                end
            end

            rows[idx:end] = constrhessI
            cols[idx:end] = constrhessJ
        
        else
            tic()
            hfunc(x, sub(values, 1:length(hI)), nldata.nlobj)
            scale!(sub(values, 1:length(hI)), objscale*obj_factor)
            # quadratic objective
            idx = 1+length(hI)
            for k in 1:length(qobj.qvars1)
                if qobj.qvars1[k].col == qobj.qvars2[k].col
                    values[idx] = objscale*obj_factor*2*qobj.qcoeffs[k]
                else
                    values[idx] = objscale*obj_factor*qobj.qcoeffs[k]
                end
                idx += 1
            end
            # quadratic constraints
            for (i,c::QuadConstraint) in enumerate(m.quadconstr)
                l = lambda[length(m.linconstr)+i]
                for k in 1:length(c.terms.qvars1)
                    if c.terms.qvars1[k].col == c.terms.qvars2[k].col
                        values[idx] = l*2*c.terms.qcoeffs[k]
                    else
                        values[idx] = l*c.terms.qcoeffs[k]
                    end
                    idx += 1
                end
            end

            eval_hess!(sub(values, idx:length(values)), nldata.nlconstrlist, x, sub(lambda, (length(m.linconstr)+length(m.quadconstr)+1):length(lambda)))
            th += toq()

        end
    end
    #print("LB: ");show([linrowlb,nlrowlb]);println()
    #print("UB: ");show([linrowub,nlrowub]);println()
    tprep = toq()
    prob = createProblem(m.numCols, m.colLower, m.colUpper, length(m.linconstr)+length(m.quadconstr)+n_nlconstr,
        [linrowlb,quadrowlb,nlrowlb], [linrowub,quadrowub,nlrowub], nnz_jac, nnz_hess,
        eval_f, eval_g, eval_grad_f, eval_jac_g, eval_h)

    if !any(isnan(m.colVal))
        prob.x = m.colVal
    else
        # solve LP to find feasible point
        # do we need an iterior point?
        lpsol = linprog(zeros(m.numCols), A, linrowlb, linrowub, m.colLower, m.colUpper)
        @assert lpsol.status == :Optimal
        prob.x = lpsol.sol
    end

    # pass solver options to IPopt
    if !isempty(options)
        for (key,value) in options
            addOption(prob, key, value)
        end
    end
    status = Ipopt.ApplicationReturnStatus[solveProblem(prob)]
    m.colVal = prob.x
    m.objVal = objscale*prob.obj_val

    #println("feval $tf\nfgrad $tgf\ngeval $tg\njaceval $tjg\nhess $th")

    # translate status
    if status == :Solve_Succeeded || status == :Solved_To_Acceptable_Level
        return :Optimal
    elseif status == :Infeasible_Problem_Detected
        return :Infeasible
    else
        suppress_warnings || warn("Ipopt returned nonoptimal status: $status")
        return status
    end

end





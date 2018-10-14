#  Copyright 2017, Iain Dunning, Joey Huchette, Miles Lubin, and contributors
#  This Source Code Form is subject to the terms of the Mozilla Public
#  License, v. 2.0. If a copy of the MPL was not distributed with this
#  file, You can obtain one at http://mozilla.org/MPL/2.0/.


mutable struct NonlinearExprData
    nd::Vector{NodeData}
    const_values::Vector{Float64}
end

include("parsenlp.jl")

const NonlinearConstraint = GenericRangeConstraint{NonlinearExprData}

mutable struct NLPData
    nlobj
    nlconstr::Vector{NonlinearConstraint}
    nlexpr::Vector{NonlinearExprData}
    nlconstrDuals::Vector{Float64}
    nlparamvalues::Vector{Float64}
    user_operators::ReverseDiffSparse.UserOperatorRegistry
    largest_user_input_dimension::Int
    evaluator
end

function NonlinearExpression(m::Model,ex::NonlinearExprData)
    initNLP(m)
    nldata::NLPData = m.nlpdata
    push!(nldata.nlexpr, ex)
    return NonlinearExpression(m, length(nldata.nlexpr))
end

function newparameter(m::Model,value::Number)
    initNLP(m)
    nldata::NLPData = m.nlpdata
    push!(nldata.nlparamvalues, value)
    return NonlinearParameter(m, length(nldata.nlparamvalues))
end

getvalue(p::NonlinearParameter) = p.m.nlpdata.nlparamvalues[p.index]::Float64

setvalue(p::NonlinearParameter,v::Number) = (p.m.nlpdata.nlparamvalues[p.index] = v)

NLPData() = NLPData(nothing, NonlinearConstraint[], NonlinearExprData[], Float64[], Float64[], ReverseDiffSparse.UserOperatorRegistry(), 0, nothing)

Base.copy(::NLPData) = error("Copying nonlinear problems not yet implemented")

function initNLP(m::Model)
    if m.nlpdata === nothing
        m.nlpdata = NLPData()
    end
end

function getdual(c::ConstraintRef{Model,NonlinearConstraint})
    initNLP(c.m)
    nldata::NLPData = c.m.nlpdata
    if length(nldata.nlconstrDuals) != length(nldata.nlconstr)
        getdualwarn(c)
        NaN
    else
        nldata.nlconstrDuals[c.idx]
    end
end

mutable struct FunctionStorage
    nd::Vector{NodeData}
    adj::SparseMatrixCSC{Bool,Int}
    const_values::Vector{Float64}
    forward_storage::Vector{Float64}
    partials_storage::Vector{Float64}
    reverse_storage::Vector{Float64}
    grad_sparsity::Vector{Int}
    hess_I::Vector{Int} # nonzero pattern of hessian
    hess_J::Vector{Int}
    rinfo::Coloring.RecoveryInfo # coloring info for hessians
    seed_matrix::Matrix{Float64}
    linearity::Linearity
    dependent_subexpressions::Vector{Int} # subexpressions which this function depends on, ordered for forward pass
end

mutable struct SubexpressionStorage
    nd::Vector{NodeData}
    adj::SparseMatrixCSC{Bool,Int}
    const_values::Vector{Float64}
    forward_storage::Vector{Float64}
    partials_storage::Vector{Float64}
    reverse_storage::Vector{Float64}
    forward_storage_ϵ::Vector{Float64}
    partials_storage_ϵ::Vector{Float64}
    reverse_storage_ϵ::Vector{Float64}
    linearity::Linearity
end

mutable struct NLPEvaluator <: MathProgBase.AbstractNLPEvaluator
    m::Model
    A::SparseMatrixCSC{Float64,Int} # linear constraint matrix
    parameter_values::Vector{Float64}
    has_nlobj::Bool
    linobj::Vector{Float64}
    objective::FunctionStorage
    constraints::Vector{FunctionStorage}
    subexpressions::Vector{SubexpressionStorage}
    subexpression_order::Vector{Int}
    subexpression_forward_values::Vector{Float64}
    subexpression_reverse_values::Vector{Float64}
    subexpression_linearity::Vector{ReverseDiffSparse.Linearity}
    subexpressions_as_julia_expressions::Vector{Any}
    last_x::Vector{Float64}
    jac_storage::Vector{Float64} # temporary storage for computing jacobians
    user_output_buffer::Vector{Float64} # temporary storage for user-defined functions
    # storage for computing hessians
    # these Float64 vectors are reinterpreted to hold multiple epsilon components
    # so the length should be multiplied by the maximum number of epsilon components
    disable_2ndorder::Bool # don't offer Hess or HessVec
    want_hess::Bool
    forward_storage_ϵ::Vector{Float64} # (longest expression)
    partials_storage_ϵ::Vector{Float64} # (longest expression)
    reverse_storage_ϵ::Vector{Float64} # (longest expression)
    input_ϵ::Vector{Float64} # (number of variables)
    output_ϵ::Vector{Float64}# (number of variables)
    subexpression_forward_values_ϵ::Vector{Float64} # (number of subexpressions)
    subexpression_reverse_values_ϵ::Vector{Float64} # (number of subexpressions)
    # hessian sparsity pattern
    hess_I::Vector{Int}
    hess_J::Vector{Int}
    max_chunk::Int # chunk size for which we've allocated storage
    # init flags
    eval_f_init::Bool
    eval_g_init::Bool
    eval_grad_f_init::Bool
    eval_jac_g_init::Bool
    eval_hesslag_init::Bool
    function NLPEvaluator(m::Model)
        d = new(m)
        numVar = m.numCols
        d.A = prepConstrMatrix(m)
        # check if we have any user-defined operators, in which case we need to
        # disable hessians.
        if isa(m.nlpdata,NLPData)
            nldata::NLPData = m.nlpdata
            has_nlobj = isa(nldata.nlobj, NonlinearExprData)
            has_user_mv_operator = false
            for nlexpr in nldata.nlexpr
                has_user_mv_operator |= ReverseDiffSparse.has_user_multivariate_operators(nlexpr.nd)
            end
            if has_nlobj

                has_user_mv_operator |= ReverseDiffSparse.has_user_multivariate_operators(nldata.nlobj.nd)
            end
            for nlconstr in nldata.nlconstr
                has_user_mv_operator |= ReverseDiffSparse.has_user_multivariate_operators(nlconstr.terms.nd)
            end
            d.disable_2ndorder = has_user_mv_operator
            d.user_output_buffer = Array{Float64}(undef, m.nlpdata.largest_user_input_dimension)
            d.jac_storage = Array{Float64}(undef, max(numVar,m.nlpdata.largest_user_input_dimension))
        else
            d.disable_2ndorder = false
            d.user_output_buffer = Array{Float64}(undef, 0)
            d.jac_storage = Array{Float64}(undef, numVar)
        end

        d.eval_f_init = false
        d.eval_g_init = false
        d.eval_grad_f_init = false
        d.eval_jac_g_init = false
        d.eval_hesslag_init = false
        return d
    end
end

function simplify_expression(nd::Vector{NodeData}, const_values, subexpression_linearity, fixed_variables, parameter_values, x_values, subexpression_values)

    adj = adjmat(nd)
    forward_storage = zeros(length(nd))
    partials_storage = zeros(length(nd))
    linearity = classify_linearity(nd, adj, subexpression_linearity, fixed_variables)
    forward_eval(forward_storage, partials_storage, nd, adj, const_values, parameter_values, x_values, subexpression_values)
    nd_new = simplify_constants(forward_storage, nd, adj, const_values, linearity)
    return nd_new, forward_storage[1], linearity[1]
end

function FunctionStorage(nd::Vector{NodeData}, const_values,numVar, coloring_storage, want_hess::Bool, subexpr::Vector{Vector{NodeData}}, dependent_subexpressions, subexpression_linearity, subexpression_edgelist, subexpression_variables, fixed_variables)

    adj = adjmat(nd)
    forward_storage = zeros(length(nd))
    partials_storage = zeros(length(nd))
    reverse_storage = zeros(length(nd))
    empty!(coloring_storage)
    compute_gradient_sparsity!(coloring_storage, nd)

    for k in dependent_subexpressions
        compute_gradient_sparsity!(coloring_storage,subexpr[k])
    end
    grad_sparsity = sort!(collect(coloring_storage))
    empty!(coloring_storage)

    if want_hess
        # compute hessian sparsity
        linearity = classify_linearity(nd, adj, subexpression_linearity, fixed_variables)
        edgelist = compute_hessian_sparsity(nd, adj, linearity, coloring_storage, subexpression_edgelist, subexpression_variables)
        hess_I, hess_J, rinfo = Coloring.hessian_color_preprocess(edgelist, numVar, coloring_storage)
        seed_matrix = Coloring.seed_matrix(rinfo)
        if linearity[1] == NONLINEAR
            @assert length(hess_I) > 0
        end
    else
        hess_I = hess_J = Int[]
        rinfo = Coloring.RecoveryInfo()
        seed_matrix = Array{Float64}(undef, 0,0)
        linearity = [NONLINEAR]
    end

    return FunctionStorage(nd, adj, const_values, forward_storage, partials_storage, reverse_storage, grad_sparsity, hess_I, hess_J, rinfo, seed_matrix, linearity[1],dependent_subexpressions)

end

function SubexpressionStorage(nd::Vector{NodeData}, const_values,numVar, fixed_variables,subexpression_linearity)

    adj = adjmat(nd)
    forward_storage = zeros(length(nd))
    partials_storage = zeros(length(nd))
    reverse_storage = zeros(length(nd))
    linearity = classify_linearity(nd, adj, subexpression_linearity, fixed_variables)

    empty_arr = Array{Float64}(undef, 0)

    return SubexpressionStorage(nd, adj, const_values, forward_storage, partials_storage, reverse_storage, empty_arr, empty_arr, empty_arr, linearity[1])

end

function MathProgBase.initialize(d::NLPEvaluator, requested_features::Vector{Symbol})
    for feat in requested_features
        if !(feat in MathProgBase.features_available(d))
            error("Unsupported feature $feat")
            # TODO: implement Jac-vec products
            # for solvers that need them
        end
    end
    if d.eval_f_init
        # we've already been initialized
        # assume no new features are being requested.
        return
    end

    initNLP(d.m) #in case the problem is purely linear/quadratic thus far
    nldata::NLPData = d.m.nlpdata

    d.constraints = FunctionStorage[]
    d.last_x = fill(NaN, d.m.numCols)

    SIMPLIFY = d.m.simplify_nonlinear_expressions

    fixed_variables = SIMPLIFY ? d.m.colLower .== d.m.colUpper : d.m.colLower .== NaN
    d.m.colVal[fixed_variables] = d.m.colLower[fixed_variables]

    d.parameter_values = nldata.nlparamvalues

    d.linobj = prepAffObjective(d.m)
    linrowlb, linrowub = prepConstrBounds(d.m)
    numVar = length(d.linobj)

    d.want_hess = (:Hess in requested_features)
    want_hess_storage = (:HessVec in requested_features) || d.want_hess
    coloring_storage = ReverseDiffSparse.Coloring.IndexedSet(numVar)

    d.has_nlobj = isa(nldata.nlobj, NonlinearExprData)
    max_expr_length = 0
    main_expressions = Array{Vector{NodeData}}(undef, 0)
    subexpr = Array{Vector{NodeData}}(undef, 0)
    for nlexpr in nldata.nlexpr
        push!(subexpr, nlexpr.nd)
    end
    if d.has_nlobj
        push!(main_expressions,nldata.nlobj.nd)
    end
    for nlconstr in nldata.nlconstr
        push!(main_expressions,nlconstr.terms.nd)
    end
    d.subexpression_order, individual_order = order_subexpressions(main_expressions,subexpr)

    d.subexpression_linearity = Array{Linearity}(undef, length(nldata.nlexpr))
    subexpression_variables = Array{Vector{Int}}(undef, length(nldata.nlexpr))
    subexpression_edgelist = Array{Set{Tuple{Int,Int}}}(undef, length(nldata.nlexpr))
    d.subexpressions = Array{SubexpressionStorage}(undef, length(nldata.nlexpr))
    d.subexpression_forward_values = Array{Float64}(undef, length(d.subexpressions))
    d.subexpression_reverse_values = Array{Float64}(undef, length(d.subexpressions))

    empty_edgelist = Set{Tuple{Int,Int}}()
    for k in d.subexpression_order # only load expressions which actually are used
        if SIMPLIFY
            nd_new, forward_value, simplified_linearity = simplify_expression(nldata.nlexpr[k].nd, nldata.nlexpr[k].const_values,d.subexpression_linearity, fixed_variables, d.parameter_values, d.m.colVal, d.subexpression_forward_values)
        else
            nd_new = nldata.nlexpr[k].nd
            forward_value = NaN
            simplified_linearity = NONLINEAR
        end
        d.subexpression_forward_values[k] = forward_value
        if simplified_linearity != CONSTANT
            d.subexpressions[k] = SubexpressionStorage(nd_new, nldata.nlexpr[k].const_values, numVar, fixed_variables, d.subexpression_linearity)
            subex = d.subexpressions[k]
            d.subexpression_linearity[k] = subex.linearity
            @assert !SIMPLIFY || subex.linearity != CONSTANT
            if d.want_hess
                empty!(coloring_storage)
                compute_gradient_sparsity!(coloring_storage,subex.nd)
                # union with all dependent expressions
                for idx in list_subexpressions(subex.nd)
                    union!(coloring_storage, subexpression_variables[idx])
                end
                subexpression_variables[k] = collect(coloring_storage)
                empty!(coloring_storage)
                linearity = classify_linearity(subex.nd, subex.adj, d.subexpression_linearity, fixed_variables)
                edgelist = compute_hessian_sparsity(subex.nd, subex.adj, linearity,coloring_storage,subexpression_edgelist, subexpression_variables)
                subexpression_edgelist[k] = edgelist
            end
        else
            d.subexpression_linearity[k] = simplified_linearity
            subexpression_edgelist[k] = empty_edgelist
        end

    end

    if :ExprGraph in requested_features
        d.subexpressions_as_julia_expressions = Array{Any}(undef, length(subexpr))
        for k in d.subexpression_order
            if d.subexpression_linearity[k] != CONSTANT || !SIMPLIFY
                ex = d.subexpressions[k]
                d.subexpressions_as_julia_expressions[k] = tapeToExpr(d.m, 1, ex.nd, ex.adj, ex.const_values, d.parameter_values, d.subexpressions_as_julia_expressions, nldata.user_operators, true, true)
            else
                d.subexpressions_as_julia_expressions[k] = d.subexpression_forward_values[k]
            end
        end
    end

    if SIMPLIFY
        main_expressions = Array{Vector{NodeData}}(undef, 0)

        # simplify objective and constraint expressions
        if d.has_nlobj
            nd_new, forward_value = simplify_expression(nldata.nlobj.nd, nldata.nlobj.const_values,d.subexpression_linearity, fixed_variables, d.parameter_values, d.m.colVal, d.subexpression_forward_values)
            push!(main_expressions,nd_new)
        end
        for k in 1:length(nldata.nlconstr)
            nd_new, forward_value = simplify_expression(nldata.nlconstr[k].terms.nd, nldata.nlconstr[k].terms.const_values,d.subexpression_linearity, fixed_variables, d.parameter_values, d.m.colVal, d.subexpression_forward_values)
            push!(main_expressions,nd_new)
        end
        # recompute dependencies after simplification
        d.subexpression_order, individual_order = order_subexpressions(main_expressions,subexpr)

        subexpr = Array{Vector{NodeData}}(undef, length(d.subexpressions))
        for k in d.subexpression_order
            subexpr[k] = d.subexpressions[k].nd
        end
    end

    max_chunk = 1


    if d.has_nlobj
        @assert length(d.m.obj.qvars1) == 0 && length(d.m.obj.aff.vars) == 0
        nd = main_expressions[1]
        d.objective = FunctionStorage(nd, nldata.nlobj.const_values, numVar, coloring_storage, d.want_hess, subexpr, individual_order[1], d.subexpression_linearity, subexpression_edgelist, subexpression_variables, fixed_variables)
        max_expr_length = max(max_expr_length, length(d.objective.nd))
        max_chunk = max(max_chunk, size(d.objective.seed_matrix,2))
    end

    for k in 1:length(nldata.nlconstr)
        nlconstr = nldata.nlconstr[k]
        idx = (d.has_nlobj) ? k+1 : k
        nd = main_expressions[idx]
        push!(d.constraints, FunctionStorage(nd, nlconstr.terms.const_values, numVar, coloring_storage, d.want_hess, subexpr, individual_order[idx], d.subexpression_linearity, subexpression_edgelist, subexpression_variables, fixed_variables))
        max_expr_length = max(max_expr_length, length(d.constraints[end].nd))
        max_chunk = max(max_chunk, size(d.constraints[end].seed_matrix,2))
    end

    max_chunk = min(max_chunk, 10) # 10 is hardcoded upper bound to avoid excess memory allocation

    if d.want_hess || want_hess_storage # storage for Hess or HessVec
        d.input_ϵ = Array{Float64}(undef, max_chunk*d.m.numCols)
        d.output_ϵ = Array{Float64}(undef, max_chunk*d.m.numCols)
        d.forward_storage_ϵ = Array{Float64}(undef, max_chunk*max_expr_length)
        d.partials_storage_ϵ = Array{Float64}(undef, max_chunk*max_expr_length)
        d.reverse_storage_ϵ = Array{Float64}(undef, max_chunk*max_expr_length)
        d.subexpression_forward_values_ϵ = Array{Float64}(undef, max_chunk*length(d.subexpressions))
        d.subexpression_reverse_values_ϵ = Array{Float64}(undef, max_chunk*length(d.subexpressions))
        for k in d.subexpression_order
            subex = d.subexpressions[k]
            subex.forward_storage_ϵ = zeros(Float64,max_chunk*length(subex.nd))
            subex.partials_storage_ϵ = zeros(Float64,max_chunk*length(subex.nd))
            subex.reverse_storage_ϵ = zeros(Float64,max_chunk*length(subex.nd))
        end
        d.max_chunk = max_chunk
        if d.want_hess
            d.hess_I, d.hess_J = _hesslag_structure(d)
            # JIT warm-up
            MathProgBase.eval_hesslag(d, Array{Float64}(undef, length(d.hess_I)), d.m.colVal, 1.0, ones(MathProgBase.numconstr(d.m)))
        end
    end

    # JIT warm-up
    if :Grad in requested_features
        MathProgBase.eval_grad_f(d, zeros(numVar), d.m.colVal)
        MathProgBase.eval_g(d, zeros(MathProgBase.numconstr(d.m)), d.m.colVal)
    end

    #tprep = toq()
    #println("Prep time: $tprep")

    # init flags
    d.eval_f_init = false
    d.eval_grad_f_init = false
    d.eval_g_init = false
    d.eval_jac_g_init = false
    d.eval_hesslag_init = false

    nothing
end

function MathProgBase.features_available(d::NLPEvaluator)
    features = [:Grad, :Jac, :ExprGraph]
    if !d.disable_2ndorder
        push!(features,:Hess)
        push!(features,:HessVec)
    end
    return features
end

function forward_eval_all(d::NLPEvaluator,x)
    # do a forward pass on all expressions at x
    subexpr_values = d.subexpression_forward_values
    user_operators = d.m.nlpdata.user_operators::ReverseDiffSparse.UserOperatorRegistry
    user_input_buffer = d.jac_storage
    user_output_buffer = d.user_output_buffer
    SIMPLIFY = d.m.simplify_nonlinear_expressions
    for k in d.subexpression_order
        if SIMPLIFY && d.subexpression_linearity[k] == CONSTANT
            continue
        end
        ex = d.subexpressions[k]
        subexpr_values[k] = forward_eval(ex.forward_storage,ex.partials_storage,ex.nd,ex.adj,ex.const_values,d.parameter_values,x,subexpr_values,user_input_buffer,user_output_buffer, user_operators=user_operators)
    end
    if d.has_nlobj
        ex = d.objective
        forward_eval(ex.forward_storage,ex.partials_storage,ex.nd,ex.adj,ex.const_values,d.parameter_values,x,subexpr_values,user_input_buffer,user_output_buffer,user_operators=user_operators)
    end
    for ex in d.constraints
        forward_eval(ex.forward_storage,ex.partials_storage,ex.nd,ex.adj,ex.const_values,d.parameter_values,x,subexpr_values,user_input_buffer,user_output_buffer,user_operators=user_operators)
    end
end

function reverse_eval_all(d::NLPEvaluator,x)
    # do a reverse pass on all expressions at x
    subexpr_reverse_values = d.subexpression_reverse_values
    subexpr_values = d.subexpression_forward_values
    grad_storage = d.jac_storage
    SIMPLIFY = d.m.simplify_nonlinear_expressions
    for k in d.subexpression_order
        if SIMPLIFY && d.subexpression_linearity[k] == CONSTANT
            continue
        end
        ex = d.subexpressions[k]
        reverse_eval(ex.reverse_storage,ex.partials_storage,ex.nd,ex.adj)
    end
    if d.has_nlobj
        ex = d.objective
        reverse_eval(ex.reverse_storage,ex.partials_storage,ex.nd,ex.adj)
    end
    for ex in d.constraints
        reverse_eval(ex.reverse_storage,ex.partials_storage,ex.nd,ex.adj)
    end
    copyto!(d.last_x,x)
end

function MathProgBase.eval_f(d::NLPEvaluator, x)
    if d.last_x != x
        forward_eval_all(d,x)
        reverse_eval_all(d,x)
    end
    val = zero(eltype(x))
    if d.has_nlobj
        val = d.objective.forward_storage[1]
    else
        qobj = d.m.obj::QuadExpr
        val = dot(x,d.linobj) + qobj.aff.constant
        for k in 1:length(qobj.qvars1)
            val += qobj.qcoeffs[k]*x[qobj.qvars1[k].col]*x[qobj.qvars2[k].col]
        end
    end
    d.eval_f_init = true
    return val
end

function MathProgBase.eval_grad_f(d::NLPEvaluator, g, x)
    if d.last_x != x
        forward_eval_all(d,x)
        reverse_eval_all(d,x)
    end
    SIMPLIFY = d.m.simplify_nonlinear_expressions
    if d.has_nlobj
        fill!(g,0.0)
        ex = d.objective
        subexpr_reverse_values = d.subexpression_reverse_values
        subexpr_reverse_values[ex.dependent_subexpressions] .= 0.0
        reverse_extract(g,ex.reverse_storage,ex.nd,ex.adj,subexpr_reverse_values,1.0)
        for i in length(ex.dependent_subexpressions):-1:1
            k = ex.dependent_subexpressions[i]
            if SIMPLIFY && d.subexpression_linearity[k] == CONSTANT
                continue
            end
            subexpr = d.subexpressions[k]
            reverse_extract(g,subexpr.reverse_storage,subexpr.nd,subexpr.adj,subexpr_reverse_values,subexpr_reverse_values[k])

        end
    else
        copyto!(g,d.linobj)
        qobj::QuadExpr = d.m.obj
        for k in 1:length(qobj.qvars1)
            coef = qobj.qcoeffs[k]
            g[qobj.qvars1[k].col] += coef*x[qobj.qvars2[k].col]
            g[qobj.qvars2[k].col] += coef*x[qobj.qvars1[k].col]
        end
    end
    d.eval_grad_f_init = true
    return
end

function MathProgBase.eval_g(d::NLPEvaluator, g, x)
    if d.last_x != x
        forward_eval_all(d,x)
        reverse_eval_all(d,x)
    end
    A = d.A
    for i in 1:size(A,1); g[i] = 0.0; end
    #fill!(view(g,1:size(A,1)), 0.0)
    if VERSION >= v"0.7-"
        g[1:size(A,1)] .= A * x
    else
        A_mul_B!(view(g,1:size(A,1)),A,x)
    end
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
    for ex in d.constraints
        g[idx] = ex.forward_storage[1]
        idx += 1
    end

    d.eval_g_init = true
    #print("x = ");show(x);println()
    #println(size(A,1), " g(x) = ");show(g);println()
    return
end

function MathProgBase.eval_jac_g(d::NLPEvaluator, J, x)
    if d.last_x != x
        forward_eval_all(d,x)
        reverse_eval_all(d,x)
    end
    fill!(J,0.0)
    idx = 1
    A = d.A
    for col = 1:size(A,2)
        for pos = nzrange(A,col)
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
    grad_storage = d.jac_storage
    subexpr_reverse_values = d.subexpression_reverse_values
    SIMPLIFY = d.m.simplify_nonlinear_expressions
    for ex in d.constraints
        nzidx = ex.grad_sparsity
        grad_storage[nzidx] .= 0.0
        subexpr_reverse_values[ex.dependent_subexpressions] .= 0.0

        reverse_extract(grad_storage,ex.reverse_storage,ex.nd,ex.adj,subexpr_reverse_values,1.0)
        for i in length(ex.dependent_subexpressions):-1:1
            k = ex.dependent_subexpressions[i]
            if SIMPLIFY && d.subexpression_linearity[k] == CONSTANT
                continue
            end
            subexpr = d.subexpressions[k]
            reverse_extract(grad_storage,subexpr.reverse_storage,subexpr.nd,subexpr.adj,subexpr_reverse_values,subexpr_reverse_values[k])
        end

        for k in 1:length(nzidx)
            J[idx+k-1] = grad_storage[nzidx[k]]
        end
        idx += length(nzidx)
    end

    d.eval_jac_g_init = true
    #print("x = ");show(x);println()
    #print("V ");show(J);println()
    return
end



function MathProgBase.eval_hesslag_prod(
    d::NLPEvaluator,
    h::AbstractVector{Float64}, # output vector
    x::AbstractVector{Float64}, # current solution
    v::AbstractVector{Float64}, # rhs vector
    σ::Float64,                 # multiplier for objective
    μ::AbstractVector{Float64}) # multipliers for each constraint

    nldata = d.m.nlpdata::NLPData

    if d.last_x != x
        forward_eval_all(d,x)
        reverse_eval_all(d,x)
    end

    fill!(h, 0.0)

    # quadratic objective
    qobj::QuadExpr = d.m.obj
    for k in 1:length(qobj.qvars1)
        col1 = qobj.qvars1[k].col
        col2 = qobj.qvars2[k].col
        coef = qobj.qcoeffs[k]
        if col1 == col2
            h[col1] += σ*2*coef*v[col1]
        else
            h[col1] += σ*coef*v[col2]
            h[col2] += σ*coef*v[col1]
        end
    end

    # quadratic constraints
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

    input_ϵ = reinterpret(ForwardDiff.Partials{1,Float64}, d.input_ϵ)
    output_ϵ = reinterpret(ForwardDiff.Partials{1,Float64}, d.output_ϵ)
    for i in 1:length(x)
        input_ϵ[i] = ForwardDiff.Partials((v[i],))
    end

    # forward evaluate all subexpressions once
    subexpr_forward_values_ϵ = reinterpret(ForwardDiff.Partials{1,Float64},d.subexpression_forward_values_ϵ)
    subexpr_reverse_values_ϵ = reinterpret(ForwardDiff.Partials{1,Float64},d.subexpression_reverse_values_ϵ)
    forward_storage_ϵ = reinterpret(ForwardDiff.Partials{1,Float64},d.forward_storage_ϵ)
    reverse_storage_ϵ = reinterpret(ForwardDiff.Partials{1,Float64},d.reverse_storage_ϵ)
    partials_storage_ϵ = reinterpret(ForwardDiff.Partials{1,Float64},d.partials_storage_ϵ)
    SIMPLIFY = d.m.simplify_nonlinear_expressions
    for expridx in d.subexpression_order
        if SIMPLIFY && d.subexpression_linearity[expridx] == CONSTANT
            continue
        end
        subexpr = d.subexpressions[expridx]
        sub_forward_storage_ϵ = reinterpret(ForwardDiff.Partials{1,Float64},subexpr.forward_storage_ϵ)
        sub_partials_storage_ϵ = reinterpret(ForwardDiff.Partials{1,Float64},subexpr.partials_storage_ϵ)
        subexpr_forward_values_ϵ[expridx] = forward_eval_ϵ(subexpr.forward_storage,sub_forward_storage_ϵ, subexpr.partials_storage, sub_partials_storage_ϵ, subexpr.nd, subexpr.adj, input_ϵ, subexpr_forward_values_ϵ,user_operators=nldata.user_operators)
    end
    # we only need to do one reverse pass through the subexpressions as well
    zero_ϵ = zero(ForwardDiff.Partials{1,Float64})
    fill!(subexpr_reverse_values_ϵ,zero_ϵ)
    fill!(d.subexpression_reverse_values,0.0)
    fill!(reverse_storage_ϵ,zero_ϵ)
    fill!(output_ϵ,zero_ϵ)
    if d.has_nlobj
        ex = d.objective
        forward_eval_ϵ(ex.forward_storage, forward_storage_ϵ, ex.partials_storage, partials_storage_ϵ, ex.nd,ex.adj,input_ϵ, subexpr_forward_values_ϵ,user_operators=nldata.user_operators)
        reverse_eval_ϵ(output_ϵ,ex.reverse_storage, reverse_storage_ϵ,ex.partials_storage, partials_storage_ϵ,ex.nd,ex.adj,d.subexpression_reverse_values,subexpr_reverse_values_ϵ, σ, zero_ϵ)
    end


    for i in 1:length(d.constraints)
        ex = d.constraints[i]
        l = μ[row]
        forward_eval_ϵ(ex.forward_storage, forward_storage_ϵ, ex.partials_storage, partials_storage_ϵ, ex.nd,ex.adj, input_ϵ,subexpr_forward_values_ϵ,user_operators=nldata.user_operators)
        reverse_eval_ϵ(output_ϵ, ex.reverse_storage, reverse_storage_ϵ, ex.partials_storage, partials_storage_ϵ, ex.nd,ex.adj,d.subexpression_reverse_values,subexpr_reverse_values_ϵ, l, zero_ϵ)
        row += 1
    end

    for i in length(d.subexpression_order):-1:1
        expridx = d.subexpression_order[i]
        if SIMPLIFY && d.subexpression_linearity[expridx] == CONSTANT
            continue
        end
        subexpr = d.subexpressions[expridx]
        sub_reverse_storage_ϵ = reinterpret(ForwardDiff.Partials{1,Float64},subexpr.reverse_storage_ϵ)
        sub_partials_storage_ϵ = reinterpret(ForwardDiff.Partials{1,Float64},subexpr.partials_storage_ϵ)
        reverse_eval_ϵ(output_ϵ,subexpr.reverse_storage,sub_reverse_storage_ϵ, subexpr.partials_storage, sub_partials_storage_ϵ,subexpr.nd,subexpr.adj,d.subexpression_reverse_values,subexpr_reverse_values_ϵ,d.subexpression_reverse_values[expridx],subexpr_reverse_values_ϵ[expridx])
    end

    for i in 1:length(x)
        h[i] += output_ϵ[i].values[1]
    end

end

function MathProgBase.eval_hesslag(
    d::NLPEvaluator,
    H::AbstractVector{Float64},         # Sparse hessian entry vector
    x::AbstractVector{Float64},         # Current solution
    obj_factor::Float64,                # Lagrangian multiplier for objective
    lambda::AbstractVector{Float64})    # Multipliers for each constraint

    qobj = d.m.obj::QuadExpr
    nldata = d.m.nlpdata::NLPData

    d.want_hess || error("Hessian computations were not requested on the call to MathProgBase.initialize.")

    if d.last_x != x
        forward_eval_all(d,x)
        reverse_eval_all(d,x)
    end

    # quadratic objective
    nzcount = 1
    for k in 1:length(qobj.qvars1)
        if qobj.qvars1[k].col == qobj.qvars2[k].col
            H[nzcount] = obj_factor*2*qobj.qcoeffs[k]
        else
            H[nzcount] = obj_factor*qobj.qcoeffs[k]
        end
        nzcount += 1
    end
    # quadratic constraints
    quadconstr = d.m.quadconstr::Vector{QuadConstraint}
    for i in 1:length(quadconstr)
        c = quadconstr[i]
        l = lambda[length(d.m.linconstr)+i]
        for k in 1:length(c.terms.qvars1)
            if c.terms.qvars1[k].col == c.terms.qvars2[k].col
                H[nzcount] = l*2*c.terms.qcoeffs[k]
            else
                H[nzcount] = l*c.terms.qcoeffs[k]
            end
            nzcount += 1
        end
    end

    fill!(d.input_ϵ,0.0)
    recovery_tmp_storage = d.output_ϵ
    nzcount -= 1

    if d.has_nlobj
        ex = d.objective
        chunk = min(size(ex.seed_matrix,2),d.max_chunk)
        if chunk == 1
            # skip dynamic dispatch
            nzthis = hessian_slice(d, ex, x, H, obj_factor, nzcount, recovery_tmp_storage, Val{1})::Int
        else
            nzthis = hessian_slice(d, ex, x, H, obj_factor, nzcount, recovery_tmp_storage, Val{chunk})::Int
        end
        nzcount += nzthis
    end

    for i in 1:length(d.constraints)
        ex = d.constraints[i]
        chunk = min(size(ex.seed_matrix,2),d.max_chunk)
        if chunk == 1
            nzthis = hessian_slice(d, ex, x, H, lambda[i+length(quadconstr)+length(d.m.linconstr)], nzcount, recovery_tmp_storage, Val{1})::Int
        else
            nzthis = hessian_slice(d, ex, x, H, lambda[i+length(quadconstr)+length(d.m.linconstr)], nzcount, recovery_tmp_storage, Val{chunk})::Int
        end
        nzcount += nzthis
    end

    d.eval_hesslag_init = true
    return

end

function hessian_slice_inner(d, ex, R, input_ϵ, output_ϵ, ::Type{Val{CHUNK}}) where CHUNK

    subexpr_forward_values_ϵ = reinterpret_unsafe(ForwardDiff.Partials{CHUNK,Float64},d.subexpression_forward_values_ϵ)
    subexpr_reverse_values_ϵ = reinterpret_unsafe(ForwardDiff.Partials{CHUNK,Float64},d.subexpression_reverse_values_ϵ)
    forward_storage_ϵ = reinterpret_unsafe(ForwardDiff.Partials{CHUNK,Float64},d.forward_storage_ϵ)
    reverse_storage_ϵ = reinterpret_unsafe(ForwardDiff.Partials{CHUNK,Float64},d.reverse_storage_ϵ)
    partials_storage_ϵ = reinterpret_unsafe(ForwardDiff.Partials{CHUNK,Float64},d.partials_storage_ϵ)
    zero_ϵ = zero(ForwardDiff.Partials{CHUNK,Float64})

    user_operators = d.m.nlpdata.user_operators::ReverseDiffSparse.UserOperatorRegistry
    SIMPLIFY = d.m.simplify_nonlinear_expressions
    # do a forward pass
    for expridx in ex.dependent_subexpressions
        if SIMPLIFY && d.subexpression_linearity[expridx] == CONSTANT
            continue
        end
        subexpr = d.subexpressions[expridx]
        sub_forward_storage_ϵ = reinterpret_unsafe(ForwardDiff.Partials{CHUNK,Float64},subexpr.forward_storage_ϵ)
        sub_partials_storage_ϵ = reinterpret_unsafe(ForwardDiff.Partials{CHUNK,Float64},subexpr.partials_storage_ϵ)
        subexpr_forward_values_ϵ[expridx] = forward_eval_ϵ(subexpr.forward_storage,sub_forward_storage_ϵ,subexpr.partials_storage,sub_partials_storage_ϵ, subexpr.nd, subexpr.adj, input_ϵ, subexpr_forward_values_ϵ,user_operators=user_operators)
    end
    forward_eval_ϵ(ex.forward_storage,forward_storage_ϵ,ex.partials_storage, partials_storage_ϵ,ex.nd,ex.adj,input_ϵ, subexpr_forward_values_ϵ,user_operators=user_operators)

    # do a reverse pass
    subexpr_reverse_values_ϵ[ex.dependent_subexpressions] = zero_ϵ
    d.subexpression_reverse_values[ex.dependent_subexpressions] .= 0.0

    reverse_eval_ϵ(output_ϵ, ex.reverse_storage, reverse_storage_ϵ,ex.partials_storage, partials_storage_ϵ,ex.nd,ex.adj,d.subexpression_reverse_values,subexpr_reverse_values_ϵ, 1.0, zero_ϵ)
    for i in length(ex.dependent_subexpressions):-1:1
        expridx = ex.dependent_subexpressions[i]
        if SIMPLIFY && d.subexpression_linearity[expridx] == CONSTANT
            continue
        end
        subexpr = d.subexpressions[expridx]
        sub_reverse_storage_ϵ = reinterpret_unsafe(ForwardDiff.Partials{CHUNK,Float64},subexpr.reverse_storage_ϵ)
        sub_partials_storage_ϵ = reinterpret_unsafe(ForwardDiff.Partials{CHUNK,Float64},subexpr.partials_storage_ϵ)
        reverse_eval_ϵ(output_ϵ, subexpr.reverse_storage, sub_reverse_storage_ϵ,subexpr.partials_storage,sub_partials_storage_ϵ,subexpr.nd,subexpr.adj,d.subexpression_reverse_values,subexpr_reverse_values_ϵ,d.subexpression_reverse_values[expridx],subexpr_reverse_values_ϵ[expridx])
    end
end

function hessian_slice(d, ex, x, H, scale, nzcount, recovery_tmp_storage,::Type{Val{CHUNK}}) where CHUNK

    nzthis = length(ex.hess_I)
    if ex.linearity == LINEAR
        @assert nzthis == 0
        return 0
    end
    R = ex.seed_matrix
    Coloring.prepare_seed_matrix!(R,ex.rinfo)
    local_to_global_idx = ex.rinfo.local_indices

    zero_ϵ = zero(ForwardDiff.Partials{CHUNK,Float64})

    input_ϵ_raw = d.input_ϵ
    output_ϵ_raw = d.output_ϵ
    input_ϵ = reinterpret_unsafe(ForwardDiff.Partials{CHUNK,Float64}, input_ϵ_raw)
    output_ϵ = reinterpret_unsafe(ForwardDiff.Partials{CHUNK,Float64}, output_ϵ_raw)


    # compute hessian-vector products
    num_products = size(R,2) # number of hessian-vector products
    num_chunks = div(num_products, CHUNK)
    @assert size(R,1) == length(local_to_global_idx)
    numVar = length(x)

    for k in 1:CHUNK:CHUNK*num_chunks

        for r in 1:length(local_to_global_idx)
            # set up directional derivatives
            @inbounds idx = local_to_global_idx[r]
            # load up R[r,k,k+1,...,k+CHUNK-1] into input_ϵ
            for s in 1:CHUNK
                input_ϵ_raw[(idx-1)*CHUNK + s] = R[r,k+s-1]
            end
            @inbounds output_ϵ[idx] = zero_ϵ
        end

        hessian_slice_inner(d, ex, R, input_ϵ, output_ϵ, Val{CHUNK})

        # collect directional derivatives
        for r in 1:length(local_to_global_idx)
            idx = local_to_global_idx[r]
            # load output_ϵ into R[r,k,k+1,...,k+CHUNK-1]
            for s in 1:CHUNK
                R[r,k+s-1] = output_ϵ_raw[(idx-1)*CHUNK + s]
            end
            @inbounds input_ϵ[idx] = zero_ϵ
        end

    end

    # leftover chunk
    remaining = num_products - CHUNK*num_chunks
    if remaining > 0
        k = CHUNK*num_chunks+1
        for r in 1:length(local_to_global_idx)
            # set up directional derivatives
            @inbounds idx = local_to_global_idx[r]
            # load up R[r,k,k+1,...,k+remaining-1] into input_ϵ
            for s in 1:remaining
                # leave junk in the unused components
                input_ϵ_raw[(idx-1)*CHUNK + s] = R[r,k+s-1]
            end
            @inbounds output_ϵ[idx] = zero_ϵ
        end

        hessian_slice_inner(d, ex, R, input_ϵ, output_ϵ, Val{CHUNK})

        # collect directional derivatives
        for r in 1:length(local_to_global_idx)
            idx = local_to_global_idx[r]
            # load output_ϵ into R[r,k,k+1,...,k+remaining-1]
            for s in 1:remaining
                R[r,k+s-1] = output_ϵ_raw[(idx-1)*CHUNK + s]
            end
            @inbounds input_ϵ[idx] = zero_ϵ
        end
    end

    # Output is in R, now recover

    #output_slice = view(H, (nzcount+1):(nzcount+nzthis))
    output_slice = VectorView(nzcount, nzthis, pointer(H))
    Coloring.recover_from_matmat!(output_slice, R, ex.rinfo, recovery_tmp_storage)
    rmul!(output_slice, scale)
    return nzthis

end

MathProgBase.isobjlinear(d::NLPEvaluator) = !d.has_nlobj && (length(d.m.obj.qvars1) == 0)
# interpret quadratic to include purely linear
MathProgBase.isobjquadratic(d::NLPEvaluator) = !d.has_nlobj

MathProgBase.isconstrlinear(d::NLPEvaluator, i::Integer) = (i <= length(d.m.linconstr))

function MathProgBase.jac_structure(d::NLPEvaluator)
    # Jacobian structure
    jac_I = Int[]
    jac_J = Int[]
    A = d.A
    for col = 1:size(A,2)
        for pos = nzrange(A,col)
            push!(jac_I, A.rowval[pos])
            push!(jac_J, col)
        end
    end
    rowoffset = size(A,1)+1
    for c::QuadConstraint in d.m.quadconstr
        aff = c.terms.aff
        for k in 1:length(aff.vars)
            push!(jac_I, rowoffset)
            push!(jac_J, aff.vars[k].col)
        end
        for k in 1:length(c.terms.qvars1)
            push!(jac_I, rowoffset)
            push!(jac_I, rowoffset)
            push!(jac_J, c.terms.qvars1[k].col)
            push!(jac_J, c.terms.qvars2[k].col)
        end
        rowoffset += 1
    end
    for ex in d.constraints
        idx = ex.grad_sparsity
        for i in 1:length(idx)
            push!(jac_I, rowoffset)
            push!(jac_J, idx[i])
        end
        rowoffset += 1
    end
    return jac_I, jac_J
end
function MathProgBase.hesslag_structure(d::NLPEvaluator)
    d.want_hess || error("Hessian computations were not requested on the call to MathProgBase.initialize.")
    return d.hess_I,d.hess_J
end
function _hesslag_structure(d::NLPEvaluator)
    hess_I = Int[]
    hess_J = Int[]

    qobj::QuadExpr = d.m.obj
    for k in 1:length(qobj.qvars1)
        qidx1 = qobj.qvars1[k].col
        qidx2 = qobj.qvars2[k].col
        if qidx2 > qidx1
            qidx1, qidx2 = qidx2, qidx1
        end
        push!(hess_I, qidx1)
        push!(hess_J, qidx2)
    end
    # quadratic constraints
    for c::QuadConstraint in d.m.quadconstr
        for k in 1:length(c.terms.qvars1)
            qidx1 = c.terms.qvars1[k].col
            qidx2 = c.terms.qvars2[k].col
            if qidx2 > qidx1
                qidx1, qidx2 = qidx2, qidx1
            end
            push!(hess_I, qidx1)
            push!(hess_J, qidx2)
        end
    end

    if d.has_nlobj
        append!(hess_I, d.objective.hess_I)
        append!(hess_J, d.objective.hess_J)
    end
    for ex in d.constraints
        append!(hess_I, ex.hess_I)
        append!(hess_J, ex.hess_J)
    end

    return hess_I, hess_J
end

# currently don't merge duplicates (this isn't required by MPB standard)
function affToExpr(aff::AffExpr, constant::Bool)
    if length(aff.vars) == 0 && !constant
        return 0
    end
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

mutable struct VariablePrintWrapper
    v::Variable
    mode
end
Base.show(io::IO,v::VariablePrintWrapper) = print(io,var_str(v.mode,v.v))
mutable struct ParameterPrintWrapper
    idx::Int
    mode
end
function Base.show(io::IO,p::ParameterPrintWrapper)
    if p.mode == IJuliaMode
        print(io,"parameter_{$(p.idx)}")
    else
        print(io,"parameter[$(p.idx)]")
    end
end
mutable struct SubexpressionPrintWrapper
    idx::Int
    mode
end
function Base.show(io::IO,s::SubexpressionPrintWrapper)
    if s.mode == IJuliaMode
        print(io,"subexpression_{$(s.idx)}")
    else
        print(io,"subexpression[$(s.idx)]")
    end
end

# we splat in the subexpressions (for now)
function tapeToExpr(m::Model, k, nd::Vector{NodeData}, adj, const_values, parameter_values, subexpressions::Vector{Any},user_operators::ReverseDiffSparse.UserOperatorRegistry, generic_variable_names::Bool, splat_subexpressions::Bool, print_mode=REPLMode)

    children_arr = rowvals(adj)

    nod = nd[k]
    if nod.nodetype == VARIABLE
        if generic_variable_names
            return Expr(:ref,:x,nod.index)
        else
            # mode only matters when generic_variable_names == false
            return VariablePrintWrapper(Variable(m,nod.index),print_mode)
        end
    elseif nod.nodetype == VALUE
        return const_values[nod.index]
    elseif nod.nodetype == SUBEXPRESSION
        if splat_subexpressions
            return subexpressions[nod.index]
        else
            return SubexpressionPrintWrapper(nod.index,print_mode)
        end
    elseif nod.nodetype == PARAMETER
        if splat_subexpressions
            return parameter_values[nod.index]
        else
            return ParameterPrintWrapper(nod.index,print_mode)
        end
    elseif nod.nodetype == CALL
        op = nod.index
        opsymbol = :error
        if op < ReverseDiffSparse.USER_OPERATOR_ID_START
            opsymbol = operators[op]
        else
            for (key,value) in user_operators.multivariate_operator_to_id
                if value == op - ReverseDiffSparse.USER_OPERATOR_ID_START + 1
                    opsymbol = key
                end
            end
        end
        @assert opsymbol != :error
        children_idx = nzrange(adj,k)
        if opsymbol == :+ && length(children_idx) == 0
            return 0
        elseif opsymbol == :* && length(children_idx) == 0
            return 1
        end
        ex = Expr(:call,opsymbol)
        for cidx in children_idx
            push!(ex.args, tapeToExpr(m, children_arr[cidx], nd, adj, const_values, parameter_values, subexpressions, user_operators, generic_variable_names, splat_subexpressions, print_mode))
        end
        return ex
    elseif nod.nodetype == CALLUNIVAR
        op = nod.index
        opsymbol = :error
        if op < ReverseDiffSparse.USER_UNIVAR_OPERATOR_ID_START
            opsymbol = univariate_operators[op]
        else
            for (key,value) in user_operators.univariate_operator_to_id
                if value == op - ReverseDiffSparse.USER_UNIVAR_OPERATOR_ID_START + 1
                    opsymbol = key
                end
            end
        end
        @assert opsymbol != :error
        cidx = first(nzrange(adj,k))
        return Expr(:call,opsymbol,tapeToExpr(m, children_arr[cidx], nd, adj, const_values, parameter_values, subexpressions, user_operators, generic_variable_names, splat_subexpressions, print_mode))
    elseif nod.nodetype == COMPARISON
        op = nod.index
        opsymbol = comparison_operators[op]
        children_idx = nzrange(adj,k)
        if length(children_idx) > 2
            ex = Expr(:comparison)
            for cidx in children_idx
                push!(ex.args, tapeToExpr(m, children_arr[cidx], nd, adj, const_values, parameter_values, subexpressions, user_operators, generic_variable_names, splat_subexpressions, print_mode))
                push!(ex.args, opsymbol)
            end
            pop!(ex.args)
        else
            ex = Expr(:call, opsymbol)
            push!(ex.args, tapeToExpr(m, children_arr[children_idx[1]], nd, adj, const_values, parameter_values, subexpressions, user_operators, generic_variable_names, splat_subexpressions, print_mode))
            push!(ex.args, tapeToExpr(m, children_arr[children_idx[2]], nd, adj, const_values, parameter_values, subexpressions, user_operators, generic_variable_names, splat_subexpressions, print_mode))
        end
        return ex
    elseif nod.nodetype == LOGIC
        op = nod.index
        opsymbol = logic_operators[op]
        children_idx = nzrange(adj,k)
        lhs = tapeToExpr(m, children_arr[first(children_idx)], nd, adj, const_values, parameter_values, subexpressions, user_operators, generic_variable_names, splat_subexpressions, print_mode)
        rhs = tapeToExpr(m, children_arr[last(children_idx)], nd, adj, const_values, parameter_values, subexpressions, user_operators, generic_variable_names, splat_subexpressions, print_mode)
        return Expr(opsymbol, lhs, rhs)
    end
    error()


end


function MathProgBase.obj_expr(d::NLPEvaluator)
    if d.has_nlobj
        # expressions are simplified if requested
        ex = d.objective
        return tapeToExpr(d.m, 1, ex.nd, ex.adj, ex.const_values, d.parameter_values, d.subexpressions_as_julia_expressions,d.m.nlpdata.user_operators, true, true)
    else
        return quadToExpr(d.m.obj, true)
    end
end

function MathProgBase.constr_expr(d::NLPEvaluator,i::Integer)
    nlin = length(d.m.linconstr)
    nquad = length(d.m.quadconstr)
    if i <= nlin
        constr = d.m.linconstr[i]
        ex = affToExpr(constr.terms, false)
        if sense(constr) == :range
            return Expr(:comparison, constr.lb, :(<=), ex, :(<=), constr.ub)
        else
            return Expr(:call, sense(constr), ex, rhs(constr))
        end
    elseif i > nlin && i <= nlin + nquad
        i -= nlin
        qconstr = d.m.quadconstr[i]
        return Expr(:call, qconstr.sense, quadToExpr(qconstr.terms, true), 0)
    else
        i -= nlin + nquad
        ex = d.constraints[i] # may be simplified
        julia_expr = tapeToExpr(d.m, 1, ex.nd, ex.adj, ex.const_values, d.parameter_values, d.subexpressions_as_julia_expressions,d.m.nlpdata.user_operators, true, true)
        constr = d.m.nlpdata.nlconstr[i]
        if sense(constr) == :range
            return Expr(:comparison, constr.lb, :(<=), julia_expr, :(<=), constr.ub)
        else
            return Expr(:call, sense(constr), julia_expr, rhs(constr))
        end
    end
end

function EnableNLPResolve()
    warn_once("NLP resolve is now enabled by default. The EnableNLPResolve() method will be removed in a future release.")
    return
end
function DisableNLPResolve()
    warn_once("NLP resolve is now enabled by default. The DisableNLPResolve() method will be removed in a future release.")
    return
end
export EnableNLPResolve, DisableNLPResolve

function _buildInternalModel_nlp(m::Model, traits)

    linobj = prepAffObjective(m)
    linrowlb, linrowub = prepConstrBounds(m)

    nldata::NLPData = m.nlpdata
    if m.internalModelLoaded && !m.simplify_nonlinear_expressions
        @assert isa(nldata.evaluator, NLPEvaluator)
        d = nldata.evaluator
        fill!(d.last_x, NaN)
    else
        d = NLPEvaluator(m)
        nldata.evaluator = d
    end

    nlp_lb, nlp_ub = constraintbounds(m)
    numConstr = length(nlp_lb)

    m.internalModel = MathProgBase.NonlinearModel(m.solver)

    MathProgBase.loadproblem!(m.internalModel, m.numCols, numConstr, m.colLower, m.colUpper, nlp_lb, nlp_ub, m.objSense, d)
    if traits.int
        if applicable(MathProgBase.setvartype!, m.internalModel, m.colCat)
            MathProgBase.setvartype!(m.internalModel, vartypes_without_fixed(m))
        else
            error("Solver does not support discrete variables")
        end
    end

    MathProgBase.setwarmstart!(m.internalModel, tidy_warmstart(m))

    registercallbacks(m)

    m.internalModelLoaded = true

    nothing
end


function solvenlp(m::Model, traits; suppress_warnings=false)

    @assert m.internalModelLoaded

    MathProgBase.optimize!(m.internalModel)
    stat = MathProgBase.status(m.internalModel)

    if stat != :Infeasible && stat != :Unbounded
        m.objVal = MathProgBase.getobjval(m.internalModel)
        m.colVal = MathProgBase.getsolution(m.internalModel)
        try
            objBound = MathProgBase.getobjbound(m.internalModel) + m.obj.aff.constant
            # Don't corrupt objBound if the above call fails
            m.objBound = objBound
        catch
            nothing
        end
    end

    if stat != :Optimal
        suppress_warnings || Compat.@warn("Not solved to optimality, status: $stat")
    end
    if stat == :Optimal && !traits.int
        if applicable(MathProgBase.getconstrduals, m.internalModel) && applicable(MathProgBase.getreducedcosts, m.internalModel)
            nlduals = MathProgBase.getconstrduals(m.internalModel)
            m.linconstrDuals = nlduals[1:length(m.linconstr)]
            # quadratic duals currently not available, formulate as nonlinear constraint if needed
            m.nlpdata.nlconstrDuals = nlduals[length(m.linconstr)+length(m.quadconstr)+1:end]
            m.redCosts = MathProgBase.getreducedcosts(m.internalModel)
        else
            suppress_warnings || warn_once("Nonlinear solver does not provide dual solutions")
        end
    end

    #d = m.nlpdata.evaluator
    #println("feval $(d.eval_f_timer)\nfgrad $(d.eval_grad_f_timer)\ngeval $(d.eval_g_timer)\njaceval $(d.eval_jac_g_timer)\nhess $(d.eval_hesslag_timer)")

    return stat::Symbol

end

# getvalue for nonlinear subexpressions
getvalue(x::NonlinearExpression) = _getValue(x)
function _getValue(x::NonlinearExpression)
    m = x.m
    # recompute EVERYTHING here
    # could be smarter and cache

    nldata::NLPData = m.nlpdata
    subexpr = Array{Vector{NodeData}}(undef, 0)
    for nlexpr in nldata.nlexpr
        push!(subexpr, nlexpr.nd)
    end

    this_subexpr = nldata.nlexpr[x.index]

    max_len = length(this_subexpr.nd)

    subexpression_order, individual_order = order_subexpressions(Vector{NodeData}[this_subexpr.nd],subexpr)

    subexpr_values = Array{Float64}(undef, length(subexpr))

    for k in subexpression_order
        max_len = max(max_len, length(nldata.nlexpr[k].nd))
    end

    forward_storage = Array{Float64}(undef, max_len)
    partials_storage = Array{Float64}(undef, max_len)
    user_input_buffer = zeros(nldata.largest_user_input_dimension)
    user_output_buffer = zeros(nldata.largest_user_input_dimension)

    for k in subexpression_order # compute value of dependent subexpressions
        ex = nldata.nlexpr[k]
        adj = adjmat(ex.nd)
        subexpr_values[k] = forward_eval(forward_storage,partials_storage,ex.nd,adj,ex.const_values,nldata.nlparamvalues,m.colVal,subexpr_values,user_input_buffer,user_output_buffer)
    end

    adj = adjmat(this_subexpr.nd)

    return forward_eval(forward_storage,partials_storage,this_subexpr.nd,adj,this_subexpr.const_values,nldata.nlparamvalues,m.colVal,subexpr_values,user_input_buffer,user_output_buffer)
end

mutable struct UserFunctionEvaluator <: MathProgBase.AbstractNLPEvaluator
    f
    ∇f
    len::Int
end

function MathProgBase.eval_f(d::UserFunctionEvaluator,x)
    @assert length(x) == d.len
    return d.f(x)::eltype(x)
end
function MathProgBase.eval_grad_f(d::UserFunctionEvaluator,grad,x)
    d.∇f(grad,x)
    nothing
end

function UserAutoDiffEvaluator(dimension::Integer, f::Function, ::Type{T} = Float64) where T
    g = x -> f(x...)
    cfg = ForwardDiff.GradientConfig(g, zeros(T, dimension))
    ∇f = (out, y) -> ForwardDiff.gradient!(out, g, y, cfg)
    return UserFunctionEvaluator(g, ∇f, dimension)
end

function register(m::Model, s::Symbol, dimension::Integer, f::Function; autodiff::Bool=false)
    autodiff == true || error("If only the function is provided, must set autodiff=true")
    initNLP(m)

    if dimension == 1
        fprime = x -> ForwardDiff.derivative(f, x)
        fprimeprime = x -> ForwardDiff.derivative(fprime, x)
        ReverseDiffSparse.register_univariate_operator!(m.nlpdata.user_operators, s, f, fprime, fprimeprime)
    else
        m.nlpdata.largest_user_input_dimension = max(m.nlpdata.largest_user_input_dimension,dimension)
        ReverseDiffSparse.register_multivariate_operator!(m.nlpdata.user_operators, s, UserAutoDiffEvaluator(dimension, f))
    end

end

function register(m::Model, s::Symbol, dimension::Integer, f::Function, ∇f::Function; autodiff::Bool=false)
    initNLP(m)
    if dimension == 1
        autodiff == true || error("Currently must provide 2nd order derivatives of univariate functions. Try setting autodiff=true.")
        fprimeprime = x -> ForwardDiff.derivative(∇f, x)
        ReverseDiffSparse.register_univariate_operator!(m.nlpdata.user_operators, s, f, ∇f, fprimeprime)
    else
        autodiff == false || warn_once("autodiff=true ignored since gradient is already provided.")
        m.nlpdata.largest_user_input_dimension = max(m.nlpdata.largest_user_input_dimension,dimension)
        d = UserFunctionEvaluator(x -> f(x...), (g,x)->∇f(g,x...), dimension)
        ReverseDiffSparse.register_multivariate_operator!(m.nlpdata.user_operators, s, d)
    end

end

function register(m::Model, s::Symbol, dimension::Integer, f::Function, ∇f::Function, ∇²f::Function)
    dimension == 1 || error("Providing hessians for multivariate functions is not yet supported")
    initNLP(m)
    ReverseDiffSparse.register_univariate_operator!(m.nlpdata.user_operators, s, f, ∇f, ∇²f)
end

register(s::Symbol, dimension::Integer, f::Function; autodiff::Bool=false) = error("Function registration is now local to JuMP models. JuMP.register requires the model object as the first argument.")
register(s::Symbol, dimension::Integer, f::Function, ∇f::Function; autodiff::Bool=false) = error("Function registration is now local to JuMP models. JuMP.register requires the model object as the first argument.")
register(s::Symbol, dimension::Integer, f::Function, ∇f::Function, ∇²f::Function) = error("Function registration is now local to JuMP models. JuMP.register requires the model object as the first argument.")


# Ex: setNLobjective(m, :Min, :($x + $y^2))
function setNLobjective(m::Model, sense::Symbol, x)
    initNLP(m)
    setobjectivesense(m, sense)
    ex = NonlinearExprData(m, x)
    m.nlpdata.nlobj = ex
    m.obj = zero(QuadExpr)
    m.internalModelLoaded = false
    return
end

# Ex: addNLconstraint(m, :($x + $y^2 <= 1))
function addNLconstraint(m::Model, ex::Expr)
    initNLP(m)
    m.internalModelLoaded = false
    if isexpr(ex, :call) # one-sided constraint
        # Simple comparison - move everything to the LHS
        op = ex.args[1]
        if op == :(==)
            lb = 0.0
            ub = 0.0
        elseif op == :(<=) || op == :(≤)
            lb = -Inf
            ub = 0.0
        elseif op == :(>=) || op == :(≥)
            lb = 0.0
            ub = Inf
        else
            error("in addNLconstraint ($ex): expected comparison operator (<=, >=, or ==).")
        end
        lhs = :($(ex.args[2]) - $(ex.args[3]))
        c = NonlinearConstraint(NonlinearExprData(m, lhs), lb, ub)
        push!(m.nlpdata.nlconstr, c)
        return ConstraintRef{Model,NonlinearConstraint}(m, length(m.nlpdata.nlconstr))
    elseif isexpr(ex, :comparison)
        # ranged row
        if (ex.args[2] != :<= && ex.args[2] != :≤) || (ex.args[4] != :<= && ex.args[4] != :≤)
            error("in addNLconstraint ($ex): only ranged rows of the form lb <= expr <= ub are supported.")
        end
        lb = ex.args[1]
        ub = ex.args[5]
        if !isa(lb,Number)
            error(string("in addNLconstraint (",ex,"): expected ",lb," to be a number."))
        elseif !isa(ub,Number)
            error(string("in addNLconstraint (",ex,"): expected ",ub," to be a number."))
        end
        c = NonlinearConstraint(NonlinearExprData(m, ex.args[3]), lb, ub)
        push!(m.nlpdata.nlconstr, c)
        return ConstraintRef{Model,NonlinearConstraint}(m, length(m.nlpdata.nlconstr))
    else
        # Unknown
        error("in addNLconstraint ($ex): constraints must be in one of the following forms:\n" *
              "       expr1 <= expr2\n" * "       expr1 >= expr2\n" *
              "       expr1 == expr2")
    end
end

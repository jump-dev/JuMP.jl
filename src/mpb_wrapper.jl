
import MathProgBase
# Implements a MathProgBase solver using the expression graph interface


immutable RDSSolver <: MathProgBase.AbstractMathProgSolver
    realsolver::MathProgBase.AbstractMathProgSolver
end

export RDSSolver

type RDSModel <: MathProgBase.AbstractNonlinearModel
    realmodel::MathProgBase.AbstractNonlinearModel
    numVar::Int
    numConstr::Int
    nlp_eval
end

MathProgBase.NonlinearModel(s::RDSSolver) = RDSModel(MathProgBase.NonlinearModel(s.realsolver), 0, 0, nothing)

type FunctionStorage
    nd::Vector{NodeData}
    adj::SparseMatrixCSC{Bool,Int}
    const_values::Vector{Float64}
    forward_storage::Vector{Float64}
    reverse_storage::Vector{Float64}
    grad_sparsity::Vector{Int}
    hess_I::Vector{Int} # nonzero pattern of hessian
    hess_J::Vector{Int}
    rinfo::Coloring.RecoveryInfo # coloring info for hessians
    seed_matrix::Matrix{Float64}
    linearity::Linearity
end

type RDSNLPEvaluator <: MathProgBase.AbstractNLPEvaluator
    numVar::Int
    numConstr::Int
    nlp_eval
    expressions::Vector{FunctionStorage} # objective is in the first slot
    last_x::Vector{Float64}
    jac_storage::Vector{Float64} # temporary storage for computing jacobians
    # storage for computing hessians
    forward_storage_hess::Vector{Dual{Float64}} # length is of the longest expression
    reverse_storage_hess::Vector{Dual{Float64}} # length is of the longest expression
    forward_input_vector::Vector{Dual{Float64}} # length is number of variables
    reverse_output_vector::Vector{Dual{Float64}}# length is number of variables
end

function MathProgBase.loadproblem!(m::RDSModel, numVar, numConstr, l, u, lb, ub, sense, d::MathProgBase.AbstractNLPEvaluator)

    m.numVar = numVar
    m.numConstr = numConstr

    nlp = RDSNLPEvaluator(numVar, numConstr, d, FunctionStorage[], Array(Float64,numVar), Array(Float64, numVar),Vector{Dual{Float64}}(0),Vector{Dual{Float64}}(0),
    Array(Dual{Float64},numVar),Array(Dual{Float64},numVar))

    MathProgBase.loadproblem!(m.realmodel, numVar, numConstr, l, u, lb, ub, sense, nlp)
end


MathProgBase.features_available(d::RDSNLPEvaluator) = [:Grad,:Jac,:Hess]

# convert expression into flat tree format
function FunctionStorage(ex,numVar, want_hess::Bool)
    
    nd,const_values = expr_to_nodedata(ex)
    adj = adjmat(nd)
    forward_storage = zeros(length(nd))
    reverse_storage = zeros(length(nd))
    grad_sparsity = compute_gradient_sparsity(nd, adj)

    if want_hess
        # compute hessian sparsity
        linearity = classify_linearity(nd, adj)
        edgelist = compute_hessian_sparsity(nd, adj, linearity)
        hess_I, hess_J, rinfo = Coloring.hessian_color_preprocess(edgelist, numVar)
        seed_matrix = Coloring.seed_matrix(rinfo)
    else
        hess_I = hess_J = Int[]
        rinfo = Coloring.RecoveryInfo()
        seed_matrix = Array(Float64,0,0)
        linearity = [NONLINEAR]
    end

    return FunctionStorage(nd, adj, const_values, forward_storage, reverse_storage, grad_sparsity, hess_I, hess_J, rinfo, seed_matrix, linearity[1])

end

function MathProgBase.initialize(d::RDSNLPEvaluator, requested_features)

    MathProgBase.initialize(d.nlp_eval, [:ExprGraph]) # Jac and Grad for debugging only

    want_hess = (:Hess in requested_features)

    obj = FunctionStorage(MathProgBase.obj_expr(d.nlp_eval), d.numVar, want_hess)
    #@show MathProgBase.obj_expr(d.nlp_eval)
    max_expr_length = length(obj.nd)
    push!(d.expressions, obj)

    for i in 1:d.numConstr
        constr = MathProgBase.constr_expr(d.nlp_eval,i)
        @assert constr.head == :comparison
        if length(constr.args) == 3
            ex = constr.args[1]
        else
            ex = constr.args[3]
        end
        push!(d.expressions, FunctionStorage(ex, d.numVar, want_hess))
        max_expr_length = max(max_expr_length, length(d.expressions[end].nd))
    end

    if want_hess # allocate extra storage
        d.forward_storage_hess = Array(Dual{Float64},max_expr_length)
        d.reverse_storage_hess = Array(Dual{Float64},max_expr_length)
    end
end

function forward_eval_all(d::RDSNLPEvaluator,x)
    # do a forward pass on all expressions at x
    for i in 1:d.numConstr+1
        ex = d.expressions[i]
        forward_eval(ex.forward_storage,ex.nd,ex.adj,ex.const_values,x)
    end
    copy!(d.last_x,x)
end

function MathProgBase.eval_f(d::RDSNLPEvaluator,x)
    if d.last_x != x
        forward_eval_all(d,x)
    end
    fval = d.expressions[1].forward_storage[1]

    #=
    fval_j = MathProgBase.eval_f(d.nlp_eval,x)
    if !isapprox(fval,fval_j)
        @show x
        @show d.last_x
        @show fval
        @show fval_j
        exit()
    end=#
    return fval
end

function MathProgBase.eval_g(d::RDSNLPEvaluator,g,x)
    if d.last_x != x
        forward_eval_all(d,x)
    end

    for i in 1:d.numConstr
        ex = d.expressions[i+1]
        g[i] = ex.forward_storage[1]
    end
    return
end

function MathProgBase.eval_grad_f(d::RDSNLPEvaluator,g,x)
    if d.last_x != x
        forward_eval_all(d,x)
    end

    fill!(g,0.0)
    ex = d.expressions[1]
    reverse_eval(g,ex.reverse_storage,ex.forward_storage,ex.nd,ex.adj,ex.const_values)
    return
end

function MathProgBase.jac_structure(d::RDSNLPEvaluator)
    I = Int[]
    J = Int[]

    for i in 1:d.numConstr
        idx = d.expressions[i+1].grad_sparsity
        for k in 1:length(idx)
            push!(I,i)
            push!(J,idx[k])
        end
    end

    return I,J
end

function MathProgBase.eval_jac_g(d::RDSNLPEvaluator, J, x)
    if d.last_x != x
        forward_eval_all(d,x)
    end

    grad_storage = d.jac_storage
    nz = 0

    for i in 1:d.numConstr
        ex = d.expressions[i+1]
        idx = ex.grad_sparsity
        grad_storage[idx] = 0.0
        reverse_eval(grad_storage,ex.reverse_storage,ex.forward_storage,ex.nd,ex.adj,ex.const_values)
        for k in 1:length(idx)
            J[nz+k] = grad_storage[idx[k]]
        end
        nz += length(idx)
    end

    return
end

function MathProgBase.hesslag_structure(d::RDSNLPEvaluator)

    hess_I = Int[]
    hess_J = Int[]

    for ex in d.expressions
        append!(hess_I, ex.hess_I)
        append!(hess_J, ex.hess_J)
    end
    return hess_I, hess_J
end

function MathProgBase.eval_hesslag(d::RDSNLPEvaluator, H, x, σ, μ)

    nzcount = 0

    for i in 1:d.numVar
        d.forward_input_vector[i] = Dual(x[i],0.0)
    end

    recovery_tmp_storage = reinterpret(Float64, d.reverse_output_vector)
    for i in 1:length(d.expressions)
        ex = d.expressions[i]
        nzthis = length(ex.hess_I)
        if ex.linearity == LINEAR
            @assert nzthis == 0
            continue
        end
        seed = ex.seed_matrix
        Coloring.prepare_seed_matrix!(seed,ex.rinfo)


        hessmat_eval!(seed, d.reverse_storage_hess, d.forward_storage_hess, ex.nd, ex.adj, ex.const_values, x, d.reverse_output_vector, d.forward_input_vector, ex.rinfo.local_indices)
        # Output is in seed, now recover

        output_slice = sub(H, (nzcount+1):(nzcount+nzthis))
        Coloring.recover_from_matmat!(output_slice, seed, ex.rinfo, recovery_tmp_storage)
        if i == 1
            scale!(output_slice, σ)
        else
            scale!(output_slice, μ[i-1])
        end
        nzcount += nzthis
    end

end

MathProgBase.setwarmstart!(m::RDSModel, x) = MathProgBase.setwarmstart!(m.realmodel,x)
MathProgBase.optimize!(m::RDSModel) = MathProgBase.optimize!(m.realmodel)
MathProgBase.status(m::RDSModel) = MathProgBase.status(m.realmodel)
MathProgBase.getsolution(m::RDSModel) = MathProgBase.getsolution(m.realmodel)
MathProgBase.getobjval(m::RDSModel) = MathProgBase.getobjval(m.realmodel)
MathProgBase.getreducedcosts(m::RDSModel) = MathProgBase.getreducedcosts(m.realmodel)
MathProgBase.getconstrduals(m::RDSModel) = MathProgBase.getconstrduals(m.realmodel)

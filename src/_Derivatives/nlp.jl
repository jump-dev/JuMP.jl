"""
    _replace_moi_variables(
        nd::Vector{NodeData},
        moi_index_to_consecutive_index,
    )

This function returns a copy of `nd::Vector{NodeData}`, but with the
`MOIVARIABLE` nodes replaced by a new `MOI.VariableIndex` in which the index of
the variable is the 1-indexed column of the variable in the problem.
"""
function _replace_moi_variables(
    nd::Vector{NodeData},
    moi_index_to_consecutive_index::Dict{MOI.VariableIndex,Int},
)
    new_nd = Vector{NodeData}(undef, length(nd))
    for (i, node) in enumerate(nd)
        if node.nodetype == MOIVARIABLE
            new_nd[i] = NodeData(
                VARIABLE,
                moi_index_to_consecutive_index[MOI.VariableIndex(node.index)],
                node.parent,
            )
        else
            new_nd[i] = node
        end
    end
    return new_nd
end

mutable struct _SubexpressionStorage
    nd::Vector{NodeData}
    adj::SparseArrays.SparseMatrixCSC{Bool,Int}
    const_values::Vector{Float64}
    forward_storage::Vector{Float64}
    partials_storage::Vector{Float64}
    reverse_storage::Vector{Float64}
    forward_storage_ϵ::Vector{Float64}
    partials_storage_ϵ::Vector{Float64}
    reverse_storage_ϵ::Vector{Float64}
    linearity::Linearity
end

function _SubexpressionStorage(
    nd::Vector{NodeData},
    const_values::Vector{Float64},
    num_variables,
    subexpression_linearity,
    moi_index_to_consecutive_index::Dict{MOI.VariableIndex,Int},
)
    nd = _replace_moi_variables(nd, moi_index_to_consecutive_index)
    adj = adjmat(nd)
    # classify_linearity returns a vector containing the linearity of each
    # element in `nd`. For the linearity of the subexpression, we only care
    # about the linearity of the first element.
    linearity = classify_linearity(nd, adj, subexpression_linearity)
    return _SubexpressionStorage(
        nd,
        adj,
        const_values,
        zeros(length(nd)),  # forward_storage
        zeros(length(nd)),  # partials_storage
        zeros(length(nd)),  # reverse_storage
        Float64[],  # forward_storage_ϵ
        Float64[],  # partials_storage_ϵ
        Float64[],  # reverse_storage_ϵ
        linearity[1],
    )
end

mutable struct _FunctionStorage
    nd::Vector{NodeData}
    adj::SparseArrays.SparseMatrixCSC{Bool,Int}
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

function _FunctionStorage(
    nd::Vector{NodeData},
    const_values,
    num_variables,
    coloring_storage,
    want_hess::Bool,
    subexpressions::Vector{_SubexpressionStorage},
    dependent_subexpressions,
    subexpression_linearity,
    subexpression_edgelist,
    subexpression_variables,
    moi_index_to_consecutive_index::Dict{MOI.VariableIndex,Int},
)
    nd = _replace_moi_variables(nd, moi_index_to_consecutive_index)
    adj = adjmat(nd)
    forward_storage = zeros(length(nd))
    partials_storage = zeros(length(nd))
    reverse_storage = zeros(length(nd))
    empty!(coloring_storage)
    compute_gradient_sparsity!(coloring_storage, nd)
    for k in dependent_subexpressions
        compute_gradient_sparsity!(coloring_storage, subexpressions[k].nd)
    end
    grad_sparsity = sort!(collect(coloring_storage))
    empty!(coloring_storage)
    if want_hess
        # compute hessian sparsity
        linearity = classify_linearity(nd, adj, subexpression_linearity)
        edgelist = compute_hessian_sparsity(
            nd,
            adj,
            linearity,
            coloring_storage,
            subexpression_edgelist,
            subexpression_variables,
        )
        hess_I, hess_J, rinfo = Coloring.hessian_color_preprocess(
            edgelist,
            num_variables,
            coloring_storage,
        )
        seed_matrix = Coloring.seed_matrix(rinfo)
    else
        hess_I = hess_J = Int[]
        rinfo = Coloring.RecoveryInfo()
        seed_matrix = Array{Float64}(undef, 0, 0)
        linearity = [NONLINEAR]
    end
    return _FunctionStorage(
        nd,
        adj,
        const_values,
        forward_storage,
        partials_storage,
        reverse_storage,
        grad_sparsity,
        hess_I,
        hess_J,
        rinfo,
        seed_matrix,
        linearity[1],
        dependent_subexpressions,
    )
end

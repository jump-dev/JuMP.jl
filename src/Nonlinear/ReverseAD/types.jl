#  Copyright 2017, Iain Dunning, Joey Huchette, Miles Lubin, and contributors
#  This Source Code Form is subject to the terms of the Mozilla Public
#  License, v. 2.0. If a copy of the MPL was not distributed with this
#  file, You can obtain one at https://mozilla.org/MPL/2.0/.

mutable struct _SubexpressionStorage
    nodes::Vector{Nonlinear.Node}
    adj::SparseArrays.SparseMatrixCSC{Bool,Int}
    const_values::Vector{Float64}
    forward_storage::Vector{Float64}
    partials_storage::Vector{Float64}
    reverse_storage::Vector{Float64}
    forward_storage_ϵ::Vector{Float64}
    partials_storage_ϵ::Vector{Float64}
    reverse_storage_ϵ::Vector{Float64}
    linearity::Linearity

    function _SubexpressionStorage(
        expr::Nonlinear.NonlinearExpression,
        subexpression_linearity,
        moi_index_to_consecutive_index,
    )
        nodes =
            _replace_moi_variables(expr.nodes, moi_index_to_consecutive_index)
        adj = Nonlinear._adjacency_matrix(nodes)
        N = length(nodes)
        linearity = _classify_linearity(nodes, adj, subexpression_linearity)
        return new(
            nodes,
            adj,
            expr.values,
            zeros(N),  # forward_storage,
            zeros(N),  # partials_storage,
            zeros(N),  # reverse_storage,
            Float64[],
            Float64[],
            Float64[],
            linearity[1],
        )
    end
end

mutable struct _FunctionStorage
    nodes::Vector{Nonlinear.Node}
    adj::SparseArrays.SparseMatrixCSC{Bool,Int}
    const_values::Vector{Float64}
    forward_storage::Vector{Float64}
    partials_storage::Vector{Float64}
    reverse_storage::Vector{Float64}
    grad_sparsity::Vector{Int}
    # Nonzero pattern of Hessian matrix
    hess_I::Vector{Int}
    hess_J::Vector{Int}
    rinfo::Coloring.RecoveryInfo # coloring info for hessians
    seed_matrix::Matrix{Float64}
    linearity::Linearity
    # subexpressions which this function depends on, ordered for forward pass.
    dependent_subexpressions::Vector{Int}

    function _FunctionStorage(
        nodes::Vector{Nonlinear.Node},
        const_values,
        num_variables,
        coloring_storage::Coloring.IndexedSet,
        want_hess::Bool,
        subexpressions::Vector{_SubexpressionStorage},
        dependent_subexpressions,
        subexpression_linearity,
        subexpression_edgelist,
        subexpression_variables,
        moi_index_to_consecutive_index,
    )
        nodes = _replace_moi_variables(nodes, moi_index_to_consecutive_index)
        adj = Nonlinear._adjacency_matrix(nodes)
        N = length(nodes)
        empty!(coloring_storage)
        _compute_gradient_sparsity!(coloring_storage, nodes)
        for k in dependent_subexpressions
            _compute_gradient_sparsity!(
                coloring_storage,
                subexpressions[k].nodes,
            )
        end
        grad_sparsity = sort!(collect(coloring_storage))
        empty!(coloring_storage)
        if want_hess
            linearity = _classify_linearity(nodes, adj, subexpression_linearity)
            edgelist = _compute_hessian_sparsity(
                nodes,
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
            return new(
                nodes,
                adj,
                const_values,
                zeros(N),  # forward_storage,
                zeros(N),  # partials_storage,
                zeros(N),  # reverse_storage,
                grad_sparsity,
                hess_I,
                hess_J,
                rinfo,
                seed_matrix,
                linearity[1],
                dependent_subexpressions,
            )
        else
            return new(
                nodes,
                adj,
                const_values,
                zeros(N),  # forward_storage,
                zeros(N),  # partials_storage,
                zeros(N),  # reverse_storage,
                grad_sparsity,
                Int[],
                Int[],
                Coloring.RecoveryInfo(),
                Array{Float64}(undef, 0, 0),
                NONLINEAR,
                dependent_subexpressions,
            )
        end
    end
end

"""
    NLPEvaluator(
        data::Nonlinear.NonlinearData,
        ordered_variables::Vector{MOI.VariableIndex},
    )

Return an `MOI.AbstractNLPEvaluator`. Before using, you must initialize the
evaluator using `MOI.initialize`.
"""
mutable struct NLPEvaluator <: MOI.AbstractNLPEvaluator
    data::Nonlinear.NonlinearData
    ordered_variables::Vector{MOI.VariableIndex}

    objective::Union{Nothing,_FunctionStorage}
    constraints::Vector{_FunctionStorage}
    subexpressions::Vector{_SubexpressionStorage}
    subexpression_order::Vector{Int}
    # Storage for the subexpressions in reverse-mode automatic differentiation.
    subexpression_forward_values::Vector{Float64}
    subexpression_reverse_values::Vector{Float64}
    subexpression_linearity::Vector{Linearity}

    # A cache of the last x. This is used to guide whether we need to re-run
    # reverse-mode automatic differentiation.
    last_x::Vector{Float64}

    # Temporary storage for computing jacobians. This is also used as temporary
    # storage for the input of multivariate functions.
    jac_storage::Vector{Float64}
    # Temporary storage for the gradient of multivariate functions
    user_output_buffer::Vector{Float64}

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
    hessian_sparsity::Vector{Tuple{Int64,Int64}}
    max_chunk::Int # chunk size for which we've allocated storage

    function NLPEvaluator(
        data::Nonlinear.NonlinearData,
        ordered_variables::Vector{MOI.VariableIndex},
    )
        return new(data, ordered_variables)
    end
end

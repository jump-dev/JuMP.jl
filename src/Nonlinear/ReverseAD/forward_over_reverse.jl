#  Copyright 2017, Iain Dunning, Joey Huchette, Miles Lubin, and contributors
#  This Source Code Form is subject to the terms of the Mozilla Public
#  License, v. 2.0. If a copy of the MPL was not distributed with this
#  file, You can obtain one at https://mozilla.org/MPL/2.0/.

const TAG = :ReverseAD

"""
    _eval_hessian(
        d::NLPEvaluator,
        f::_FunctionStorage,
        H::AbstractVector{Float64},
        λ::Float64,
        offset::Int,
    )::Int

Evaluate the hessian matrix of the function `f` and store the result, scaled by
`λ`, in `H`, beginning at element `offset+1`. This function assumes that
`_reverse_mode(d, x)` has already been called.

Returns the number of non-zeros in the computed Hessian, which will be used to
update the offset for the next call.
"""
function _eval_hessian(
    d::NLPEvaluator,
    f::_FunctionStorage,
    H::AbstractVector{Float64},
    λ::Float64,
    offset::Int,
)::Int
    chunk = min(size(f.seed_matrix, 2), d.max_chunk)
    # As a performance optimization, skip dynamic dispatch if the chunk is 1.
    if chunk == 1
        return _eval_hessian_inner(d, f, H, λ, offset, Val(1))
    else
        return _eval_hessian_inner(d, f, H, λ, offset, Val(chunk))
    end
end

function _eval_hessian_inner(
    d::NLPEvaluator,
    ex::_FunctionStorage,
    H::AbstractVector{Float64},
    scale::Float64,
    nzcount::Int,
    ::Val{CHUNK},
) where {CHUNK}
    if ex.linearity == LINEAR
        @assert length(ex.hess_I) == 0
        return 0
    end
    T = ForwardDiff.Partials{CHUNK,Float64}  # This is our element type.
    Coloring.prepare_seed_matrix!(ex.seed_matrix, ex.rinfo)
    local_to_global_idx = ex.rinfo.local_indices
    input_ϵ_raw, output_ϵ_raw = d.input_ϵ, d.output_ϵ
    input_ϵ = _reinterpret_unsafe(T, input_ϵ_raw)
    output_ϵ = _reinterpret_unsafe(T, output_ϵ_raw)
    # Compute hessian-vector products
    num_products = size(ex.seed_matrix, 2) # number of hessian-vector products
    num_chunks = div(num_products, CHUNK)
    @assert size(ex.seed_matrix, 1) == length(local_to_global_idx)
    for k in 1:CHUNK:CHUNK*num_chunks
        for r in 1:length(local_to_global_idx)
            # set up directional derivatives
            @inbounds idx = local_to_global_idx[r]
            # load up ex.seed_matrix[r,k,k+1,...,k+CHUNK-1] into input_ϵ
            for s in 1:CHUNK
                input_ϵ_raw[(idx-1)*CHUNK+s] = ex.seed_matrix[r, k+s-1]
            end
            @inbounds output_ϵ[idx] = zero(T)
        end
        _hessian_slice_inner(d, ex, input_ϵ, output_ϵ, T)
        # collect directional derivatives
        for r in 1:length(local_to_global_idx)
            idx = local_to_global_idx[r]
            # load output_ϵ into ex.seed_matrix[r,k,k+1,...,k+CHUNK-1]
            for s in 1:CHUNK
                ex.seed_matrix[r, k+s-1] = output_ϵ_raw[(idx-1)*CHUNK+s]
            end
            @inbounds input_ϵ[idx] = zero(T)
        end
    end
    # leftover chunk
    remaining = num_products - CHUNK * num_chunks
    if remaining > 0
        k = CHUNK * num_chunks + 1
        for r in 1:length(local_to_global_idx)
            # set up directional derivatives
            @inbounds idx = local_to_global_idx[r]
            # load up ex.seed_matrix[r,k,k+1,...,k+remaining-1] into input_ϵ
            for s in 1:remaining
                # leave junk in the unused components
                input_ϵ_raw[(idx-1)*CHUNK+s] = ex.seed_matrix[r, k+s-1]
            end
            @inbounds output_ϵ[idx] = zero(T)
        end
        _hessian_slice_inner(d, ex, input_ϵ, output_ϵ, T)
        # collect directional derivatives
        for r in 1:length(local_to_global_idx)
            idx = local_to_global_idx[r]
            # load output_ϵ into ex.seed_matrix[r,k,k+1,...,k+remaining-1]
            for s in 1:remaining
                ex.seed_matrix[r, k+s-1] = output_ϵ_raw[(idx-1)*CHUNK+s]
            end
            @inbounds input_ϵ[idx] = zero(T)
        end
    end
    output_slice = _UnsafeVectorView(nzcount, length(ex.hess_I), pointer(H))
    Coloring.recover_from_matmat!(
        output_slice,
        ex.seed_matrix,
        ex.rinfo,
        d.output_ϵ,
    )
    for i in 1:length(output_slice)
        output_slice[i] *= scale
    end
    return length(ex.hess_I)
end

function _hessian_slice_inner(d, ex, input_ϵ, output_ϵ, ::Type{T}) where {T}
    subexpr_forward_values_ϵ =
        _reinterpret_unsafe(T, d.subexpression_forward_values_ϵ)
    for i in ex.dependent_subexpressions
        subexpr = d.subexpressions[i]
        subexpr_forward_values_ϵ[i] = _forward_eval_ϵ(
            d,
            subexpr,
            _reinterpret_unsafe(T, subexpr.forward_storage_ϵ),
            _reinterpret_unsafe(T, subexpr.partials_storage_ϵ),
            input_ϵ,
            subexpr_forward_values_ϵ,
            d.data.operators,
        )
    end
    _forward_eval_ϵ(
        d,
        ex,
        _reinterpret_unsafe(T, d.forward_storage_ϵ),
        _reinterpret_unsafe(T, d.partials_storage_ϵ),
        input_ϵ,
        subexpr_forward_values_ϵ,
        d.data.operators,
    )
    # do a reverse pass
    subexpr_reverse_values_ϵ =
        _reinterpret_unsafe(T, d.subexpression_reverse_values_ϵ)
    for i in ex.dependent_subexpressions
        subexpr_reverse_values_ϵ[i] = zero(T)
        d.subexpression_reverse_values[i] = 0.0
    end
    _reverse_eval_ϵ(
        output_ϵ,
        ex,
        _reinterpret_unsafe(T, d.reverse_storage_ϵ),
        _reinterpret_unsafe(T, d.partials_storage_ϵ),
        d.subexpression_reverse_values,
        subexpr_reverse_values_ϵ,
        1.0,
        zero(T),
    )
    for i in length(ex.dependent_subexpressions):-1:1
        j = ex.dependent_subexpressions[i]
        subexpr = d.subexpressions[j]
        _reverse_eval_ϵ(
            output_ϵ,
            subexpr,
            _reinterpret_unsafe(T, subexpr.reverse_storage_ϵ),
            _reinterpret_unsafe(T, subexpr.partials_storage_ϵ),
            d.subexpression_reverse_values,
            subexpr_reverse_values_ϵ,
            d.subexpression_reverse_values[j],
            subexpr_reverse_values_ϵ[j],
        )
    end
    return
end

"""
    _forward_eval_ϵ(
        d,
        ex::Union{_FunctionStorage,_SubexpressionStorage},
        storage_ϵ::AbstractVector{ForwardDiff.Partials{N,T}},
        partials_storage_ϵ::AbstractVector{ForwardDiff.Partials{N,T}},
        x_values_ϵ,
        subexpression_values_ϵ,
        user_operators::Nonlinear.OperatorRegistry,
    ) where {N,T}

Evaluate the directional derivatives of the expression tree in `ex`.

This is equivalent to evaluating the expression tree using DualNumbers or
ForwardDiff, but instead we keep the `ex.forward_storage` for the epsilon
components separate so that we don't need to recompute the real components.

This assumes that `_reverse_model(d, x)` has already been called.
"""
function _forward_eval_ϵ(
    d,
    ex::Union{_FunctionStorage,_SubexpressionStorage},
    storage_ϵ::AbstractVector{ForwardDiff.Partials{N,T}},
    partials_storage_ϵ::AbstractVector{ForwardDiff.Partials{N,T}},
    x_values_ϵ,
    subexpression_values_ϵ,
    user_operators::Nonlinear.OperatorRegistry,
) where {N,T}
    @assert length(storage_ϵ) >= length(ex.nodes)
    @assert length(partials_storage_ϵ) >= length(ex.nodes)
    zero_ϵ = zero(ForwardDiff.Partials{N,T})
    # ex.nodes is already in order such that parents always appear before children
    # so a backwards pass through ex.nodes is a forward pass through the tree
    children_arr = SparseArrays.rowvals(ex.adj)
    for k in length(ex.nodes):-1:1
        node = ex.nodes[k]
        partials_storage_ϵ[k] = zero_ϵ
        if node.type == Nonlinear.NODE_VARIABLE
            storage_ϵ[k] = x_values_ϵ[node.index]
        elseif node.type == Nonlinear.NODE_VALUE
            storage_ϵ[k] = zero_ϵ
        elseif node.type == Nonlinear.NODE_SUBEXPRESSION
            storage_ϵ[k] = subexpression_values_ϵ[node.index]
        elseif node.type == Nonlinear.NODE_PARAMETER
            storage_ϵ[k] = zero_ϵ
        else
            @assert node.type != Nonlinear.NODE_MOI_VARIABLE
            storage_ϵ[k] = zero_ϵ
            children_indices = SparseArrays.nzrange(ex.adj, k)
            for c_idx in children_indices
                @inbounds ix = children_arr[c_idx]
                @inbounds storage_val = storage_ϵ[ix]
                # TODO: This "if" statement can take 8% of the hessian
                # evaluation time! Find a more efficient way.
                if !isfinite(ex.partials_storage[ix]) && storage_val == zero_ϵ
                    continue
                end
                storage_ϵ[k] += storage_val * ex.partials_storage[ix]
            end
            if node.type == Nonlinear.NODE_CALL_MULTIVARIATE
                n_children = length(children_indices)
                op = user_operators.multivariate_operators[node.index]
                if op == :* && n_children == 2
                    # A performance optimization: two-argument multiplications
                    # are quite common, so we specialize on them.
                    i = children_arr[children_indices[1]]
                    j = children_arr[children_indices[2]]
                    partials_storage_ϵ[i] = storage_ϵ[j]
                    partials_storage_ϵ[j] = storage_ϵ[i]
                    continue
                end
                f_input = _UnsafeVectorView(d.jac_storage, n_children)
                for (i, c) in enumerate(children_indices)
                    f_input[i] = ex.forward_storage[children_arr[c]]
                end
                H = _UnsafeHessianView(d.user_output_buffer, n_children)
                has_hessian = Nonlinear.eval_multivariate_hessian(
                    user_operators,
                    op,
                    H,
                    f_input,
                )
                if !has_hessian
                    continue
                end
                for col in 1:n_children
                    dual = zero(ForwardDiff.Partials{N,T})
                    for row in 1:n_children
                        # Make sure we get the upper-triangular component.
                        h = row > col ? H[col, row] : H[row, col]
                        # Performance optimization: hessians can be quite sparse
                        if !iszero(h)
                            i = children_arr[children_indices[row]]
                            dual += h * storage_ϵ[i]
                        end
                    end
                    i = children_arr[children_indices[col]]
                    partials_storage_ϵ[i] = dual
                end
            elseif node.type == Nonlinear.NODE_CALL_UNIVARIATE
                @inbounds child_idx = children_arr[ex.adj.colptr[k]]
                f′′ = Nonlinear.eval_univariate_hessian(
                    user_operators,
                    user_operators.univariate_operators[node.index],
                    ex.forward_storage[child_idx],
                )
                partials_storage_ϵ[child_idx] = f′′ * storage_ϵ[child_idx]
            end
        end
    end
    return storage_ϵ[1]
end

# Compute directional derivatives of the reverse pass, goes with _forward_eval_ϵ
# to compute hessian-vector products.
function _reverse_eval_ϵ(
    output_ϵ::AbstractVector{ForwardDiff.Partials{N,T}},
    ex::Union{_FunctionStorage,_SubexpressionStorage},
    reverse_storage_ϵ,
    partials_storage_ϵ,
    subexpression_output,
    subexpression_output_ϵ,
    scale::T,
    scale_ϵ::ForwardDiff.Partials{N,T},
) where {N,T}
    @assert length(reverse_storage_ϵ) >= length(ex.nodes)
    @assert length(partials_storage_ϵ) >= length(ex.nodes)
    if ex.nodes[1].type == Nonlinear.NODE_VARIABLE
        @inbounds output_ϵ[ex.nodes[1].index] += scale_ϵ
        return
    elseif ex.nodes[1].type == Nonlinear.NODE_SUBEXPRESSION
        @inbounds subexpression_output[ex.nodes[1].index] +=
            scale * ex.reverse_storage[1]
        @inbounds subexpression_output_ϵ[ex.nodes[1].index] += scale_ϵ
        return
    end
    reverse_storage_ϵ[1] = scale_ϵ
    for k in 2:length(ex.nodes)
        @inbounds node = ex.nodes[k]
        if node.type == Nonlinear.NODE_VALUE ||
           node.type == Nonlinear.NODE_LOGIC ||
           node.type == Nonlinear.NODE_COMPARISON ||
           node.type == Nonlinear.NODE_PARAMETER
            continue
        end
        parent_value = scale * ex.reverse_storage[node.parent]
        if !isfinite(ex.partials_storage[k]) && iszero(parent_value)
            reverse_storage_ϵ[k] = zero(ForwardDiff.Partials{N,T})
        else
            reverse_storage_ϵ[k] = ForwardDiff._mul_partials(
                partials_storage_ϵ[k],
                reverse_storage_ϵ[node.parent],
                parent_value,
                ex.partials_storage[k],
            )
        end
        if node.type == Nonlinear.NODE_VARIABLE
            @inbounds output_ϵ[node.index] += reverse_storage_ϵ[k]
        elseif node.type == Nonlinear.NODE_SUBEXPRESSION
            @inbounds subexpression_output[node.index] +=
                scale * ex.reverse_storage[k]
            @inbounds subexpression_output_ϵ[node.index] += reverse_storage_ϵ[k]
        end
    end
    return
end

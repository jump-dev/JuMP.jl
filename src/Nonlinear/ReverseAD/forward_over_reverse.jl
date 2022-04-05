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
    # TODO(odow): consider reverting to a view.
    output_slice = _VectorView(nzcount, length(ex.hess_I), pointer(H))
    Coloring.recover_from_matmat!(
        output_slice,
        ex.seed_matrix,
        ex.rinfo,
        d.output_ϵ,
    )
    _rmul!(output_slice, scale)
    return length(ex.hess_I)
end

function _hessian_slice_inner(d, ex, input_ϵ, output_ϵ, ::Type{T}) where {T}
    subexpr_forward_values_ϵ =
        _reinterpret_unsafe(T, d.subexpression_forward_values_ϵ)
    for i in ex.dependent_subexpressions
        subexpr = d.subexpressions[i]
        subexpr_forward_values_ϵ[i] = _forward_eval_ϵ(
            subexpr,
            _reinterpret_unsafe(T, subexpr.forward_storage_ϵ),
            _reinterpret_unsafe(T, subexpr.partials_storage_ϵ),
            input_ϵ,
            subexpr_forward_values_ϵ,
            d.data.operators,
        )
    end
    _forward_eval_ϵ(
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
            ϵtmp = zero_ϵ
            @inbounds children_idx = SparseArrays.nzrange(ex.adj, k)
            for c_idx in children_idx
                @inbounds ix = children_arr[c_idx]
                @inbounds partial = ex.partials_storage[ix]
                @inbounds storage_val = storage_ϵ[ix]
                # TODO: This "if" statement can take 8% of the hessian
                # evaluation time! Find a more efficient way.
                if !isfinite(partial) && storage_val == zero_ϵ
                    continue
                end
                ϵtmp += storage_val * ex.partials_storage[ix]
            end
            storage_ϵ[k] = ϵtmp
            if node.type == Nonlinear.NODE_CALL_MULTIVARIATE
                # TODO(odow): consider how to refactor this into Nonlinear.
                op = node.index
                n_children = length(children_idx)
                if op == 3 # :*
                    # Lazy approach for now.
                    anyzero = false
                    tmp_prod = one(ForwardDiff.Dual{TAG,T,N})
                    for c_idx in children_idx
                        ix = children_arr[c_idx]
                        sval = ex.forward_storage[ix]
                        gnum = ForwardDiff.Dual{TAG}(sval, storage_ϵ[ix])
                        tmp_prod *= gnum
                        anyzero = ifelse(sval * sval == zero(T), true, anyzero)
                    end
                    # By a quirk of floating-point numbers, we can have
                    # anyzero == true && ForwardDiff.value(tmp_prod) != zero(T)
                    if anyzero || n_children <= 2
                        for c_idx in children_idx
                            prod_others = one(ForwardDiff.Dual{TAG,T,N})
                            for c_idx2 in children_idx
                                (c_idx == c_idx2) && continue
                                ix = children_arr[c_idx2]
                                gnum = ForwardDiff.Dual{TAG}(
                                    ex.forward_storage[ix],
                                    storage_ϵ[ix],
                                )
                                prod_others *= gnum
                            end
                            partials_storage_ϵ[children_arr[c_idx]] =
                                ForwardDiff.partials(prod_others)
                        end
                    else
                        for c_idx in children_idx
                            ix = children_arr[c_idx]
                            prod_others =
                                tmp_prod / ForwardDiff.Dual{TAG}(
                                    ex.forward_storage[ix],
                                    storage_ϵ[ix],
                                )
                            partials_storage_ϵ[ix] =
                                ForwardDiff.partials(prod_others)
                        end
                    end
                elseif op == 4 # :^
                    @assert n_children == 2
                    idx1 = first(children_idx)
                    idx2 = last(children_idx)
                    @inbounds ix1 = children_arr[idx1]
                    @inbounds ix2 = children_arr[idx2]
                    @inbounds base = ex.forward_storage[ix1]
                    @inbounds base_ϵ = storage_ϵ[ix1]
                    @inbounds exponent = ex.forward_storage[ix2]
                    @inbounds exponent_ϵ = storage_ϵ[ix2]
                    base_gnum = ForwardDiff.Dual{TAG}(base, base_ϵ)
                    exponent_gnum = ForwardDiff.Dual{TAG}(exponent, exponent_ϵ)
                    if exponent == 2
                        partials_storage_ϵ[ix1] = 2 * base_ϵ
                    elseif exponent == 1
                        partials_storage_ϵ[ix1] = zero_ϵ
                    else
                        partials_storage_ϵ[ix1] = ForwardDiff.partials(
                            exponent_gnum * pow(base_gnum, exponent_gnum - 1),
                        )
                    end
                    result_gnum = ForwardDiff.Dual{TAG}(
                        ex.forward_storage[k],
                        storage_ϵ[k],
                    )
                    # TODO(odow): fix me to use NaNMath.jl instead
                    log_base_gnum = base_gnum < 0 ? NaN : log(base_gnum)
                    partials_storage_ϵ[ix2] =
                        ForwardDiff.partials(result_gnum * log_base_gnum)
                elseif op == 5 # :/
                    @assert n_children == 2
                    idx1 = first(children_idx)
                    idx2 = last(children_idx)
                    @inbounds ix1 = children_arr[idx1]
                    @inbounds ix2 = children_arr[idx2]
                    @inbounds numerator = ex.forward_storage[ix1]
                    @inbounds numerator_ϵ = storage_ϵ[ix1]
                    @inbounds denominator = ex.forward_storage[ix2]
                    @inbounds denominator_ϵ = storage_ϵ[ix2]
                    recip_denominator =
                        1 / ForwardDiff.Dual{TAG}(denominator, denominator_ϵ)
                    partials_storage_ϵ[ix1] =
                        ForwardDiff.partials(recip_denominator)
                    partials_storage_ϵ[ix2] = ForwardDiff.partials(
                        -ForwardDiff.Dual{TAG}(numerator, numerator_ϵ) *
                        recip_denominator *
                        recip_denominator,
                    )
                elseif op > 6
                    error(
                        "User-defined operators not supported for hessian " *
                        "computations",
                    )
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

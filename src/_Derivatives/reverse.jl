
# reverse-mode evaluation of an expression tree

# assumes partials_storage is already updated
# dense gradient output, assumes initialized to zero
# if subexpressions are present, must run reverse_eval on subexpression tapes afterwards
function reverse_eval(
    reverse_storage::AbstractVector{T},
    partials_storage::AbstractVector{T},
    nd::Vector{NodeData},
    adj,
) where {T}
    @assert length(reverse_storage) >= length(nd)
    @assert length(partials_storage) >= length(nd)

    # nd is already in order such that parents always appear before children
    # so a forward pass through nd is a backwards pass through the tree

    # reverse_storage[k] is the partial derivative of the output with respect to
    # the value of node k
    reverse_storage[1] = one(T)

    for k in 2:length(nd)
        @inbounds nod = nd[k]
        if nod.nodetype == VALUE ||
           nod.nodetype == LOGIC ||
           nod.nodetype == COMPARISON ||
           nod.nodetype == PARAMETER
            continue
        end
        @inbounds rev_parent = reverse_storage[nod.parent]
        @inbounds partial = partials_storage[k]
        @inbounds reverse_storage[k] = ifelse(
            rev_parent == 0.0 && !isfinite(partial),
            rev_parent,
            rev_parent * partial,
        )
        #@inbounds reverse_storage[k] = reverse_storage[nod.parent]*partials_storage[k]
    end
    #@show storage

    return nothing
end

export reverse_eval

# assume we've already run the reverse pass, now just extract the answer
# given the scaling value
function reverse_extract(
    output::AbstractVector{T},
    reverse_storage::AbstractVector{T},
    nd::Vector{NodeData},
    adj,
    subexpression_output,
    scale_value::T,
) where {T}
    @assert length(reverse_storage) >= length(nd)

    for k in 1:length(nd)
        @inbounds nod = nd[k]
        if nod.nodetype == VARIABLE
            @inbounds output[nod.index] += scale_value * reverse_storage[k]
        elseif nod.nodetype == SUBEXPRESSION
            @inbounds subexpression_output[nod.index] +=
                scale_value * reverse_storage[k]
        end
    end
    #@show storage

    return nothing
end

export reverse_extract

# Compute directional derivatives of the reverse pass, goes with forward_eval_ϵ
# to compute hessian-vector products.
function reverse_eval_ϵ(
    output_ϵ::AbstractVector{ForwardDiff.Partials{N,T}},
    reverse_storage::AbstractVector{T},
    reverse_storage_ϵ,
    partials_storage::AbstractVector{T},
    partials_storage_ϵ,
    nd::Vector{NodeData},
    adj,
    subexpression_output,
    subexpression_output_ϵ,
    scale_value::T,
    scale_value_ϵ::ForwardDiff.Partials{N,T},
) where {N,T}
    @assert length(reverse_storage_ϵ) >= length(nd)
    @assert length(partials_storage_ϵ) >= length(nd)

    if nd[1].nodetype == VARIABLE
        @inbounds output_ϵ[nd[1].index] += scale_value_ϵ
        return
    elseif nd[1].nodetype == SUBEXPRESSION
        @inbounds subexpression_output[nd[1].index] +=
            scale_value * reverse_storage[1]
        @inbounds subexpression_output_ϵ[nd[1].index] += scale_value_ϵ
        return
    end

    reverse_storage_ϵ[1] = scale_value_ϵ

    for k in 2:length(nd)
        @inbounds nod = nd[k]
        if nod.nodetype == VALUE ||
           nod.nodetype == LOGIC ||
           nod.nodetype == COMPARISON ||
           nod.nodetype == PARAMETER
            continue
        end
        # compute the value of reverse_storage[k]
        @inbounds parentval = scale_value * reverse_storage[nod.parent]
        @inbounds parentval_ϵ = reverse_storage_ϵ[nod.parent]
        @inbounds partial = partials_storage[k]
        @inbounds partial_ϵ = partials_storage_ϵ[k]

        #reverse_storage_ϵ[k] = parentval*partial_ϵ + partial*parentval_ϵ
        if !isfinite(partial) && parentval == 0.0
            reverse_storage_ϵ[k] = zero(ForwardDiff.Partials{N,T})
        else
            reverse_storage_ϵ[k] = ForwardDiff._mul_partials(
                partial_ϵ,
                parentval_ϵ,
                parentval,
                partial,
            )
        end

        if nod.nodetype == VARIABLE
            @inbounds output_ϵ[nod.index] += reverse_storage_ϵ[k]
        elseif nod.nodetype == SUBEXPRESSION
            @inbounds subexpression_output[nod.index] +=
                scale_value * reverse_storage[k]
            @inbounds subexpression_output_ϵ[nod.index] += reverse_storage_ϵ[k]
        end
    end
    #@show storage

    return nothing
end

export reverse_eval_ϵ



# reverse-mode evaluation of an expression tree

# assumes partials_storage is already updated
# dense gradient output, assumes initialized to zero
# if subexpressions are present, must run reverse_eval on subexpression tapes afterwards
function reverse_eval{T}(output::Vector{T},reverse_storage::Vector{T},partials_storage::Vector{T},nd::Vector{NodeData},adj,subexpression_output,scale_value::T)

    @assert length(reverse_storage) >= length(nd)
    @assert length(partials_storage) >= length(nd)

    # nd is already in order such that parents always appear before children
    # so a forward pass through nd is a backwards pass through the tree

    children_arr = rowvals(adj)

    if nd[1].nodetype == VARIABLE
        output[nd[1].index] += scale_value
        return # trivial case
    elseif nd[1].nodetype == SUBEXPRESSION
        subexpression_output[nd[1].index] += scale_value
        return
    end

    # reverse_storage[k] is the partial derivative of the output with respect to
    # the value of node k
    reverse_storage[1] = scale_value

    for k in 2:length(nd)
        @inbounds nod = nd[k]
        if nod.nodetype == VALUE || nod.nodetype == LOGIC || nod.nodetype == COMPARISON || nod.nodetype == PARAMETER
            continue
        end
        # compute the value of reverse_storage[k]
        @inbounds parentval = reverse_storage[nod.parent]
        reverse_storage[k] = parentval*partials_storage[k]

        if nod.nodetype == VARIABLE
            @inbounds output[nod.index] += reverse_storage[k]
        elseif nod.nodetype == SUBEXPRESSION
            @inbounds subexpression_output[nod.index] += reverse_storage[k]
        end
    end
    #@show storage

    nothing

end

export reverse_eval


# Hessian-matrix products
# forward_input_vector should already be initialized with the input x values
function hessmat_eval!{T}(R::Matrix{T},rev_storage::Vector{Dual{T}},forward_storage::Vector{Dual{T}},partials_storage,nd::Vector{NodeData},adj,const_values,x_values::Vector{T},reverse_output_vector::Vector{Dual{T}}, forward_input_vector::Vector{Dual{T}},local_to_global_idx::Vector{Int})

    num_products = size(R,2) # number of hessian-vector products
    @assert size(R,1) == length(local_to_global_idx)
    numVar = length(x_values)

    for k in 1:num_products

        for r in 1:length(local_to_global_idx)
            # set up directional derivatives
            @inbounds idx = local_to_global_idx[r]
            @inbounds forward_input_vector[idx] = Dual(x_values[idx],R[r,k])
            @inbounds reverse_output_vector[idx] = zero(Dual{T})
        end

        # do a forward pass
        forward_eval(forward_storage,partials_storage,nd,adj,const_values,[],forward_input_vector,[])
        # do a reverse pass
        reverse_eval(reverse_output_vector,rev_storage,partials_storage,nd,adj,[],Dual(1.0))

        # collect directional derivatives
        for r in 1:length(local_to_global_idx)
            idx = local_to_global_idx[r]
            R[r,k] = epsilon(reverse_output_vector[idx])
        end

    end


end

export hessmat_eval!



# reverse-mode evaluation of an expression tree

# assumes forward_storage is already updated
# dense gradient output, assumes initialized to zero
function reverse_eval{T}(output::Vector{T},rev_storage::Vector{T},forward_storage::Vector{T},nd::Vector{NodeData},adj,const_values)

    @assert length(rev_storage) >= length(nd)
    @assert length(forward_storage) >= length(nd)

    # nd is already in order such that parents always appear before children
    # so a forward pass through nd is a backwards pass through the tree

    children_arr = rowvals(adj)

    if nd[1].nodetype == VARIABLE
        output[nd[1].index] += 1
        return # trivial case
    end

    # reverse_storage[k] is the partial derivative of the output with respect to
    # the value of node k
    rev_storage[1] = 1

    for k in 2:length(nd)
        @inbounds nod = nd[k]
        (nod.nodetype == VALUE) && continue
        # compute the value of reverse_storage[k]
        parentidx = nod.parent
        @inbounds par = nd[parentidx]
        op = par.index
        if par.nodetype == CALL
            if op == 1 # :+
                @inbounds rev_storage[k] = rev_storage[parentidx]
            elseif op == 2 # :-
                @inbounds siblings_idx = nzrange(adj,parentidx)
                if nod.whichchild == 1
                    @inbounds rev_storage[k] = rev_storage[parentidx]
                else
                    @inbounds rev_storage[k] = -rev_storage[parentidx]
                end
            elseif op == 3 # :*
                # dummy version for now
                parent_val = forward_storage[parentidx]
                if parent_val == 0.0
                    @inbounds siblings_idx = nzrange(adj,parentidx)
                    # product of all other siblings
                    prod_others = one(T)
                    for r in 1:length(siblings_idx)
                        r == nod.whichchild && continue
                        sib_idx = first(siblings_idx) + r - 1
                        @inbounds prod_others *= forward_storage[children_arr[sib_idx]]
                        prod_others == 0.0 && break
                    end
                    @inbounds rev_storage[k] = rev_storage[parentidx]*prod_others
                else
                    @inbounds rev_storage[k] = rev_storage[parentidx]*(parent_val/forward_storage[k])
                end
            elseif op == 4 # :^
                @inbounds siblings_idx = nzrange(adj,parentidx)
                if nod.whichchild == 1 # base
                    @inbounds exponentidx = children_arr[last(siblings_idx)]
                    @inbounds exponent = forward_storage[exponentidx]
                    if exponent == 2
                        @inbounds rev_storage[k] = rev_storage[parentidx]*2*forward_storage[k]
                    else
                        rev_storage[k] = rev_storage[parentidx]*exponent*forward_storage[k]^(exponent-1)
                    end
                else
                    baseidx = children_arr[first(siblings_idx)]
                    base = forward_storage[baseidx]
                    rev_storage[k] = rev_storage[parentidx]*forward_storage[parentidx]*log(base)
                end
            elseif op == 5 # :/
                @inbounds siblings_idx = nzrange(adj,parentidx)
                if nod.whichchild == 1 # numerator
                    @inbounds denomidx = children_arr[last(siblings_idx)]
                    @inbounds denom = forward_storage[denomidx]
                    @inbounds rev_storage[k] = rev_storage[parentidx]/denom
                else # denominator
                    @inbounds numeratoridx = children_arr[first(siblings_idx)]
                    @inbounds numerator = forward_storage[numeratoridx]
                    @inbounds rev_storage[k] = -rev_storage[parentidx]*numerator*forward_storage[k]^(-2)
                end
            else
                error()
            end
        else
            @assert par.nodetype == CALLUNIVAR
            @inbounds this_value = forward_storage[k]
            @inbounds rev_storage[k] = rev_storage[parentidx]*univariate_deriv(op,this_value)
        end

        if nod.nodetype == VARIABLE
            @inbounds output[nod.index] += rev_storage[k]
        end
    end
    #@show storage

    nothing

end

export reverse_eval

switchblock = Expr(:block)
for i = 1:length(univariate_operators)
    deriv_expr = univariate_operator_deriv[i]
	ex = :(return $deriv_expr::T)
    push!(switchblock.args,i,ex)
end
switchexpr = Expr(:macrocall, Expr(:.,:Lazy,quot(symbol("@switch"))), :operator_id,switchblock)

@eval @inline function univariate_deriv{T}(operator_id,x::T)
    $switchexpr
end


# Hessian-matrix products

function hessmat_eval!{T}(R::Matrix{T},rev_storage::Vector{Dual{T}},forward_storage::Vector{Dual{T}},nd::Vector{NodeData},adj,const_values,x_values::Vector{T},reverse_output_vector::Vector{Dual{T}}, forward_input_vector::Vector{Dual{T}},local_to_global_idx::Vector{Int})

    num_products = size(R,2) # number of hessian-vector products
    @assert size(R,1) == length(local_to_global_idx)
    numVar = length(x_values)

    for i in 1:numVar
        forward_input_vector[i] = Dual(x_values[i],0.0)
    end

    for k in 1:num_products

        for r in 1:length(local_to_global_idx)
            # set up directional derivatives
            idx = local_to_global_idx[r]
            forward_input_vector[idx] = Dual(x_values[idx],R[r,k])
            reverse_output_vector[idx] = zero(Dual{T})
        end

        # do a forward pass
        forward_eval(forward_storage,nd,adj,const_values,forward_input_vector)
        # do a reverse pass
        reverse_eval(reverse_output_vector,rev_storage,forward_storage,nd,adj,const_values)

        # collect directional derivatives
        for r in 1:length(local_to_global_idx)
            idx = local_to_global_idx[r]
            R[r,k] = epsilon(reverse_output_vector[idx])
        end

    end


end

export hessmat_eval!



# reverse-mode evaluation of an expression tree

# assumes forward_storage is already updated
# dense gradient output, assumes initialized to zero
# if subexpressions are present, must run reverse_eval on subexpression tapes afterwards
function reverse_eval{T}(output::Vector{T},rev_storage::Vector{T},forward_storage::Vector{T},partials_storage::Vector{T},nd::Vector{NodeData},adj,subexpression_output,scale_value::T)

    @assert length(rev_storage) >= length(nd)
    @assert length(forward_storage) >= length(nd)

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
    rev_storage[1] = scale_value

    for k in 2:length(nd)
        @inbounds nod = nd[k]
        if nod.nodetype == VALUE || nod.nodetype == LOGIC || nod.nodetype == COMPARISON || nod.nodetype == PARAMETER
            continue
        end
        # compute the value of reverse_storage[k]
        parentidx = nod.parent
        @inbounds par = nd[parentidx]
        @inbounds parentval = rev_storage[parentidx]
        op = par.index
        if par.nodetype == CALL
            if op == 1 # :+
                @inbounds rev_storage[k] = parentval
            elseif op == 2 # :-
                @inbounds siblings_idx = nzrange(adj,parentidx)
                if nod.whichchild == 1
                    @inbounds rev_storage[k] = parentval
                else
                    @inbounds rev_storage[k] = -parentval
                end
            elseif op == 3 # :*
                # dummy version for now
                @inbounds siblings_idx = nzrange(adj,parentidx)
                n_siblings = length(siblings_idx)
                if n_siblings == 2
                    otheridx = ifelse(nod.whichchild == 1, last(siblings_idx),first(siblings_idx))
                    @inbounds prod_others = forward_storage[children_arr[otheridx]]
                    @inbounds rev_storage[k] = parentval*prod_others
                else
                    @inbounds parent_val = forward_storage[parentidx]
                    if parent_val == 0.0
                        # product of all other siblings
                        prod_others = one(T)
                        for r in 1:n_siblings
                            r == nod.whichchild && continue
                            sib_idx = first(siblings_idx) + r - 1
                            @inbounds prod_others *= forward_storage[children_arr[sib_idx]]
                            prod_others == 0.0 && break
                        end
                        @inbounds rev_storage[k] = parentval*prod_others
                    else
                        @inbounds rev_storage[k] = parentval*(parent_val/forward_storage[k])
                    end
                end
            elseif op == 4 # :^
                @inbounds siblings_idx = nzrange(adj,parentidx)
                if nod.whichchild == 1 # base
                    @inbounds exponentidx = children_arr[last(siblings_idx)]
                    @inbounds exponent = forward_storage[exponentidx]
                    if exponent == 2
                        @inbounds rev_storage[k] = parentval*2*forward_storage[k]
                    else
                        rev_storage[k] = parentval*exponent*pow(forward_storage[k],exponent-1)
                    end
                else
                    baseidx = children_arr[first(siblings_idx)]
                    base = forward_storage[baseidx]
                    rev_storage[k] = parentval*forward_storage[parentidx]*log(base)
                end
            elseif op == 5 # :/
                @inbounds siblings_idx = nzrange(adj,parentidx)
                if nod.whichchild == 1 # numerator
                    @inbounds denomidx = children_arr[last(siblings_idx)]
                    @inbounds denom = forward_storage[denomidx]
                    @inbounds rev_storage[k] = parentval/denom
                else # denominator
                    @inbounds numeratoridx = children_arr[first(siblings_idx)]
                    @inbounds numerator = forward_storage[numeratoridx]
                    @inbounds rev_storage[k] = -parentval*numerator*pow(forward_storage[k],-2)
                end
            elseif op == 6 # ifelse
                @inbounds siblings_idx = nzrange(adj,parentidx)
                conditionidx = first(siblings_idx)
                @inbounds condition = (forward_storage[children_arr[conditionidx]]==1)
                if (condition && nod.whichchild == 2) || (!condition && nod.whichchild == 3)
                    rev_storage[k] = parentval
                else
                    rev_storage[k] = zero(T)
                end
            else
                error()
            end
        elseif par.nodetype == LOGIC || par.nodetype == COMPARISON
            # nonlinear, but these don't go into the derivatives
            rev_storage[k] = zero(T)
        else
            @assert par.nodetype == CALLUNIVAR
            @inbounds this_fprime = partials_storage[k]
            @inbounds rev_storage[k] = parentval*this_fprime
        end

        if nod.nodetype == VARIABLE
            @inbounds output[nod.index] += rev_storage[k]
        elseif nod.nodetype == SUBEXPRESSION
            @inbounds subexpression_output[nod.index] += rev_storage[k]
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
        reverse_eval(reverse_output_vector,rev_storage,forward_storage,partials_storage,nd,adj,[],Dual(1.0))

        # collect directional derivatives
        for r in 1:length(local_to_global_idx)
            idx = local_to_global_idx[r]
            R[r,k] = epsilon(reverse_output_vector[idx])
        end

    end


end

export hessmat_eval!

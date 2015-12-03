

# reverse-mode evaluation of an expression tree

# assumes forward_storage is already updated
# dense gradient output, assumes initialized to zero
function reverse_eval{T}(output::Vector{T},rev_storage::Vector{T},forward_storage::Vector{T},nd::Vector{NodeData},adj,const_values,x_values::Vector{T})

    @assert length(forward_storage) == length(rev_storage) == length(nd)

    # nd is already in order such that parents always appear before children
    # so a forward pass through nd is a backwards pass through the tree

    children_arr = rowvals(adj)

    # reverse_storage[k] is the partial derivative of the output with respect to
    # the value of node k
    rev_storage[1] = 1
    @assert nd[1].nodetype != VARIABLE

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
                @inbounds rev_storage[k] = rev_storage[parentidx]*(forward_storage[parentidx]/forward_storage[k])
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
	ex = :(return $deriv_expr)
    push!(switchblock.args,i,ex)
end
switchexpr = Expr(:macrocall, Expr(:.,:Lazy,quot(symbol("@switch"))), :operator_id,switchblock)

@eval @inline function univariate_deriv(operator_id,x)
    $switchexpr
end

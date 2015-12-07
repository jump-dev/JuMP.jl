

# forward-mode evaluation of an expression tree
# tree is represented as Vector{NodeData} and adjacency matrix from adjmat()

function forward_eval{T}(storage::Vector{T},nd::Vector{NodeData},adj,const_values,x_values::Vector{T})

    @assert length(storage) >= length(nd)

    # nd is already in order such that parents always appear before children
    # so a backwards pass through nd is a forward pass through the tree

    children_arr = rowvals(adj)

    for k in length(nd):-1:1
        # compute the value of node k
        @inbounds nod = nd[k]
        if nod.nodetype == VARIABLE
            @inbounds storage[k] = x_values[nod.index]
        elseif nod.nodetype == VALUE
            @inbounds storage[k] = const_values[nod.index]
        elseif nod.nodetype == CALL
            op = nod.index
            @inbounds children_idx = nzrange(adj,k)
            #@show children_idx
            n_children = length(children_idx)
            if op == 1 # :+
                tmp_sum = zero(T)
                # sum over children
                for c_idx in children_idx
                    @inbounds tmp_sum += storage[children_arr[c_idx]]
                end
                storage[k] = tmp_sum
            elseif op == 2 # :-
                child1 = first(children_idx)
                @inbounds tmp_sub = storage[children_arr[child1]]
                @assert n_children == 2
                @inbounds tmp_sub -= storage[children_arr[child1+1]]
                storage[k] = tmp_sub
            elseif op == 3 # :*
                tmp_prod = one(T)
                for c_idx in children_idx
                    @inbounds tmp_prod *= storage[children_arr[c_idx]]
                end
                storage[k] = tmp_prod
            elseif op == 4 # :^
                @assert n_children == 2
                idx1 = first(children_idx)
                idx2 = last(children_idx)
                @inbounds base = storage[children_arr[idx1]]
                @inbounds exponent = storage[children_arr[idx2]]
                if exponent == 2
                    storage[k] = base*base
                else
                    storage[k] = base^exponent
                end
            elseif op == 5 # :/
                @assert n_children == 2
                idx1 = first(children_idx)
                idx2 = last(children_idx)
                @inbounds numerator = storage[children_arr[idx1]]
                @inbounds denominator = storage[children_arr[idx2]]
                storage[k] = numerator/denominator
            else
                error()
            end
        elseif nod.nodetype == CALLUNIVAR # univariate function
            op = nod.index
            @inbounds child_idx = children_arr[adj.colptr[k]]
            #@assert child_idx == children_arr[first(nzrange(adj,k))]
            child_val = storage[child_idx]
            @inbounds storage[k] = eval_univariate(op, child_val)
        end

    end
    #@show storage

    return storage[1]

end

export forward_eval


switchblock = Expr(:block)
for i = 1:length(univariate_operators)
    op = univariate_operators[i]
	ex = :(return $op(v))
    push!(switchblock.args,i,ex)
end
switchexpr = Expr(:macrocall, Expr(:.,:Lazy,quot(symbol("@switch"))), :operator_id,switchblock)

@eval @inline function eval_univariate(operator_id,v)
    $switchexpr
end

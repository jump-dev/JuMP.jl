

# forward-mode evaluation of an expression tree
# tree is represented as Vector{NodeData} and adjacency matrix from adjmat()

function forward_eval{T}(storage::Vector{T},nd::Vector{NodeData},adj,const_values,x_values::Vector{T})

    @assert length(storage) == length(nd)

    # nd is already in order such that parents always appear before children
    # so a backwards pass through nd is a forward pass through the tree

    children_arr = rowvals(adj)

    for k in length(nd):-1:1
        # compute the value of node k
        nod = nd[k]
        if nod.nodetype == VARIABLE
            storage[k] = x_values[nod.index]
        elseif nod.nodetype == VALUE
            storage[k] = const_values[nod.index]
        elseif nod.nodetype == CALL
            op = nod.index
            children_idx = nzrange(adj,k)
            #@show children_idx
            n_children = length(children_idx)
            # can generate this code later
            if op == 1 # :+
                tmp_sum = zero(T)
                # sum over children
                for c_idx in children_idx
                    @inbounds tmp_sum += storage[children_arr[c_idx]]
                end
                storage[k] = tmp_sum
            elseif op == 2 # :-
                child1 = first(children_idx)
                tmp_sub = storage[children_arr[child1]]
                if n_children == 1
                    storage[k] = -tmp_sub
                else
                    @assert n_children == 2
                    tmp_sub -= storage[children_arr[child1+1]] 
                    storage[k] = tmp_sub
                end
            elseif op == 3 # :*
                tmp_prod = one(T)
                for c_idx in children_idx
                    @inbounds tmp_prod *= storage[children_arr[c_idx]]
                end
                storage[k] = tmp_prod
            elseif op == 4 # :sin
                @assert n_children == 1
                storage[k] = sin(storage[children_arr[first(children_idx)]])
            elseif op == 5 # :cos
                @assert n_children == 1
                storage[k] = cos(storage[children_arr[first(children_idx)]])
            elseif op == 6 # :^
                @assert n_children == 2
                idx1 = children_idx[1] 
                idx2 = children_idx[2] 
                @inbounds base = storage[children_arr[idx1]]
                @inbounds exponent = storage[children_arr[idx2]]
                if exponent == 2
                    storage[k] = base*base
                else
                    storage[k] = base^exponent
                end
            else
                error()
            end
        end
    end
    #@show storage

    return storage[1]

end

export forward_eval

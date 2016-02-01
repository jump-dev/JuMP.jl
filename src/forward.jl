

# forward-mode evaluation of an expression tree
# tree is represented as Vector{NodeData} and adjacency matrix from adjmat()
# assumes values of expressions have already been computed
# partials_storage[k] is the partial derivative of nd[k].parent with respect to the value
# of node k. It's efficient to compute this at the same time as the value of the parent.
# Since we use it in reverse mode and in dual forward mode.
# Note that partials_storage makes a subtle assumption that we have a tree instead of
# a general DAG. If we have a DAG, then need to associate storage with each edge of the DAG.
function forward_eval{T}(storage::Vector{T},partials_storage::Vector{T},nd::Vector{NodeData},adj,const_values,parameter_values,x_values::Vector{T},subexpression_values)

    @assert length(storage) >= length(nd)
    @assert length(partials_storage) >= length(nd)

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
        elseif nod.nodetype == SUBEXPRESSION
            @inbounds storage[k] = subexpression_values[nod.index]
        elseif nod.nodetype == PARAMETER
            @inbounds storage[k] = parameter_values[nod.index]
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
                tmp_prod = ProductAccumulator(T)
                for c_idx in children_idx
                    @inbounds tmp_prod = add_term(tmp_prod,storage[children_arr[c_idx]])
                end
                product_value = tmp_prod.allbutone*tmp_prod.smallest
                for c_idx in children_idx
                    @inbounds ix = children_arr[c_idx]
                    @inbounds child_value = storage[ix]
                    if child_value == tmp_prod.smallest
                        @inbounds partials_storage[ix] = tmp_prod.allbutone
                    else
                        @inbounds partials_storage[ix] = product_value/child_value
                    end
                end
                @inbounds storage[k] = product_value
            elseif op == 4 # :^
                @assert n_children == 2
                idx1 = first(children_idx)
                idx2 = last(children_idx)
                @inbounds ix1 = children_arr[idx1]
                @inbounds ix2 = children_arr[idx2]
                @inbounds base = storage[ix1]
                @inbounds exponent = storage[ix2]
                if exponent == 2
                    @inbounds storage[k] = base*base
                    @inbounds partials_storage[ix1] = 2*base
                else
                    storage[k] = pow(base,exponent)
                    partials_storage[ix1] = exponent*pow(base,exponent-1)
                end
            elseif op == 5 # :/
                @assert n_children == 2
                idx1 = first(children_idx)
                idx2 = last(children_idx)
                @inbounds ix1 = children_arr[idx1]
                @inbounds ix2 = children_arr[idx2]
                @inbounds numerator = storage[ix1]
                @inbounds denominator = storage[ix2]
                # only store the partial wrt numerator because it's cheap
                # partial wrt denominator may not be needed if constant
                recip_denominator = 1/denominator
                @inbounds partials_storage[ix1] = recip_denominator
                storage[k] = numerator*recip_denominator
            elseif op == 6 # ifelse
                @assert n_children == 3
                idx1 = first(children_idx)
                @inbounds condition = storage[children_arr[idx1]]
                @inbounds lhs = storage[children_arr[idx1+1]]
                @inbounds rhs = storage[children_arr[idx1+2]]
                @inbounds partials_storage[children_arr[idx1+1]] = condition == 1
                @inbounds partials_storage[children_arr[idx1+2]] = !(condition == 1)
                storage[k] = ifelse(condition == 1, lhs, rhs)
            else
                error("Unsupported operation $(operators[op])")
            end
        elseif nod.nodetype == CALLUNIVAR # univariate function
            op = nod.index
            @inbounds child_idx = children_arr[adj.colptr[k]]
            #@assert child_idx == children_arr[first(nzrange(adj,k))]
            child_val = storage[child_idx]
            fval, fprimeval = eval_univariate(op, child_val)
            @inbounds partials_storage[child_idx] = fprimeval
            @inbounds storage[k] = fval
        elseif nod.nodetype == COMPARISON
            op = nod.index

            @inbounds children_idx = nzrange(adj,k)
            n_children = length(children_idx)
            result = true
            for r in 1:n_children-1
                cval_lhs = storage[children_arr[children_idx[r]]]
                cval_rhs = storage[children_arr[children_idx[r+1]]]
                if op == 1
                    result &= cval_lhs <= cval_rhs
                elseif op == 2
                    result &= cval_lhs == cval_rhs
                elseif op == 3
                    result &= cval_lhs >= cval_rhs
                elseif op == 4
                    result &= cval_lhs < cval_rhs
                elseif op == 5
                    result &= cval_lhs > cval_rhs
                end
            end
            storage[k] = result
        elseif nod.nodetype == LOGIC
            op = nod.index
            @inbounds children_idx = nzrange(adj,k)
            # boolean values are stored as floats
            cval_lhs = (storage[children_arr[first(children_idx)]] == 1)
            cval_rhs = (storage[children_arr[last(children_idx)]] == 1)
            if op == 1
                storage[k] = cval_lhs && cval_rhs
            elseif op == 2
                storage[k] = cval_lhs || cval_rhs
            end
        end

    end
    #@show storage

    return storage[1]

end

export forward_eval

# Evaluate directional derivatives of an expression tree.
# This is equivalent to evaluating the expression tree using DualNumbers or ForwardDiff,
# but instead we keep the storage for the epsilon components separate so that we don't
# need to recompute the real components.
# We assume that forward_eval has already been called.
function forward_eval_dual{T,N}(storage::Vector{T},storage_ϵ::Vector{ForwardDiff.PartialsTup{T,N}},nd::Vector{NodeData},adj,x_values_ϵ::Vector{ForwardDiff.PartialsTup{T,N}},subexpression_values_ϵ)

    @assert length(storage) >= length(nd)
    @assert length(storage_ϵ) >= length(nd)

    zero_ϵ = zero_partials(eltype(storage_ϵ))

    children_arr = rowvals(adj)

    for k in length(nd):-1:1
        # compute the value of node k
        @inbounds nod = nd[k]
        if nod.nodetype == VARIABLE
            @inbounds storage_ϵ[k] = x_values_ϵ[nod.index]
        elseif nod.nodetype == VALUE
            @inbounds storage_ϵ[k] = zero_ϵ
        elseif nod.nodetype == SUBEXPRESSION
            @inbounds storage_ϵ[k] = subexpression_values_ϵ[nod.index]
        elseif nod.nodetype == PARAMETER
            @inbounds storage_ϵ[k] = zero_ϵ
        elseif nod.nodetype == CALL
            op = nod.index
            @inbounds children_idx = nzrange(adj,k)
            #@show children_idx
            n_children = length(children_idx)
            if op == 1 # :+
                tmp_sum = zero_ϵ
                # sum over children
                for c_idx in children_idx
                    @inbounds tmp_sum += storage_ϵ[children_arr[c_idx]]
                end
                storage_ϵ[k] = tmp_sum
            elseif op == 2 # :-
                child1 = first(children_idx)
                @inbounds tmp_sub = storage_ϵ[children_arr[child1]]
                @assert n_children == 2
                @inbounds tmp_sub -= storage_ϵ[children_arr[child1+1]]
                storage_ϵ[k] = tmp_sub
            elseif op == 3 # :*
                # Not much to be saved by skipping multiplication of real parts
                tmp_prod = one(GradientNumber{N,T,NTuple{N,T}})
                for c_idx in children_idx
                    @inbounds ix = children_arr[c_idx]
                    @inbounds gnum = GradientNumber(storage[ix], storage_ϵ[ix])
                    @inbounds tmp_prod *= gnum
                end
                storage_ϵ[k] = grad(tmp_prod)
            elseif op == 4 # :^
                @assert n_children == 2
                idx1 = first(children_idx)
                idx2 = last(children_idx)
                @inbounds base = storage[children_arr[idx1]]
                @inbounds base_ϵ = storage_ϵ[children_arr[idx1]]
                @inbounds exponent = storage[children_arr[idx2]]
                @inbounds exponent_ϵ = storage[children_arr[idx2]]
                # lazy approach for now
                storage_ϵ[k] = grad(pow(GradientNumber(base,base_ϵ),GradientNumber(exponent,exponent_ϵ)))
            elseif op == 5 # :/
                @assert n_children == 2
                idx1 = first(children_idx)
                idx2 = last(children_idx)
                @inbounds numerator = storage[children_arr[idx1]]
                @inbounds numerator_ϵ = storage_ϵ[children_arr[idx1]]
                @inbounds denominator = storage[children_arr[idx2]]
                @inbounds denominator_ϵ = storage_ϵ[children_arr[idx2]]
                storage_ϵ[k] = grad(GradientNumber(numerator,numerator_ϵ)/GradientNumbers(denominator,denominator_ϵ))
            elseif op == 6 # ifelse
                @assert n_children == 3
                idx1 = first(children_idx)
                @inbounds condition = storage[children_arr[idx1]]
                @inbounds lhs = storage_ϵ[children_arr[idx1+1]]
                @inbounds rhs = storage_ϵ[children_arr[idx1+2]]
                storage_ϵ[k] = ifelse(condition == 1, lhs, rhs)
            else
                error("Unsupported operation $(operators[op])")
            end
        elseif nod.nodetype == CALLUNIVAR # univariate function
            op = nod.index
            @inbounds child_idx = children_arr[adj.colptr[k]]
            #@assert child_idx == children_arr[first(nzrange(adj,k))]
            child_val = storage[child_idx]
            @inbounds storage[k] = eval_univariate(op, child_val)
        elseif nod.nodetype == COMPARISON
            # nothing to do
        elseif nod.nodetype == LOGIC
            # nothing to do
        end

    end
    #@show storage

    return storage[1]

end


switchblock = Expr(:block)
for i = 1:length(univariate_operators)
    op = univariate_operators[i]
    deriv_expr = univariate_operator_deriv[i]
	ex = :(return $op(x), $deriv_expr::T)
    push!(switchblock.args,i,ex)
end
switchexpr = Expr(:macrocall, Expr(:.,:Lazy,quot(symbol("@switch"))), :operator_id,switchblock)

@eval @inline function eval_univariate{T}(operator_id,x::T)
    $switchexpr
end

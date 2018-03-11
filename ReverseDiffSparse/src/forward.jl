

# forward-mode evaluation of an expression tree
# tree is represented as Vector{NodeData} and adjacency matrix from adjmat()
# assumes values of expressions have already been computed
# partials_storage[k] is the partial derivative of nd[k].parent with respect to the value
# of node k. It's efficient to compute this at the same time as the value of the parent.
# Since we use it in reverse mode and in dual forward mode.
# Note that partials_storage makes a subtle assumption that we have a tree instead of
# a general DAG. If we have a DAG, then need to associate storage with each edge of the DAG.
# user_input_buffer and user_output_buffer are used as temporary storage
# when handling user-defined functions
function forward_eval{T}(storage::Vector{T},partials_storage::Vector{T},nd::Vector{NodeData},adj,const_values,parameter_values,x_values::Vector{T},subexpression_values,user_input_buffer=[],user_output_buffer=[];user_operators::UserOperatorRegistry=UserOperatorRegistry())

    @assert length(storage) >= length(nd)
    @assert length(partials_storage) >= length(nd)

    # nd is already in order such that parents always appear before children
    # so a backwards pass through nd is a forward pass through the tree

    children_arr = rowvals(adj)

    for k in length(nd):-1:1
        # compute the value of node k
        @inbounds nod = nd[k]
        partials_storage[k] = zero(T)
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
                    @inbounds ix = children_arr[c_idx]
                    @inbounds partials_storage[ix] = one(T)
                    @inbounds tmp_sum += storage[ix]
                end
                storage[k] = tmp_sum
            elseif op == 2 # :-
                child1 = first(children_idx)
                @assert n_children == 2
                @inbounds ix1 = children_arr[child1]
                @inbounds ix2 = children_arr[child1+1]
                @inbounds tmp_sub = storage[ix1]
                @inbounds tmp_sub -= storage[ix2]
                @inbounds partials_storage[ix1] = one(T)
                @inbounds partials_storage[ix2] = -one(T)
                storage[k] = tmp_sub
            elseif op == 3 # :*
                tmp_prod = one(T)
                for c_idx in children_idx
                    @inbounds tmp_prod *= storage[children_arr[c_idx]]
                end
                if tmp_prod == zero(T) # inefficient
                    for c_idx in children_idx
                        prod_others = one(T)
                        for c_idx2 in children_idx
                            (c_idx == c_idx2) && continue
                            prod_others *= storage[children_arr[c_idx2]]
                        end
                        partials_storage[children_arr[c_idx]] = prod_others
                    end
                else
                    for c_idx in children_idx
                        ix = children_arr[c_idx]
                        partials_storage[ix] = tmp_prod/storage[ix]
                    end
                end
                @inbounds storage[k] = tmp_prod
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
                partials_storage[ix2] = storage[k]*log(base)
            elseif op == 5 # :/
                @assert n_children == 2
                idx1 = first(children_idx)
                idx2 = last(children_idx)
                @inbounds ix1 = children_arr[idx1]
                @inbounds ix2 = children_arr[idx2]
                @inbounds numerator = storage[ix1]
                @inbounds denominator = storage[ix2]
                recip_denominator = 1/denominator
                @inbounds partials_storage[ix1] = recip_denominator
                partials_storage[ix2] = -numerator*recip_denominator*recip_denominator
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
            elseif op >= USER_OPERATOR_ID_START
                evaluator = user_operators.multivariate_operator_evaluator[op - USER_OPERATOR_ID_START+1]
                f_input = view(user_input_buffer, 1:n_children)
                grad_output = view(user_output_buffer, 1:n_children)
                r = 1
                for c_idx in children_idx
                    ix = children_arr[c_idx]
                    f_input[r] = storage[ix]
                    grad_output[r] = 0.0
                    r += 1
                end
                fval = MathProgBase.eval_f(evaluator, f_input)::T
                MathProgBase.eval_grad_f(evaluator, grad_output, f_input)
                storage[k] = fval
                r = 1
                for c_idx in children_idx
                    ix = children_arr[c_idx]
                    partials_storage[ix] = grad_output[r]
                    r += 1
                end
            else
                error("Unsupported operation $(operators[op])")
            end
        elseif nod.nodetype == CALLUNIVAR # univariate function
            op = nod.index
            @inbounds child_idx = children_arr[adj.colptr[k]]
            #@assert child_idx == children_arr[first(nzrange(adj,k))]
            child_val = storage[child_idx]
            if op >= USER_UNIVAR_OPERATOR_ID_START
                userop = op - USER_UNIVAR_OPERATOR_ID_START + 1
                f = user_operators.univariate_operator_f[userop]
                fprime = user_operators.univariate_operator_fprime[userop]
                fval = f(child_val)::T
                fprimeval = fprime(child_val)::T
            else
                fval, fprimeval = eval_univariate(op, child_val)
            end
            @inbounds partials_storage[child_idx] = fprimeval
            @inbounds storage[k] = fval
        elseif nod.nodetype == COMPARISON
            op = nod.index

            @inbounds children_idx = nzrange(adj,k)
            n_children = length(children_idx)
            result = true
            for r in 1:n_children-1
                partials_storage[children_arr[children_idx[r]]] = zero(T)
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
            partials_storage[children_arr[children_idx[n_children]]] = zero(T)
            storage[k] = result
        elseif nod.nodetype == LOGIC
            op = nod.index
            @inbounds children_idx = nzrange(adj,k)
            # boolean values are stored as floats
            partials_storage[children_arr[first(children_idx)]] = zero(T)
            partials_storage[children_arr[last(children_idx)]] = zero(T)
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
# Computes partials_storage_ϵ as well
# We assume that forward_eval has already been called.
function forward_eval_ϵ{N,T}(storage::Vector{T},storage_ϵ::DenseVector{ForwardDiff.Partials{N,T}},partials_storage::Vector{T},partials_storage_ϵ::DenseVector{ForwardDiff.Partials{N,T}},nd::Vector{NodeData},adj,x_values_ϵ,subexpression_values_ϵ;user_operators::UserOperatorRegistry=UserOperatorRegistry())

    @assert length(storage_ϵ) >= length(nd)
    @assert length(partials_storage_ϵ) >= length(nd)

    zero_ϵ = zero(ForwardDiff.Partials{N,T})

    # nd is already in order such that parents always appear before children
    # so a backwards pass through nd is a forward pass through the tree

    children_arr = rowvals(adj)

    for k in length(nd):-1:1
        # compute the value of node k
        @inbounds nod = nd[k]
        partials_storage_ϵ[k] = zero_ϵ
        if nod.nodetype == VARIABLE
            @inbounds storage_ϵ[k] = x_values_ϵ[nod.index]
        elseif nod.nodetype == VALUE
            @inbounds storage_ϵ[k] = zero_ϵ
        elseif nod.nodetype == SUBEXPRESSION
            @inbounds storage_ϵ[k] = subexpression_values_ϵ[nod.index]
        elseif nod.nodetype == PARAMETER
            @inbounds storage_ϵ[k] = zero_ϵ
        else
            ϵtmp = zero_ϵ
            @inbounds children_idx = nzrange(adj,k)
            for c_idx in children_idx
                ix = children_arr[c_idx]
                @inbounds partial = partials_storage[ix]
                if !isfinite(partials_storage[ix]) && storage_ϵ[ix] == zero_ϵ
                    continue
                end
                ϵtmp += storage_ϵ[ix]*partials_storage[ix]
            end
            storage_ϵ[k] = ϵtmp

            if nod.nodetype == CALL
                op = nod.index
                n_children = length(children_idx)
                if op == 3 # :*
                    # Lazy approach for now
                    anyzero = false
                    tmp_prod = one(ForwardDiff.Dual{TAG,T,N})
                    for c_idx in children_idx
                        ix = children_arr[c_idx]
                        sval = storage[ix]
                        gnum = ForwardDiff.Dual{TAG}(sval,storage_ϵ[ix])
                        tmp_prod *= gnum
                        anyzero = ifelse(sval*sval == zero(T), true, anyzero)
                    end
                    # By a quirk of floating-point numbers, we can have
                    # anyzero == true && ForwardDiff.value(tmp_prod) != zero(T)
                    if anyzero # inefficient
                        for c_idx in children_idx
                            prod_others = one(ForwardDiff.Dual{TAG,T,N})
                            for c_idx2 in children_idx
                                (c_idx == c_idx2) && continue
                                ix = children_arr[c_idx2]
                                gnum = ForwardDiff.Dual{TAG}(storage[ix],storage_ϵ[ix])
                                prod_others *= gnum
                            end
                            partials_storage_ϵ[children_arr[c_idx]] = ForwardDiff.partials(prod_others)
                        end
                    else
                        for c_idx in children_idx
                            ix = children_arr[c_idx]
                            prod_others = tmp_prod/ForwardDiff.Dual{TAG}(storage[ix],storage_ϵ[ix])
                            partials_storage_ϵ[ix] = ForwardDiff.partials(prod_others)
                        end
                    end
                elseif op == 4 # :^
                    @assert n_children == 2
                    idx1 = first(children_idx)
                    idx2 = last(children_idx)
                    @inbounds ix1 = children_arr[idx1]
                    @inbounds ix2 = children_arr[idx2]
                    @inbounds base = storage[ix1]
                    @inbounds base_ϵ = storage_ϵ[ix1]
                    @inbounds exponent = storage[ix2]
                    @inbounds exponent_ϵ = storage_ϵ[ix2]
                    base_gnum = ForwardDiff.Dual{TAG}(base,base_ϵ)
                    exponent_gnum = ForwardDiff.Dual{TAG}(exponent,exponent_ϵ)
                    if exponent == 2
                        partials_storage_ϵ[ix1] = 2*base_ϵ
                    else
                        partials_storage_ϵ[ix1] = ForwardDiff.partials(exponent_gnum*pow(base_gnum,exponent_gnum-1))
                    end
                    result_gnum = ForwardDiff.Dual{TAG}(storage[k],storage_ϵ[k])
                    partials_storage_ϵ[ix2] = ForwardDiff.partials(result_gnum*log(base_gnum))
                elseif op == 5 # :/
                    @assert n_children == 2
                    idx1 = first(children_idx)
                    idx2 = last(children_idx)
                    @inbounds ix1 = children_arr[idx1]
                    @inbounds ix2 = children_arr[idx2]
                    @inbounds numerator = storage[ix1]
                    @inbounds numerator_ϵ = storage_ϵ[ix1]
                    @inbounds denominator = storage[ix2]
                    @inbounds denominator_ϵ = storage_ϵ[ix2]
                    recip_denominator = 1/ForwardDiff.Dual{TAG}(denominator,denominator_ϵ)
                    partials_storage_ϵ[ix1] = ForwardDiff.partials(recip_denominator)
                    partials_storage_ϵ[ix2] = ForwardDiff.partials(-ForwardDiff.Dual{TAG}(numerator,numerator_ϵ)*recip_denominator*recip_denominator)
                elseif op >= USER_OPERATOR_ID_START
                    error("User-defined operators not supported for hessian computations")
                end
            elseif nod.nodetype == CALLUNIVAR # univariate function
                op = nod.index
                @inbounds child_idx = children_arr[adj.colptr[k]]
                child_val = storage[child_idx]
                if op >= USER_UNIVAR_OPERATOR_ID_START
                    userop = op - USER_UNIVAR_OPERATOR_ID_START + 1
                    fprimeprime = user_operators.univariate_operator_fprimeprime[userop](child_val)::T
                else
                    fprimeprime = eval_univariate_2nd_deriv(op, child_val,storage[k])
                end
                partials_storage_ϵ[child_idx] = fprimeprime*storage_ϵ[child_idx]
            end
        end

    end

    return storage_ϵ[1]

end

export forward_eval_ϵ


exprs = Expr[]
for i = 1:length(univariate_operators)
    op = univariate_operators[i]
    deriv_expr = univariate_operator_deriv[i]
    ex = :(return $op(x), $deriv_expr::T)
    push!(exprs, ex)
end

function binaryswitch(ids, exprs)
    if length(exprs) <= 3
        out = Expr(:if, Expr(:call, :(==), :operator_id, ids[1]), exprs[1])
        if length(exprs) > 1
            push!(out.args, binaryswitch(ids[2:end], exprs[2:end]))
        end
        return out
    else
        mid = length(exprs) >>> 1
        return Expr(:if, Expr(:call, :(<=), :operator_id, ids[mid]),
            binaryswitch(ids[1:mid], exprs[1:mid]),
            binaryswitch(ids[mid+1:end], exprs[mid+1:end]))
    end
end
switchexpr = binaryswitch(1:length(exprs), exprs)

@eval @inline function eval_univariate{T}(operator_id,x::T)
    $switchexpr
    error("No match for operator_id")
end

# TODO: optimize sin/cos/exp
ids = Int[]
exprs = Expr[]
for i = 1:length(univariate_operators)
    op = univariate_operators[i]
    if op == :asec || op == :acsc || op == :asecd || op == :acscd || op == :acsch || op == :trigamma
        # there's an abs in the derivative that Calculus can't symbolically differentiate
        continue
    end
    if i in 1:3 # :+, :-, :abs
        deriv_expr = :(zero(T))
    elseif op == :sin || op == :cos
        deriv_expr = :(-fval)
    elseif op == :exp
        deriv_expr = :(fval)
    else
        deriv_expr = Calculus.differentiate(univariate_operator_deriv[i],:x)
    end
    ex = :(return $deriv_expr::T)
    push!(ids, i)
    push!(exprs, ex)
end
switchexpr = binaryswitch(ids, exprs)

@eval @inline function eval_univariate_2nd_deriv{T}(operator_id,x::T,fval::T)
    $switchexpr
    error("No match for operator_id")
end

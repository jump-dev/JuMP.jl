

# returns the list of subexpressions which a given tape depends on directly
function list_subexpressions(nd::Vector{NodeData})

    # TODO: code overlap with compute_gradient_sparsity
    indices = Set{Int}()
    for k in 1:length(nd)
        nod = nd[k]
        if nod.nodetype == SUBEXPRESSION
            push!(indices, nod.index)
        end
    end

    return sort(collect(indices))
end

# order the subexpressions which main_expressions depend on
# such that we can run forward mode in this order
function order_subexpressions(main_expressions::Vector{Vector{NodeData}},subexpressions::Vector{Vector{NodeData}})
    nsub = length(subexpressions)
    computed = falses(nsub)
    dependencies = Dict{Int,Vector{Int}}()
    to_visit = collect(nsub+1:nsub+length(main_expressions))
    depended_on_by = [Set{Int}() for i in 1:nsub]

    while !isempty(to_visit)
        idx = pop!(to_visit)
        if idx > nsub
            li = list_subexpressions(main_expressions[idx-nsub])
        else
            computed[idx] && continue
            li = list_subexpressions(subexpressions[idx])
            computed[idx] = true
        end
        dependencies[idx] = li
        for k in li
            if idx > nsub
                push!(depended_on_by[k], idx-nsub)
            else
                union!(depended_on_by[k], depended_on_by[idx])
            end
            push!(to_visit,k)
        end
    end

    # now order dependencies
    I = Int[]
    J = Int[]
    for (idx,li) in dependencies
        for k in li
            push!(I,idx)
            push!(J,k)
        end
    end
    N = nsub+length(main_expressions)
    sp = sparse(I,J,ones(length(I)),N,N)
    cmap = Array{Int}(N)
    order = reverse(Coloring.reverse_topological_sort_by_dfs(sp.rowval,sp.colptr,N,cmap)[1])
    # remove the subexpressions which never appear anywhere
    # and the indices of the main expressions
    order_filtered = collect(filter(idx -> (idx <= nsub && computed[idx]), order))
    # also generate an individual order for each main expression
    individual_order = [Int[] for i in 1:length(main_expressions)]
    for o in order_filtered
        for i in depended_on_by[o]
            push!(individual_order[i],o)
        end
    end
    
    return order_filtered, individual_order
end

export list_subexpressions, order_subexpressions


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

# Order the subexpressions which main_expressions depend on such that we can
# run forward mode in this order.
function order_subexpressions(
    main_expressions::Vector{Vector{NodeData}},
    subexpressions::Vector{Vector{NodeData}},
)
    num_sub = length(subexpressions)
    computed = falses(num_sub)
    dependencies = Dict{Int,Vector{Int}}()
    to_visit = collect(num_sub+1:num_sub+length(main_expressions))
    # For each subexpression k, the indices of the main expressions that depend
    # on k, possibly transitively.
    depended_on_by = [Set{Int}() for i in 1:num_sub]
    while !isempty(to_visit)
        idx = pop!(to_visit)
        if idx > num_sub
            subexpr = list_subexpressions(main_expressions[idx-num_sub])
        else
            computed[idx] && continue
            subexpr = list_subexpressions(subexpressions[idx])
            computed[idx] = true
        end
        dependencies[idx] = subexpr
        for k in subexpr
            if idx > num_sub
                push!(depended_on_by[k], idx - num_sub)
            end
            push!(to_visit, k)
        end
    end

    # Now order dependencies.
    I = Int[]
    J = Int[]
    for (idx, subexpr) in dependencies
        for k in subexpr
            push!(I, idx)
            push!(J, k)
        end
    end
    N = num_sub + length(main_expressions)
    sp = sparse(I, J, ones(length(I)), N, N)
    cmap = Vector{Int}(undef, N)

    order = reverse(Coloring.reverse_topological_sort_by_dfs(
        sp.rowval,
        sp.colptr,
        N,
        cmap,
    )[1])
    # Remove the subexpressions which never appear anywhere and the indices of
    # the main expressions.
    condition(idx) = idx <= num_sub && computed[idx]
    order_filtered = collect(filter(condition, order))

    # Propagate transitive dependencies.
    for o in Iterators.reverse(order_filtered)
        @assert !isempty(depended_on_by[o])
        for k in list_subexpressions(subexpressions[o])
            union!(depended_on_by[k], depended_on_by[o])
        end
    end

    # Generate an individual order for each main expression.
    individual_order = [Int[] for i in 1:length(main_expressions)]
    for o in order_filtered
        for i in depended_on_by[o]
            push!(individual_order[i], o)
        end
    end

    return order_filtered, individual_order
end

# An implementation of depth-first-search for topologically sorting the
# subexpressions needed to compute main_expression.
function order_subexpression(
    main_expression::Vector{NodeData},
    subexpressions::Vector{Vector{NodeData}},
)
    L = Int[]
    marked_nodes = Set{Int}()
    unmarked_nodes = Set{Int}(
        node.index for node in main_expression if node.nodetype == SUBEXPRESSION
    )
    function visit(n::Int)
        if n in marked_nodes
            return
        end
        for m in subexpressions[n]
            if m.nodetype == SUBEXPRESSION
                visit(m.index)
            end
        end
        push!(marked_nodes, n)
        push!(L, n)
        return
    end
    while !isempty(unmarked_nodes)
        n = pop!(unmarked_nodes)
        visit(n)
    end
    return L
end

export list_subexpressions, order_subexpressions, order_subexpression

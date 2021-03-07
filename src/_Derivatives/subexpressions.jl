
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

"""
Return the list of dependent subexpressions. Similar to `list_subexpressions`,
but without the unnecessary `Set{Int}` -> `Vector{Int}` -> `sort`.
"""
function _subexpressions(nd::Vector{NodeData})
    return Int[n.index for n in nd if n.nodetype == SUBEXPRESSION]
end

"""
Return `G[node]`, lazily building it from `subexpressions` if it is not yet
assigned.
"""
function _children(
    G::Vector{Vector{Int}},
    subexpressions::Vector{Vector{NodeData}},
    node::Int,
)
    if !isassigned(G, node)
        G[node] = _subexpressions(subexpressions[node])
    end
    return G[node]
end

"""
    _topological_sort(starts, G::Vector{Vector{Int}})

Return a topologically sorted list of nodes in `G` that need to be evaluated
to compute the nodes contained in `starts`.

## Notes

* It is important to not use recursion here, because expressions may have
  arbitrary levels of nesting!
* This function assumes `G` is acyclic.
"""
function _topological_sort(
    starts,
    G::Vector{Vector{Int}},
    subexpressions::Vector{Vector{NodeData}},
)
    ordered = Int[]
    visited = fill(false, length(G))
    stack = Tuple{Int,Bool}[]
    for s in starts
        if visited[s]
            continue  # This node has already been visited.
        else
            push!(stack, (s, true))
            visited[s] = true
        end
        while !isempty(stack)
            # Get a new node from the top of the stack.
            node, needs_checking = pop!(stack)
            # If !needs_checking, we must be returning to this node for a second
            # time, and we have already checked all of the children. Therefore,
            # we can add it to the set of ordered nodes.
            if !needs_checking
                push!(ordered, node)
                continue
            end
            # Re-add the node to the stack, but set the `false` flag this time
            # so next time we visit it, it will go on the `ordered` list
            # instead.
            push!(stack, (node, false))
            # Now check all of the children.
            for child in _children(G, subexpressions, node)
                if visited[child]
                    continue  # This node has already been visited.
                else
                    # Add the node to the stack. Because we haven't visited it
                    # before, we need to check all of the children.
                    push!(stack, (child, true))
                    visited[child] = true
                end
            end
        end
    end
    return ordered
end

"""
    order_subexpressions(
        main_expressions::Vector{Vector{NodeData}},
        subexpressions::Vector{Vector{NodeData}};
    )

Topologically sort the subexpressions needed to evaluate `main_expressions`.

Returns two things:

 * A `Vector{Int64}` containing the ordered list of subexpression-indices that
   need to be evaluated to compute all `main_expressions`
 * A `Vector{Vector{Int64}}`, containing a list of ordered lists of
   subexpression-indices that need to be evaluated to compute
   `main_expressions[i]`.

**Warning:** This doesn't handle cyclic expressions! But this should be fine
because we can't compute them in JuMP anyway.
"""
function order_subexpressions(
    main_expressions::Vector{Vector{NodeData}},
    subexpressions::Vector{Vector{NodeData}},
)
    # The graph of node dependencies. Constructed lazily.
    G = Vector{Vector{Int}}(undef, length(subexpressions))
    # Node dependencies of the main expressions.
    starts = Set{Int}()
    individual_sorts = Vector{Int}[]
    for expression in main_expressions
        s = _subexpressions(expression)
        union!(starts, s)
        push!(individual_sorts, _topological_sort(s, G, subexpressions))
    end
    full_sort = _topological_sort(starts, G, subexpressions)
    return full_sort, individual_sorts
end

"""
    order_subexpressions(
        main_expression::Vector{NodeData},
        subexpressions::Vector{Vector{NodeData}},
    )

Return a topologically ordered list of subexpression-indices needed to evaluate
`main_expression`.
"""
function order_subexpressions(
    main_expression::Vector{NodeData},
    subexpressions::Vector{Vector{NodeData}},
)
    s = _subexpressions(main_expression)
    G = Vector{Vector{Int}}(undef, length(subexpressions))
    return _topological_sort(s, G, subexpressions)
end

export list_subexpressions, order_subexpressions

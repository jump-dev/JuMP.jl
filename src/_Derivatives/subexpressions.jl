
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
    _subexpressions(nd::Vector{NodeData})::Vector{Int}

Return the integer-indices of subexpressions which the tape `nd` depends on
directly.

This list may contain duplicates. If you want a list of unique indices, use
`list_subexpressions` instead.
"""
function _subexpressions(nd::Vector{NodeData})
    return Int[n.index for n in nd if n.nodetype == SUBEXPRESSION]
end

"""
    _topological_sort(
        starts,
        subexpressions::Vector{Vector{NodeData}},
        G::Vector{Vector{Int}} = Vector{Vector{Int}}(undef, subexpressions),
    )

Return a topologically sorted list of the integer subexpresssion indices that
need to be computed to evalute `subexpressions[s]` for all `s in starts`.
(`starts` may contain duplicates, and just needs to be iterable.)

`G[i]` is a lazily-computed list of "out" edges from node `i`, in terms of the
integer-valued subexpression index (i.e., `node.index`). This list may contain
duplicates.

If calling `_topological_sort` a single time, you may omit the `G` argument.

However, if calling `_topological_sort` multiple times on the _same_ vector of
subexpresssions, you should create `G` once (either as the uninitialized vector,
or by explicitly computing the full `G`), and pass it in.

## Notes

* It is important to not use recursion here, because expressions may have
  arbitrary levels of nesting!
* This function assumes `subexpressions` is acyclic.
"""
function _topological_sort(
    starts,
    subexpressions::Vector{Vector{NodeData}},
    G::Vector{Vector{Int}} = Vector{Vector{Int}}(undef, subexpressions),
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
            if !isassigned(G, node)
                G[node] = _subexpressions(subexpressions[node])
            end
            for child in G[node]
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

Topologically sort the subexpression needed to evaluate `main_expressions`.

Returns two things:

 * A `Vector{Int}` containing the ordered list of subexpression-indices that
   need to be evaluated to compute all `main_expressions`
 * A `Vector{Vector{Int}}`, containing a list of ordered lists of
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
        push!(individual_sorts, _topological_sort(s, subexpressions, G))
    end
    full_sort = _topological_sort(starts, subexpressions, G)
    return full_sort, individual_sorts
end

"""
    order_subexpressions(
        main_expression::Vector{NodeData},
        subexpressions::Vector{Vector{NodeData}},
    )

Return a `Vector{Int}` containing the topologically ordered list of
subexpression-indices that need to be evaluated to compute `main_expression`.
"""
function order_subexpressions(
    main_expression::Vector{NodeData},
    subexpressions::Vector{Vector{NodeData}},
)
    s = _subexpressions(main_expression)
    return _topological_sort(s, subexpressions)
end

export list_subexpressions, order_subexpressions

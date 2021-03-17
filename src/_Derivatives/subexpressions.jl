
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
    _topological_sort(
        starts::Vector{Int},
        subexpressions::Vector{Vector{NodeData}},
        subexpression_dependency_graph::Vector{Vector{Int}} =
            Vector{Vector{Int}}(undef, length(subexpressions)),
    )

Return a topologically sorted list of the integer subexpresssion indices that
need to be computed to evalute `subexpressions[s]` for all `s in starts`.

`starts` should be ordered, and not contain duplicates.

`subexpression_dependency_graph[i]` is a lazily-computed list of "out" edges
from node `i`, in terms of the integer-valued subexpression index (i.e.,
`node.index`). This list should be unique and ordered.

If calling `_topological_sort` a single time, you may omit the
`subexpression_dependency_graph` argument.

However, if calling `_topological_sort` multiple times on the _same_ vector of
subexpresssions, you should create `subexpression_dependency_graph` once (either
as the uninitialized vector, or by explicitly computing the full
`subexpression_dependency_graph`), and pass it in.

## Notes

* It is important to not use recursion here, because expressions may have
  arbitrary levels of nesting!
* This function assumes `subexpressions` is acyclic.
"""
function _topological_sort(
    starts,
    subexpressions::Vector{Vector{NodeData}},
    subexpression_dependency_graph::Vector{Vector{Int}} = Vector{Vector{Int}}(
        undef,
        length(subexpressions),
    ),
)
    ordered = Int[]
    in_order = fill(false, length(subexpressions))
    stack = Tuple{Int,Bool}[]
    for s in starts
        if in_order[s]
            continue  # s is already in `ordered`.
        end
        push!(stack, (s, true))
        while !isempty(stack)
            node, needs_checking = pop!(stack)
            if !needs_checking
                # We must be returning to this node for a second time, and we
                # have already checked all of the children. Therefore, we can
                # add it to the set of ordered nodes.
                push!(ordered, node)
                in_order[node] = true
                continue
            elseif in_order[node]
                continue  # This node has already been added to `ordered`.
            end
            # Re-add the node to the stack, but set the `false` flag this time
            # so next time we visit it, it will go on the `ordered` list
            # instead.
            push!(stack, (node, false))
            if !isassigned(subexpression_dependency_graph, node)
                subexpression_dependency_graph[node] =
                    list_subexpressions(subexpressions[node])
            end
            for child in subexpression_dependency_graph[node]
                if !in_order[child]
                    push!(stack, (child, true))
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
    subexpression_dependency_graph =
        Vector{Vector{Int}}(undef, length(subexpressions))
    # Node dependencies of the main expressions.
    starts = Set{Int}()
    individual_sorts = Vector{Int}[]
    for expression in main_expressions
        s = list_subexpressions(expression)
        union!(starts, s)
        push!(
            individual_sorts,
            _topological_sort(
                s,
                subexpressions,
                subexpression_dependency_graph,
            ),
        )
    end
    full_sort = _topological_sort(
        starts,
        subexpressions,
        subexpression_dependency_graph,
    )
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
    s = list_subexpressions(main_expression)
    return _topological_sort(s, subexpressions)
end

export list_subexpressions, order_subexpressions

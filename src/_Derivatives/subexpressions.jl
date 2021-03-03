
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

function _visit_nodes(
    subexpressions::Vector{Vector{NodeData}},
    n::Int64,
    individual::Tuple{Vector{Int64},Set{Int64}},
    full::Union{Nothing,Tuple{Vector{Int64},Set{Int64}}},
)
    # First, check if we have already visited node `n`. If so, bail from the
    # recursion.
    if n in individual[2]
        return
    end
    # Else, recurse on each SUBEXPRESSION node in subexpression n.
    for m in subexpressions[n]
        if m.nodetype == SUBEXPRESSION
            _visit_nodes(subexpressions, m.index, individual, full)
        end
    end
    # If we made it here, it means that all the children subexpressions (if they
    # exist) are added, so we're free to append subexpression `n` to the list.
    push!(individual[1], n)
    push!(individual[2], n)
    # In addition, if we are storing the full list, add `n` here too.
    if full !== nothing && !(n in full[2])
        push!(full[1], n)
        push!(full[2], n)
    end
    return
end

"""
    order_subexpressions(
        main_expression::Vector{NodeData},
        subexpressions::Vector{Vector{NodeData}},
        full::Union{Nothing,Tuple{Vector{Int64},Set{Int64}}} = nothing,
    )

Return a topologically ordered list of subexpression-indices needed to evaluate
`main_expression`.

`full` is an argument used by `order_subexpressions` to collate a list of
expressions needed to evaluate multiple `main_expressions`.

**Warning:** This doesn't handle cyclic expressions! But this should be fine
because we can't compute them in JuMP anyway.
"""
function order_subexpressions(
    main_expression::Vector{NodeData},
    subexpressions::Vector{Vector{NodeData}},
    full::Union{Nothing,Tuple{Vector{Int64},Set{Int64}}} = nothing,
)
    # `order` is a topologically sorted list of nodes. Nodes at the start of the
    # list need to be evaluated first.
    # `tagged` is a set to efficiently keep track of nodes that we have added to
    # `order`. We could just use `n in order`, but a `Set` is probably faster.
    order, tagged = Int64[], Set{Int64}()
    # The meat of DFS is this recursive call to `_visit_node` for each
    # SUBEXPRESSION node in `main_expression`.
    for node in main_expression
        if node.nodetype == SUBEXPRESSION
            _visit_nodes(subexpressions, node.index, (order, tagged), full)
        end
    end
    return order
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
    full = (Int64[], Set{Int64}())
    individual = Vector{Int64}[
        order_subexpressions(expr, subexpressions, full)
        for expr in main_expressions
    ]
    return full[1], individual
end

export list_subexpressions, order_subexpressions

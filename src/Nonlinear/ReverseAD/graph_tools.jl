#  Copyright 2017, Iain Dunning, Joey Huchette, Miles Lubin, and contributors
#  This Source Code Form is subject to the terms of the Mozilla Public
#  License, v. 2.0. If a copy of the MPL was not distributed with this
#  file, You can obtain one at https://mozilla.org/MPL/2.0/.

"""
    _replace_moi_variables(
        nodes::Vector{Nonlinear.Node},
        moi_index_to_consecutive_index::Dict{MOI.VariableIndex,Int},
    )

Return a new `Vector{Nonlinear.Node}` where all occurances of
`NODE_MOI_VARIABLE` are replaced by `NODE_VARIABLE` that is 1-indexed and
ordered.
"""
function _replace_moi_variables(
    nodes::Vector{Nonlinear.Node},
    moi_index_to_consecutive_index::Dict{MOI.VariableIndex,Int},
)
    new_nodes = Vector{Nonlinear.Node}(undef, length(nodes))
    for (i, node) in enumerate(nodes)
        if node.type == Nonlinear.NODE_MOI_VARIABLE
            new_nodes[i] = Nonlinear.Node(
                Nonlinear.NODE_VARIABLE,
                moi_index_to_consecutive_index[MOI.VariableIndex(node.index)],
                node.parent,
            )
        else
            new_nodes[i] = node
        end
    end
    return new_nodes
end

@enum(Linearity, CONSTANT, LINEAR, PIECEWISE_LINEAR, NONLINEAR)

"""
    _classify_linearity(
        nodes::Vector{Nonlinear.Node},
        adj::SparseArrays.SparseMatrixCSC,
        subexpression_linearity::Vector{Linearity},
    )

Classify the nodes in a tree as constant, linear, or nonlinear with respect to
the input.
"""
function _classify_linearity(
    nodes::Vector{Nonlinear.Node},
    adj::SparseArrays.SparseMatrixCSC,
    subexpression_linearity::Vector{Linearity},
)
    linearity = Array{Linearity}(undef, length(nodes))
    children_arr = SparseArrays.rowvals(adj)
    for k in length(nodes):-1:1
        node = nodes[k]
        if node.type == Nonlinear.NODE_VARIABLE
            linearity[k] = LINEAR
            continue
        elseif node.type == Nonlinear.NODE_VALUE
            linearity[k] = CONSTANT
            continue
        elseif node.type == Nonlinear.NODE_PARAMETER
            linearity[k] = CONSTANT
            continue
        elseif node.type == Nonlinear.NODE_SUBEXPRESSION
            linearity[k] = subexpression_linearity[node.index]
            continue
        end
        children_idx = SparseArrays.nzrange(adj, k)
        num_constant_children, any_nonlinear = 0, false
        for r in children_idx
            if linearity[children_arr[r]] == NONLINEAR
                any_nonlinear = true
                break
            elseif linearity[children_arr[r]] == CONSTANT
                num_constant_children += 1
            end
        end
        if any_nonlinear
            # If any children are nonlinear, then we're nonlinear...
            linearity[k] = NONLINEAR
            # ...except in the case of ifelse. If the operands are linear then
            # we're piecewise linear.
            op = get(
                Nonlinear.DEFAULT_MULTIVARIATE_OPERATORS,
                node.index,
                nothing,
            )
            if (
                node.type == Nonlinear.NODE_CALL_MULTIVARIATE &&
                op == :ifelse &&
                linearity[children_arr[children_idx[2]]] == LINEAR &&
                linearity[children_arr[children_idx[3]]] == LINEAR
            )
                linearity[k] = PIECEWISE_LINEAR
            end
            continue
        elseif num_constant_children == length(children_idx)
            # If all children are constant, then we're constant.
            linearity[k] = CONSTANT
            continue
        end
        # By this point, some children are constant and some are linear, so if
        # the operator is nonlinear, then we're nonlinear.
        if node.type == Nonlinear.NODE_CALL_UNIVARIATE
            op =
                get(Nonlinear.DEFAULT_UNIVARIATE_OPERATORS, node.index, nothing)
            if op == :+ || op == :-
                linearity[k] = LINEAR
            else
                linearity[k] = NONLINEAR
            end
        elseif node.type == Nonlinear.NODE_CALL_MULTIVARIATE
            op = get(
                Nonlinear.DEFAULT_MULTIVARIATE_OPERATORS,
                node.index,
                nothing,
            )
            if op == :+
                linearity[k] = LINEAR
            elseif op == :-
                linearity[k] = LINEAR
            elseif op == :*
                # Multiplication is linear if there is one non-constant term.
                one_op = num_constant_children == length(children_idx) - 1
                linearity[k] = one_op ? LINEAR : NONLINEAR
            elseif op == :^
                linearity[k] = NONLINEAR
            elseif op == :/
                if linearity[children_arr[children_idx[2]]] == CONSTANT
                    # If the denominator is constant, we're linear.
                    linearity[k] = LINEAR
                else
                    linearity[k] = NONLINEAR
                end
            elseif op == :ifelse
                linearity[k] = NONLINEAR
            else  # User-defined functions
                linearity[k] = NONLINEAR
            end
        elseif node.type == Nonlinear.NODE_LOGIC
            linearity[k] = NONLINEAR
        else
            @assert node.type == Nonlinear.NODE_COMPARISON
            linearity[k] = NONLINEAR
        end
    end
    return linearity
end

"""
    _compute_gradient_sparsity!(
        indices::Coloring.IndexedSet,
        nodes::Vector{Nonlinear.Node},
    )

Compute the sparsity pattern of the gradient of an expression (i.e., a list of
which variable indices are present).
"""
function _compute_gradient_sparsity!(
    indices::Coloring.IndexedSet,
    nodes::Vector{Nonlinear.Node},
)
    for node in nodes
        if node.type == Nonlinear.NODE_VARIABLE
            push!(indices, node.index)
        elseif node.type == Nonlinear.NODE_MOI_VARIABLE
            error(
                "Internal error: Invalid to compute sparsity if Nonlinear.NODE_MOI_VARIABLE " *
                "nodes are present.",
            )
        end
    end
    return
end

"""
    _compute_hessian_sparsity(
        nodes::Vector{Nonlinear.Node},
        adj,
        input_linearity::Vector{Linearity},
        indexedset::Coloring.IndexedSet,
        subexpression_edgelist::Vector{Set{Tuple{Int,Int}}},
        subexpression_variables::Vector{Vector{Int}},
    )

Compute the sparsity pattern the Hessian of an expression.

 * `input_linearity` is the linearity with respect to the input, computed by
   `_classify_linearity`
 * `subexpression_edgelist` is the edge_list of each subexpression
 * `subexpression_variables` is the list of all variables which appear in a
   subexpression (including recursively).

Idea: consider the (non)linearity of a node *with respect to the output*. The
children of any node which is nonlinear with respect to the output should have
nonlinear interactions, hence nonzeros in the hessian. This is not true in
general, but holds for everything we consider.

A counter example is `f(x, y, z) = x + y * z`, but we don't have any functions
like that. By "nonlinear with respect to the output", we mean that the output
depends nonlinearly on the value of the node, regardless of how the node itself
depends on the input.
"""
function _compute_hessian_sparsity(
    nodes::Vector{Nonlinear.Node},
    adj,
    input_linearity::Vector{Linearity},
    indexedset::Coloring.IndexedSet,
    subexpression_edgelist::Vector{Set{Tuple{Int,Int}}},
    subexpression_variables::Vector{Vector{Int}},
)
    # So start at the root of the tree and classify the linearity wrt the output.
    # For each nonlinear node, do a mini DFS and collect the list of children.
    # Add a nonlinear interaction between all children of a nonlinear node.
    edge_list = Set{Tuple{Int,Int}}()
    nonlinear_wrt_output = fill(false, length(nodes))
    children_arr = SparseArrays.rowvals(adj)
    stack = Int[]
    stack_ignore = Bool[]
    nonlinear_group = indexedset
    if length(nodes) == 1 && nodes[1].type == Nonlinear.NODE_SUBEXPRESSION
        # Subexpression comes in linearly, so append edge_list
        for ij in subexpression_edgelist[nodes[1].index]
            push!(edge_list, ij)
        end
    end
    for k in 2:length(nodes)
        nod = nodes[k]
        if nod.type == Nonlinear.NODE_MOI_VARIABLE
            error(
                "Internal error: Invalid to compute sparsity if Nonlinear.NODE_MOI_VARIABLE " *
                "nodes are present.",
            )
        end
        if nonlinear_wrt_output[k]
            continue # already seen this node one way or another
        elseif input_linearity[k] == CONSTANT
            continue # definitely not nonlinear
        end
        @assert !nonlinear_wrt_output[nod.parent]
        # check if the parent depends nonlinearly on the value of this node
        par = nodes[nod.parent]
        if par.type == Nonlinear.NODE_CALL_UNIVARIATE
            op = get(Nonlinear.DEFAULT_UNIVARIATE_OPERATORS, par.index, nothing)
            if op === nothing || (op != :+ && op != :-)
                nonlinear_wrt_output[k] = true
            end
        elseif par.type == Nonlinear.NODE_CALL_MULTIVARIATE
            op = get(
                Nonlinear.DEFAULT_MULTIVARIATE_OPERATORS,
                par.index,
                nothing,
            )
            if op === nothing
                nonlinear_wrt_output[k] = true
            elseif op in (:+, :-, :ifelse)
                # pass
            elseif op == :*
                # check if all siblings are constant
                sibling_idx = SparseArrays.nzrange(adj, nod.parent)
                if !all(
                    i ->
                        input_linearity[children_arr[i]] == CONSTANT ||
                            children_arr[i] == k,
                    sibling_idx,
                )
                    # at least one sibling isn't constant
                    nonlinear_wrt_output[k] = true
                end
            elseif op == :/
                # check if denominator is nonconstant
                sibling_idx = SparseArrays.nzrange(adj, nod.parent)
                if input_linearity[children_arr[last(sibling_idx)]] != CONSTANT
                    nonlinear_wrt_output[k] = true
                end
            else
                nonlinear_wrt_output[k] = true
            end
        end
        if nod.type == Nonlinear.NODE_SUBEXPRESSION && !nonlinear_wrt_output[k]
            # subexpression comes in linearly, so append edge_list
            for ij in subexpression_edgelist[nod.index]
                push!(edge_list, ij)
            end
        end
        if !nonlinear_wrt_output[k]
            continue
        end
        # do a DFS from here, including all children
        @assert isempty(stack)
        @assert isempty(stack_ignore)
        sibling_idx = SparseArrays.nzrange(adj, nod.parent)
        for sidx in sibling_idx
            push!(stack, children_arr[sidx])
            push!(stack_ignore, false)
        end
        empty!(nonlinear_group)
        while length(stack) > 0
            r = pop!(stack)
            should_ignore = pop!(stack_ignore)
            nonlinear_wrt_output[r] = true
            if nodes[r].type == Nonlinear.NODE_LOGIC ||
               nodes[r].type == Nonlinear.NODE_COMPARISON
                # don't count the nonlinear interactions inside
                # logical conditions or comparisons
                should_ignore = true
            end
            children_idx = SparseArrays.nzrange(adj, r)
            for cidx in children_idx
                push!(stack, children_arr[cidx])
                push!(stack_ignore, should_ignore)
            end
            if should_ignore
                continue
            end
            if nodes[r].type == Nonlinear.NODE_VARIABLE
                push!(nonlinear_group, nodes[r].index)
            elseif nodes[r].type == Nonlinear.NODE_SUBEXPRESSION
                # append all variables in subexpression
                union!(nonlinear_group, subexpression_variables[nodes[r].index])
            end
        end
        for i_ in 1:nonlinear_group.nnz
            i = nonlinear_group.nzidx[i_]
            for j_ in 1:nonlinear_group.nnz
                j = nonlinear_group.nzidx[j_]
                if j > i
                    continue  # Only lower triangle.
                end
                push!(edge_list, (i, j))
            end
        end
    end
    return edge_list
end

"""
    _list_subexpressions(nodes::Vector{Nonlinear.Node})

Returns the list of subexpressions which a given tape depends on directly
"""
function _list_subexpressions(nodes::Vector{Nonlinear.Node})
    indices = Set{Int}(
        n.index for n in nodes if n.type == Nonlinear.NODE_SUBEXPRESSION
    )
    return sort(collect(indices))
end

"""
    _topological_sort(
        starts::Vector{Int},
        subexpressions::Vector{Vector{Nonlinear.Node}},
        subexpression_dependency_graph::Vector{Vector{Int}} =
            Vector{Vector{Int}}(undef, length(subexpressions)),
    )

Return a topologically sorted list of the integer subexpresssion indices that
need to be computed to evaluate `subexpressions[s]` for all `s in starts`.

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
    subexpressions::Vector{Vector{Nonlinear.Node}},
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
                    _list_subexpressions(subexpressions[node])
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
    _order_subexpressions(
        main_expressions::Vector{Vector{Nonlinear.Node}},
        subexpressions::Vector{Vector{Nonlinear.Node}};
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
function _order_subexpressions(
    main_expressions::Vector{Vector{Nonlinear.Node}},
    subexpressions::Vector{Vector{Nonlinear.Node}},
)
    # The graph of node dependencies. Constructed lazily.
    subexpression_dependency_graph =
        Vector{Vector{Int}}(undef, length(subexpressions))
    # Node dependencies of the main expressions.
    starts = Set{Int}()
    individual_sorts = Vector{Int}[]
    for expression in main_expressions
        s = _list_subexpressions(expression)
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
    _order_subexpressions(
        main_expression::Vector{Nonlinear.Node},
        subexpressions::Vector{Vector{Nonlinear.Node}},
    )

Return a `Vector{Int}` containing the topologically ordered list of
subexpression-indices that need to be evaluated to compute `main_expression`.
"""
function _order_subexpressions(
    main_expression::Vector{Nonlinear.Node},
    subexpressions::Vector{Vector{Nonlinear.Node}},
)
    s = _list_subexpressions(main_expression)
    return _topological_sort(s, subexpressions)
end

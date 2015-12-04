
# compute the sparsity pattern of the gradient of an expression
# (i.e., a list of which indices are present)

function compute_gradient_sparsity(nd::Vector{NodeData},adj)

    # this is easy, just build a list of which variable indices appear
    indices = Set{Int}()
    for k in 1:length(nd)
        nod = nd[k]
        if nod.nodetype == VARIABLE
            push!(indices, nod.index)
        end
    end

    return sort(collect(indices))
end

export compute_gradient_sparsity

# compute the sparsity pattern the hessian of an expression

# input_linearity is the linearity with respect to the input
# computed by classify_linearity
function compute_hessian_sparsity(nd::Vector{NodeData},adj,input_linearity::Vector{Linearity})

    # idea: consider the linearity/nonlinearity of a node *with respect to the output*
    # The children of any node which is nonlinear with respect to the output
    # should have nonlinear interactions, hence nonzeros in the hessian.
    # This is not true in general but holds for everything we consider.
    # Counterexample: f(x,y,z) = x + y*z, but we don't have any functions like that.
    # By "nonlinear with respect to the output", we mean that the output
    # depends nonlinearly on the value of the node, regardless of
    # how the node itself depends on the input.
   
    # So start at the root of the tree and classify the linearity wrt the output.
    # For each nonlinear node, do a mini DFS and collect the list of children.
    # Add a nonlinear interaction between all children of a nonlinear node.
    
    edgelist = Set{Tuple{Int,Int}}()
    nonlinear_wrt_output = fill(false,length(nd))
    
    children_arr = rowvals(adj)

    stack = Int[]
    nonlinear_group = Set{Int}() # TODO: replace with indexed set
    
    for k in 2:length(nd)
        nod = nd[k]
        nonlinear_wrt_output[k] && continue # already seen this node one way or another
        input_linearity[k] == CONSTANT && continue # definitely not nonlinear

        @assert !nonlinear_wrt_output[nod.parent]
        # check if the parent depends nonlinearly on the value of this node
        par = nd[nod.parent]
        if par.nodetype == CALLUNIVAR
            op = par.index
            if univariate_operators[op] != :+ && univariate_operators[op] != :-
                nonlinear_wrt_output[k] = true
            end
        elseif par.nodetype == CALL
            op = par.index
            if operators[op] == :+ || operators[op] == :-
                continue
            elseif operators[op] == :*
                # check if all siblings are constant
                sibling_idx = nzrange(adj,nod.parent)
                if !all(i -> input_linearity[children_arr[i]] == CONSTANT || children_arr[i] == k, sibling_idx)
                    # at least one sibling isn't constant
                    nonlinear_wrt_output[k] = true
                end
            else # TODO: division
                nonlinear_wrt_output[k] = true
            end
        end
        
        nonlinear_wrt_output[k] || continue
        
        # do a DFS from here, including all children
        @assert isempty(stack)
        sibling_idx = nzrange(adj,nod.parent)
        for sidx in sibling_idx
            push!(stack, children_arr[sidx])
        end
        empty!(nonlinear_group)
        while length(stack) > 0
            r = pop!(stack)
            nonlinear_wrt_output[r] = true
            children_idx = nzrange(adj,r)
            for cidx in children_idx
                push!(stack, children_arr[cidx])
            end
            if nd[r].nodetype == VARIABLE
                push!(nonlinear_group, nd[r].index)
            end
        end
        for i in nonlinear_group
            for j in nonlinear_group
                j <= i || continue # only lower triangle
                push!(edgelist,(i,j))
            end
        end

    end

    return edgelist
    #=
    nz = length(edgelist)
    I = Int[]
    J = Int[]
    for (i,j) in edgelist
        push!(I,i)
        push!(J,j)
    end
    return I,J
    =#

end

export compute_hessian_sparsity

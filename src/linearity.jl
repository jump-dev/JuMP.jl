

# classify the nodes in a tree as constant, linear, or nonlinear with respect to the input

@enum Linearity CONSTANT LINEAR NONLINEAR

export CONSTANT, LINEAR, NONLINEAR

function classify_linearity(nd::Vector{NodeData},adj)

    linearity = Array(Linearity,length(nd))

    # do a forward pass through the graph, which is reverse order of nd

    children_arr = rowvals(adj)

    for k in length(nd):-1:1
        # compute the value of node k
        nod = nd[k]
        if nod.nodetype == VARIABLE
            linearity[k] = LINEAR
        elseif nod.nodetype == VALUE
            linearity[k] = CONSTANT
        else
            @assert nod.nodetype == CALL || nod.nodetype == CALLUNIVAR
            children_idx = nzrange(adj,k)
            num_constant_children = 0
            any_nonlinear = false
            for r in children_idx
                if linearity[children_arr[r]] == NONLINEAR
                    any_nonlinear = true
                    break
                elseif linearity[children_arr[r]] == CONSTANT
                    num_constant_children += 1
                end
            end
            if any_nonlinear
                # if any children are nonlinear, then we're nonlinear
                linearity[k] = NONLINEAR
                continue
            elseif num_constant_children == length(children_idx)
                # if all children are constant, then we're constant
                linearity[k] = CONSTANT
                continue
            end
            # some children are constant and some are linear
            # if operator is nonlinear, then we're nonlinear
            op = nod.index
            if nod.nodetype == CALLUNIVAR
                if univariate_operators[op] == :+ || univariate_operators[op] == :-
                    linearity[k] = LINEAR
                else
                    linearity[k] = NONLINEAR
                end
            else
                # operator with more than 1 argument
                if operators[op] == :+ || operators[op] == :-
                    linearity[k] = LINEAR
                elseif operators[op] == :* && num_constant_children == length(children_idx) - 1
                    linearity[k] = LINEAR
                elseif operators[op] == :/
                    if linearity[children_arr[last(children_idx)]] == CONSTANT
                        # denominator is constant, we're linear
                        linearity[k] = LINEAR
                    else
                        # denominator is linear
                        linearity[k] = NONLINEAR
                    end
                else # all other operators are nonlinear
                    linearity[k] = NONLINEAR
                end
            end
        end
    end

    return linearity

end

export Linearity,classify_linearity

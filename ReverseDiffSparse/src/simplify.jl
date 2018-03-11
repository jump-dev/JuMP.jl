

# returns a new Vector{NodeData} with constants simplified
# basically, lob off all of the children of constant nodes
function simplify_constants(forward_storage::Vector{T},nd::Vector{NodeData},adj,const_values,linearity) where T

    new_nd = NodeData[]
    old_index_to_new_index = fill(-1,length(nd))

    for k in 1:length(nd) # process parents first
        nod = nd[k]
        new_parentidx = nod.parent != -1 ? old_index_to_new_index[nod.parent] : -1
        if new_parentidx != -1 && linearity[nod.parent] == CONSTANT
            # if parent is constant, we skip this
            old_index_to_new_index[k] = 0
            continue
        elseif linearity[k] != CONSTANT || nod.nodetype == VALUE || nod.nodetype == PARAMETER 
            # put it back on the tape
            new_nod = NodeData(nod.nodetype,nod.index,new_parentidx)
            push!(new_nd, new_nod)
            old_index_to_new_index[k] = length(new_nd)
        else 
            # this node is constant and needs to be replaced
            push!(const_values, forward_storage[k])
            idx = length(const_values)
            new_nod = NodeData(VALUE, idx, new_parentidx)
            push!(new_nd, new_nod)
            old_index_to_new_index[k] = length(new_nd)
        end
    end

    return new_nd
end

export simplify_constants 


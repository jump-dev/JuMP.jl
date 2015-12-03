
@enum NodeType CALL VARIABLE VALUE PARAMETER

immutable NodeData
    nodetype::NodeType
    index::Int
    parent::Int
    whichchild::Int # which number child is this node
end

export NodeData

# Only need to store parents, then transpose adjacency matrix!

# for CALL, index is list of operators
# for VARIABLE, index is variable index
# for VALUE, index is into list of constants
# for PARAMETER, index is into list of parameters

const operators = [:+,:-,:*,:sin,:cos,:^]

const operator_to_id = Dict{Symbol,Int}()
for i in 1:length(operators)
    operator_to_id[operators[i]] = i
end



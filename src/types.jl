
@enum NodeType CALL CALLUNIVAR VARIABLE VALUE PARAMETER

immutable NodeData
    nodetype::NodeType
    index::Int
    parent::Int
    whichchild::Int # which number child is this node
end

export NodeData

# Only need to store parents, then transpose adjacency matrix!

# for CALL, index is into list of operators
# for CALLUNIVAR, index is into list of univariate operators
# for VARIABLE, index is variable index
# for VALUE, index is into list of constants
# for PARAMETER, index is into list of parameters

const operators = [:+,:-,:*,:^]

const operator_to_id = Dict{Symbol,Int}()
for i in 1:length(operators)
    operator_to_id[operators[i]] = i
end

const univariate_operators = Symbol[:+,:-]
const univariate_operator_to_id = Dict{Symbol,Int}(:+ => 1, :- => 2)
const univariate_operator_deriv = Any[:x,:(-x)]

for (op, deriv) in Calculus.symbolic_derivatives_1arg()
    push!(univariate_operators, op)
    push!(univariate_operator_deriv,deriv)
    univariate_operator_to_id[op] = length(univariate_operators)
end

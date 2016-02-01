
@enum NodeType CALL CALLUNIVAR VARIABLE VALUE PARAMETER SUBEXPRESSION LOGIC COMPARISON

export CALL, CALLUNIVAR, VARIABLE, VALUE, PARAMETER, SUBEXPRESSION, LOGIC, COMPARISON

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
# for SUBEXPRESSION, index is into list of subexpressions
# for LOGIC, index is into list of logical operators (inputs and outputs are boolean)
# for COMPARISON, index is into lost of comparison operators

const operators = [:+,:-,:*,:^,:/,:ifelse,:max,:min]

const operator_to_id = Dict{Symbol,Int}()
for i in 1:length(operators)
    operator_to_id[operators[i]] = i
end
export operator_to_id, operators

const univariate_operators = Symbol[:+,:-,:abs]
const univariate_operator_to_id = Dict{Symbol,Int}(:+ => 1, :- => 2, :abs => 3)
const univariate_operator_deriv = Any[:(one(x)),:(-one(x)),:(ifelse(x >= 0, one(x),-one(x)))]
export univariate_operator_to_id, univariate_operators

for (op, deriv) in Calculus.symbolic_derivatives_1arg()
    push!(univariate_operators, op)
    push!(univariate_operator_deriv,deriv)
    univariate_operator_to_id[op] = length(univariate_operators)
end

const logic_operators = [:&&,:||]
const logic_operator_to_id = Dict{Symbol,Int}()
for i in 1:length(logic_operators)
    logic_operator_to_id[logic_operators[i]] = i
end
export logic_operator_to_id, logic_operators

const comparison_operators = [:<=,:(==),:>=,:<,:>]
const comparison_operator_to_id = Dict{Symbol,Int}()
for i in 1:length(comparison_operators)
    comparison_operator_to_id[comparison_operators[i]] = i
end
export comparison_operator_to_id, comparison_operators

# accumulator for prod{} terms
# we need to handle this specially so that we can compute
# derivatives when one of the terms is zero.
# (If two are zero, all derivatives are zero.)
immutable ProductAccumulator{T}
    allbutone::T # product of all but the smallest value
    smallest::T
end

ProductAccumulator{T}(::Type{T}) = ProductAccumulator(one(T),one(T))

@inline function add_term{T}(p::ProductAccumulator{T},value::T)
    if -1e-4 <= value <= 1e-4 && abs(value) < abs(p.smallest)
        return ProductAccumulator(p.allbutone*p.smallest,value)
    else
        return ProductAccumulator(p.allbutone*value,p.smallest)
    end
end

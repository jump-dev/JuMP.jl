
@enum NodeType CALL CALLUNIVAR MOIVARIABLE VARIABLE VALUE PARAMETER SUBEXPRESSION LOGIC COMPARISON EXTRA

export CALL,
    CALLUNIVAR,
    MOIVARIABLE,
    VARIABLE,
    VALUE,
    PARAMETER,
    SUBEXPRESSION,
    LOGIC,
    COMPARISON

struct NodeData
    nodetype::NodeType
    index::Int64
    parent::Int64
end

export NodeData

# Only need to store parents, then transpose adjacency matrix!

# for CALL, index is into list of operators
# for CALLUNIVAR, index is into list of univariate operators
# for MOIVARIABLE, index is MOI variable index (not consecutive)
# for VARIABLE, index is internal (consecutive) variable index (see variable_order in MOI.initialize!)
# for VALUE, index is into list of constants
# for PARAMETER, index is into list of parameters
# for SUBEXPRESSION, index is into list of subexpressions
# for LOGIC, index is into list of logical operators (inputs and outputs are boolean)
# for COMPARISON, index is into lost of comparison operators
# for EXTRA, index is extension specific

const operators = [:+, :-, :*, :^, :/, :ifelse, :max, :min]
const USER_OPERATOR_ID_START = length(operators) + 1

const operator_to_id = Dict{Symbol,Int}()
for i in 1:length(operators)
    operator_to_id[operators[i]] = i
end
export operator_to_id, operators

const univariate_operators = Symbol[:+, :-, :abs]
const univariate_operator_to_id = Dict{Symbol,Int}(:+ => 1, :- => 2, :abs => 3)
const univariate_operator_deriv =
    Any[:(one(x)), :(-one(x)), :(ifelse(x >= 0, one(x), -one(x)))]
export univariate_operator_to_id, univariate_operators

for (op, deriv) in Calculus.symbolic_derivatives_1arg()
    push!(univariate_operators, op)
    push!(univariate_operator_deriv, deriv)
    univariate_operator_to_id[op] = length(univariate_operators)
end
const USER_UNIVAR_OPERATOR_ID_START = length(univariate_operators) + 1

const logic_operators = [:&&, :||]
const logic_operator_to_id = Dict{Symbol,Int}()
for i in 1:length(logic_operators)
    logic_operator_to_id[logic_operators[i]] = i
end
export logic_operator_to_id, logic_operators

const comparison_operators = [:<=, :(==), :>=, :<, :>]
const comparison_operator_to_id = Dict{Symbol,Int}()
for i in 1:length(comparison_operators)
    comparison_operator_to_id[comparison_operators[i]] = i
end
export comparison_operator_to_id, comparison_operators

# user-provided operators
struct UserOperatorRegistry
    multivariate_operator_to_id::Dict{Symbol,Int}
    multivariate_operator_evaluator::Vector{MOI.AbstractNLPEvaluator}
    univariate_operator_to_id::Dict{Symbol,Int}
    univariate_operator_f::Vector{Any}
    univariate_operator_fprime::Vector{Any}
    univariate_operator_fprimeprime::Vector{Any}
end

function UserOperatorRegistry()
    return UserOperatorRegistry(
        Dict{Symbol,Int}(),
        MOI.AbstractNLPEvaluator[],
        Dict{Symbol,Int}(),
        [],
        [],
        [],
    )
end

# we use the MathOptInterface NLPEvaluator interface, where the
# operator takes the place of the objective function.
# users should implement eval_f and eval_grad_f for now.
# we will eventually support hessians too
function register_multivariate_operator!(
    r::UserOperatorRegistry,
    s::Symbol,
    f::MOI.AbstractNLPEvaluator,
)
    if haskey(r.multivariate_operator_to_id, s)
        error("""The multivariate operator $s has already been defined.

        User-defined functions can not be re-registered. If you want to
        re-register a function between solves (e.g., to modify the user-defined
        function), rebuild the model or use a different function name.""")
    end
    id = length(r.multivariate_operator_evaluator) + 1
    r.multivariate_operator_to_id[s] = id
    push!(r.multivariate_operator_evaluator, f)
    return
end

export register_multivariate_operator!

# for univariate operators, just take in functions to evaluate
# zeroth, first, and second order derivatives
function register_univariate_operator!(
    r::UserOperatorRegistry,
    s::Symbol,
    f,
    fprime,
    fprimeprime,
)
    if haskey(r.univariate_operator_to_id, s)
        error("""The univariate operator $s has already been defined.

        User-defined functions can not be re-registered. If you want to
        re-register a function between solves (e.g., to modify the user-defined
        function), rebuild the model or use a different function name.""")
    end
    id = length(r.univariate_operator_f) + 1
    r.univariate_operator_to_id[s] = id
    push!(r.univariate_operator_f, f)
    push!(r.univariate_operator_fprime, fprime)
    push!(r.univariate_operator_fprimeprime, fprimeprime)
    return
end

export register_univariate_operator!

function has_user_multivariate_operators(nd::Vector{NodeData})
    for k in 1:length(nd)
        nod = nd[k]
        if nod.nodetype == CALL && nod.index >= USER_OPERATOR_ID_START
            return true
        end
    end
    return false
end

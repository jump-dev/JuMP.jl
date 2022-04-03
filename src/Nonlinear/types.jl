#  Copyright 2017, Iain Dunning, Joey Huchette, Miles Lubin, and contributors
#  This Source Code Form is subject to the terms of the Mozilla Public
#  License, v. 2.0. If a copy of the MPL was not distributed with this
#  file, You can obtain one at https://mozilla.org/MPL/2.0/.

@enum(
    NodeType,
    # Index into multivariate operators
    NODE_CALL_MULTIVARIATE,
    # Index into univariate operators
    NODE_CALL_UNIVARIATE,
    # Index into logic operators
    NODE_LOGIC,
    # Index into comparison operators
    NODE_COMPARISON,
    # Index is the value of `MOI.VariableIndex`. This is from the original
    # model, and is not consecutive.
    NODE_MOI_VARIABLE,
    # Index of the internal, consecutive, and ordered `MOI.VariableIndex`.
    NODE_VARIABLE,
    # Index is into the list of constants
    NODE_VALUE,
    # Index is into the list of parameters
    NODE_PARAMETER,
    # Index is into the list of subexpressions
    NODE_SUBEXPRESSION,
)

struct Node
    type::NodeType
    index::Int
    parent::Int
end

struct NonlinearExpression
    nodes::Vector{Node}
    values::Vector{Float64}
    NonlinearExpression() = new(Node[], Float64[])
end

function Base.:(==)(x::NonlinearExpression, y::NonlinearExpression)
    return x.nodes == y.nodes && x.values == y.values
end

struct NonlinearConstraint
    expression::NonlinearExpression
    set::Union{
        MOI.LessThan{Float64},
        MOI.GreaterThan{Float64},
        MOI.EqualTo{Float64},
        MOI.Interval{Float64},
    }
end

"""
    ParameterIndex
"""
struct ParameterIndex
    value::Int
end

"""
    ExpressionIndex
"""
struct ExpressionIndex
    value::Int
end

"""
    ConstraintIndex
"""
struct ConstraintIndex
    value::Int
end

"""
    NonlinearData()

Construct a new [`MOI.AbstractNLPEvaluator`](@ref) object.
"""
mutable struct NonlinearData <: MOI.AbstractNLPEvaluator
    inner::Union{Nothing,MOI.AbstractNLPEvaluator}
    objective::Union{Nothing,NonlinearExpression}
    expressions::Vector{NonlinearExpression}
    constraints::OrderedDict{ConstraintIndex,NonlinearConstraint}
    parameters::Vector{Float64}
    operators::OperatorRegistry
    last_constraint_index::Int64
    # Fields for initialize
    julia_expressions::Vector{Any}  # Any because expressions may be constants
    ordered_constraints::Vector{ConstraintIndex}
    constraint_dual::Vector{Float64}
    # Timers
    eval_objective_timer::Float64
    eval_constraint_timer::Float64
    eval_objective_gradient_timer::Float64
    eval_constraint_jacobian_timer::Float64
    eval_hessian_lagrangian_timer::Float64

    function NonlinearData()
        return new(
            nothing,
            nothing,
            NonlinearExpression[],
            OrderedDict{ConstraintIndex,NonlinearConstraint}(),
            Float64[],
            OperatorRegistry(),
            0,
            Expr[],
            ConstraintIndex[],
            Float64[],
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
        )
    end
end

"""
    AbstractAutomaticDifferentiation

An abstract type for extending
"""
abstract type AbstractAutomaticDifferentiation end

"""
    set_differentiation_backend(
        data::NonlinearData,
        backend::AbstractAutomaticDifferentiation,
        ordered_variables::Vector{MOI.VariableIndex}
    )

Set the automatic differentiationn backend of `data` to `backend`.

`ordered_variables` is an ordered vector containing all of the variables that
are present in the model. This order corresponds to the order of the primal
solution vector `x` that is passed to the various functions in MOI's nonlinear
API.
"""
function set_differentiation_backend end

"""
    Default() <: AbstractAutomaticDifferentiation

The default implementation of `AbstractAutomaticDifferentiation`. This only
supported `:ExprGraph`.
"""
struct Default <: AbstractAutomaticDifferentiation end

function set_differentiation_backend(
    data::NonlinearData,
    ::Default,
    ::Vector{MOI.VariableIndex},
)
    data.inner = nothing
    return
end

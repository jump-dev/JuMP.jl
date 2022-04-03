#  Copyright 2017, Iain Dunning, Joey Huchette, Miles Lubin, and contributors
#  This Source Code Form is subject to the terms of the Mozilla Public
#  License, v. 2.0. If a copy of the MPL was not distributed with this
#  file, You can obtain one at https://mozilla.org/MPL/2.0/.

function parse_expression(data::NonlinearData, input)
    expr = NonlinearExpression()
    _parse_expression(data, expr, input, -1)
    return expr
end

function _parse_expression(
    data::NonlinearData,
    expr::NonlinearExpression,
    x::Expr,
    parent_index::Int,
)
    if isexpr(x, :call)
        if length(x.args) == 2 && !isexpr(x.args[2], :...)
            _parse_univariate_expression(data, expr, x, parent_index)
        else
            _parse_multivariate_expression(data, expr, x, parent_index)
        end
    elseif isexpr(x, :comparison)
        _parse_comparison_expression(data, expr, x, parent_index)
    elseif isexpr(x, :...)
        _parse_splat_expression(data, expr, x, parent_index)
    elseif isexpr(x, :&&) || isexpr(x, :||)
        _parse_logic_expression(data, expr, x, parent_index)
    else
        error("Unsupported expression: $x")
    end
    return
end

function _parse_splat_expression(data, expr, x, parent_index)
    @assert isexpr(x, :...) && length(x.args) == 1
    if parent_index == -1
        error(
            "Unsupported use of the splatting operator. This is only " *
            "supported in the arguments of a function call.",
        )
    elseif x.args[1] isa Expr
        error(
            "Unsupported use of the splatting operator. JuMP supports " *
            "splatting only symbols. For example, `x...` is ok, but " *
            "`(x + 1)...`, `[x; y]...` and `g(f(y)...)` are not.",
        )
    end
    for xi in x.args[1]
        _parse_expression(data, expr, xi, parent_index)
    end
    return
end

function _parse_univariate_expression(
    data::NonlinearData,
    expr::NonlinearExpression,
    x::Expr,
    parent_index::Union{Symbol,Int},
)
    @assert isexpr(x, :call, 2)
    id = get(data.operators.univariate_operator_to_id, x.args[1], nothing)
    if id === nothing
        # It may also be a multivariate operator like * with one argument.
        if haskey(data.operators.multivariate_operator_to_id, x.args[1])
            _parse_multivariate_expression(data, expr, x, parent_index)
            return
        end
        error("Unable to parse: $x")
    end
    push!(expr.nodes, Node(NODE_CALL_UNIVARIATE, id, parent_index))
    _parse_expression(data, expr, x.args[2], length(expr.nodes))
    return
end

function _parse_multivariate_expression(
    data::NonlinearData,
    expr::NonlinearExpression,
    x::Expr,
    parent_index::Union{Symbol,Int},
)
    @assert isexpr(x, :call)
    id = get(data.operators.multivariate_operator_to_id, x.args[1], nothing)
    if id === nothing
        @assert x.args[1] in data.operators.comparison_operators
        _parse_inequality_expression(data, expr, x, parent_index)
        return
    end
    push!(expr.nodes, Node(NODE_CALL_MULTIVARIATE, id, parent_index))
    parent_var = length(expr.nodes)
    for i in 2:length(x.args)
        _parse_expression(data, expr, x.args[i], parent_var)
    end
    return
end

# This function parses single inequalities like `a <= b`. It's not to be
# confused with `_parse_comparison_expression`, which handles things like
# `a <= b <= c`.
function _parse_inequality_expression(
    data::NonlinearData,
    expr::NonlinearExpression,
    x::Expr,
    parent_index::Union{Symbol,Int},
)
    operator_id = data.operators.comparison_operator_to_id[x.args[1]]
    push!(expr.nodes, Node(NODE_COMPARISON, operator_id, parent_index))
    parent_var = length(expr.nodes)
    for i in 2:length(x.args)
        _parse_expression(data, expr, x.args[i], parent_var)
    end
    return
end

# This function parses double inequalities like `a <= b <= c`. It's not to be
# confused with `_parse_inequality_expression`, which handles things like
# `a <= b`.
function _parse_comparison_expression(
    data::NonlinearData,
    expr::NonlinearExpression,
    x::Expr,
    parent_index::Union{Symbol,Int},
)
    for k in 2:2:length(x.args)-1
        @assert x.args[k] == x.args[2] # don't handle a <= b >= c
    end
    operator_id = data.operators.comparison_operator_to_id[x.args[2]]
    push!(expr.nodes, Node(NODE_COMPARISON, operator_id, parent_index))
    parent_var = length(expr.nodes)
    for i in 1:2:length(x.args)
        _parse_expression(data, expr, x.args[i], parent_var)
    end
    return
end

function _parse_logic_expression(
    data::NonlinearData,
    expr::NonlinearExpression,
    x::Expr,
    parent_index::Union{Symbol,Int},
)
    id = data.operators.logic_operator_to_id[x.head]
    push!(expr.nodes, Node(NODE_LOGIC, id, parent_index))
    parent_var = length(expr.nodes)
    _parse_expression(data, expr, x.args[1], parent_var)
    _parse_expression(data, expr, x.args[2], parent_var)
    return
end

function _parse_expression(
    ::NonlinearData,
    expr::NonlinearExpression,
    x::MOI.VariableIndex,
    parent_index::Int,
)
    push!(expr.nodes, Node(NODE_MOI_VARIABLE, x.value, parent_index))
    return
end

function _parse_expression(
    ::NonlinearData,
    expr::NonlinearExpression,
    x::Real,
    parent_index::Int,
)
    push!(expr.values, convert(Float64, x)::Float64)
    push!(expr.nodes, Node(NODE_VALUE, length(expr.values), parent_index))
    return
end

function _parse_expression(
    ::NonlinearData,
    ::NonlinearExpression,
    x::Any,
    ::Int,
)
    return error(
        "Unexpected object $x of type $(typeof(x)) in nonlinear expression.",
    )
end

function _parse_expression(
    ::NonlinearData,
    expr::NonlinearExpression,
    x::ParameterIndex,
    parent_index::Int,
)
    push!(expr.nodes, Node(NODE_PARAMETER, x.value, parent_index))
    return
end

function _parse_expression(
    ::NonlinearData,
    expr::NonlinearExpression,
    x::ExpressionIndex,
    parent_index::Int,
)
    push!(expr.nodes, Node(NODE_SUBEXPRESSION, x.value, parent_index))
    return
end

function _parse_expression(
    ::NonlinearData,
    expr::NonlinearExpression,
    x::AbstractArray,
    parent_index::Int,
)
    return error(
        "Unexpected array $x in nonlinear expression. Nonlinear expressions " *
        "may contain only scalar expressions.",
    )
end

function _normalize_constraint_expr(lhs::Real, body, rhs::Real)
    return Float64(lhs), body, Float64(rhs)
end

function _normalize_constraint_expr(lhs, body, rhs)
    return error(
        "Interval constraint contains non-constant left- or right-hand " *
        "sides. Reformulate as two separate constraints, or move all " *
        "variables into the central term.",
    )
end

_normalize_constraint_expr(lhs, rhs::Real) = lhs, Float64(rhs)

_normalize_constraint_expr(lhs, rhs) = Expr(:call, :-, lhs, rhs), 0.0

function _expr_to_constraint(expr::Expr)
    if isexpr(expr, :comparison)
        @assert expr.args[2] == expr.args[4]
        @assert expr.args[2] in (:<=, :>=)
        lhs, body, rhs =
            _normalize_constraint_expr(expr.args[1], expr.args[3], expr.args[5])
        return body, MOI.Interval(lhs, rhs)
    end
    lhs, rhs = _normalize_constraint_expr(expr.args[2], expr.args[3])
    if expr.args[1] == :<=
        return :($lhs - $rhs), MOI.LessThan(0.0)
    elseif expr.args[1] == :>=
        return :($lhs - $rhs), MOI.GreaterThan(0.0)
    else
        @assert expr.args[1] == :(==)
        return :($lhs - $rhs), MOI.EqualTo(0.0)
    end
end

function _to_expr(
    data::NonlinearData,
    expr::NonlinearExpression;
    expand_subexpressions::Bool = true,
    use_x_ref::Bool = false,
)
    tree = Any[]
    for node in expr.nodes
        node_expr = if node.type == NODE_CALL_UNIVARIATE
            Expr(:call, data.operators.univariate_operators[node.index])
        elseif node.type == NODE_CALL_MULTIVARIATE
            Expr(:call, data.operators.multivariate_operators[node.index])
        elseif node.type == NODE_COMPARISON
            Expr(:call, data.operators.comparison_operators[node.index])
        elseif node.type == NODE_LOGIC
            Expr(data.operators.logic_operators[node.index])
        elseif node.type == NODE_MOI_VARIABLE
            x = MOI.VariableIndex(node.index)
            use_x_ref ? :(x[$x]) : x
        elseif node.type == NODE_PARAMETER
            if expand_subexpressions
                data.parameters[node.index]
            else
                ParameterIndex(node.index)
            end
        elseif node.type == NODE_SUBEXPRESSION
            if expand_subexpressions
                data.julia_expressions[node.index]
            else
                ExpressionIndex(node.index)
            end
        else
            @assert node.type == NODE_VALUE
            expr.values[node.index]
        end
        if 1 <= node.parent <= length(tree)
            push!(tree[node.parent].args, node_expr)
        end
        push!(tree, node_expr)
    end
    return _replace_comparison(tree[1])
end

# TODO(odow): NODE_COMPARISON is a bit of a pain to deal with...
_replace_comparison(expr) = expr

function _replace_comparison(expr::Expr)
    if isexpr(expr, :call, 4) && expr.args[1] in (:<=, :>=, :(==), :<, :>)
        return Expr(
            :comparison,
            expr.args[2],
            expr.args[1],
            expr.args[3],
            expr.args[1],
            expr.args[4],
        )
    end
    for i in 1:length(expr.args)
        expr.args[i] = _replace_comparison(expr.args[i])
    end
    return expr
end

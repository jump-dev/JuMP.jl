#  Copyright 2017, Iain Dunning, Joey Huchette, Miles Lubin, and contributors
#  This Source Code Form is subject to the terms of the Mozilla Public
#  License, v. 2.0. If a copy of the MPL was not distributed with this
#  file, You can obtain one at https://mozilla.org/MPL/2.0/.

struct NonlinearExpr <: AbstractJuMPScalar
    head::Symbol
    args::Vector{Any}
end

# We include this method so that we can refactor the internal representation of
# NonlinearExpr without having to rewrite the method overloads.
function NonlinearExpr(head::Symbol, args...)
    return NonlinearExpr(head, Any[args...])
end

Base.length(x::NonlinearExpr) = length(x.args)
Base.getindex(x::NonlinearExpr, i::Int) = x.args[i]

function function_string(::MIME"text/plain", x::NonlinearExpr)
    io, stack, is_open = IOBuffer(), Any[x], true
    while !isempty(stack)
        arg = pop!(stack)
        if !is_open && arg != ')'
            print(io, ", ")
        end
        if arg isa NonlinearExpr
            print(io, arg.head, "(")
            push!(stack, ')')
            for i in length(arg):-1:1
                push!(stack, arg[i])
            end
        else
            print(io, arg)
        end
        is_open = arg isa NonlinearExpr
    end
    seekstart(io)
    return read(io, String)
end

function function_string(::MIME"text/latex", expr::NonlinearExpr)
    return "\\textsf{$(function_string(MIME("text/plain"), expr))}"
end

_isequal(x, y) = x == y
_isequal(x::T, y::T) where {T<:AbstractJuMPScalar} = isequal_canonical(x, y)

function isequal_canonical(x::NonlinearExpr, y::NonlinearExpr)
    return x.head == y.head &&
           length(x) == length(y) &&
           all(i -> _isequal(x[i], y[i]), 1:length(x))
end

function MOI.Nonlinear.parse_expression(
    data::MOI.Nonlinear.Model,
    expr::MOI.Nonlinear.Expression,
    x::NonlinearExpr,
    parent::Int,
)
    stack = Tuple{Int,Any}[(parent, x)]
    while !isempty(stack)
        parent_node, arg = pop!(stack)
        if arg isa NonlinearExpr
            _parse_without_recursion_inner(stack, data, expr, arg, parent_node)
        else
            # We can use recursion here, because NonlinearExpr only occur in
            # other NonlinearExpr.
            MOI.Nonlinear.parse_expression(data, expr, arg, parent_node)
        end
    end
    return
end

function _parse_without_recursion_inner(stack, data, expr, x, parent)
    id = get(data.operators.univariate_operator_to_id, x.head, nothing)
    node_type = if length(x) == 1 && id !== nothing
        MOI.Nonlinear.NODE_CALL_UNIVARIATE
    else
        id = get(data.operators.multivariate_operator_to_id, x.head, nothing)
        @assert id !== nothing
        MOI.Nonlinear.NODE_CALL_MULTIVARIATE
    end
    push!(expr.nodes, MOI.Nonlinear.Node(node_type, id, parent))
    parent = length(expr.nodes)
    for i in length(x):-1:1  # Args need to be pushed onto the stack in reverse
        push!(stack, (parent, x[i]))
    end
    return
end

# Method definitions

Base.zero(::Type{NonlinearExpr}) = NonlinearExpr(:+, 0.0)

Base.one(::Type{NonlinearExpr}) = NonlinearExpr(:+, 1.0)

# Univariate operators

for f in MOI.Nonlinear.DEFAULT_UNIVARIATE_OPERATORS
    op = Meta.quot(f)
    if f == :+
        continue  # We don't need this.
    elseif f == :-
        @eval Base.:-(x::NonlinearExpr) = NonlinearExpr(:-, x)
    elseif isdefined(Base, f)
        @eval Base.$(f)(x::AbstractJuMPScalar) = NonlinearExpr($op, x)
    elseif isdefined(MOI.Nonlinear, :SpecialFunctions)
        # The operator is defined in some other package.
        SF = MOI.Nonlinear.SpecialFunctions
        if isdefined(SF, f)
            @eval $(SF).$(f)(x::AbstractJuMPScalar) = NonlinearExpr($op, x)
        end
    end
end

# Multivariate operators

# The multivariate operators in MOI are +, -, *, ^, /, ifelse, atan
#
# However, ifelse is a builtin, so we can't add methods to it.

# We need only very generic fallbacks for these, because all other cases are
# caught with more specific methods.
for f in (:+, :-, :*, :^, :/, :atan)
    op = Meta.quot(f)
    @eval begin
        function Base.$(f)(x::AbstractJuMPScalar, y::_Constant)
            rhs = convert(Float64, _constant_to_number(y))
            return NonlinearExpr($op, x, rhs)
        end
        function Base.$(f)(x::_Constant, y::AbstractJuMPScalar)
            lhs = convert(Float64, _constant_to_number(x))
            return NonlinearExpr($op, lhs, y)
        end
        function Base.$(f)(x::AbstractJuMPScalar, y::AbstractJuMPScalar)
            return NonlinearExpr($op, x, y)
        end
    end
end

# JuMP interop

# TODO
check_belongs_to_model(::NonlinearExpr, ::Model) = true

function moi_function(f::NonlinearExpr)
    ret = MOI.ScalarNonlinearFunction{Float64}(f.head, Any[])
    stack = Tuple{MOI.ScalarNonlinearFunction{Float64},Any}[]
    for arg in reverse(f.args)
        push!(stack, (ret, arg))
    end
    while !isempty(stack)
        parent, arg = pop!(stack)
        if arg isa NonlinearExpr
            new_ret = MOI.ScalarNonlinearFunction{Float64}(arg.head, Any[])
            push!(parent.args, new_ret)
            for child in reverse(arg.args)
                push!(stack, (new_ret, child))
            end
        elseif arg isa Number
            push!(parent.args, arg)
        else
            push!(parent.args, moi_function(arg))
        end
    end
    return ret
end

function jump_function(model::Model, f::MOI.ScalarNonlinearFunction)
    ret = NonlinearExpr(f.head, Any[])
    stack = Tuple{NonlinearExpr,Any}[]
    for arg in reverse(f.args)
        push!(stack, (ret, arg))
    end
    while !isempty(stack)
        parent, arg = pop!(stack)
        if arg isa MOI.ScalarNonlinearFunction
            new_ret = NonlinearExpr(arg.head, Any[])
            push!(parent.args, new_ret)
            for child in reverse(arg.args)
                push!(stack, (new_ret, child))
            end
        elseif arg isa Number
            push!(parent.args, arg)
        else
            push!(parent.args, jump_function(model, arg))
        end
    end
    return ret
end

function jump_function_type(::Model, ::Type{<:MOI.ScalarNonlinearFunction})
    return NonlinearExpr
end

function moi_function_type(::Type{NonlinearExpr})
    return MOI.ScalarNonlinearFunction{Float64}
end

function constraint_object(c::NonlinearConstraintRef)
    nlp = nonlinear_model(c.model)
    data = nlp.constraints[index(c)]
    return ScalarConstraint(jump_function(c.model, data.expression), data.set)
end

function jump_function(model::Model, expr::MOI.Nonlinear.Expression)
    nlp = nonlinear_model(model)
    parsed = Vector{Any}(undef, length(expr.nodes))
    adj = MOI.Nonlinear.adjacency_matrix(expr.nodes)
    rowvals = SparseArrays.rowvals(adj)
    for i in length(expr.nodes):-1:1
        node = expr.nodes[i]
        parsed[i] = if node.type == MOI.Nonlinear.NODE_CALL_UNIVARIATE
            NonlinearExpr(
                nlp.operators.univariate_operators[node.index],
                parsed[rowvals[SparseArrays.nzrange(adj, i)[1]]],
            )
        elseif node.type == MOI.Nonlinear.NODE_CALL_MULTIVARIATE
            NonlinearExpr(
                nlp.operators.multivariate_operators[node.index],
                Any[parsed[rowvals[j]] for j in SparseArrays.nzrange(adj, i)],
            )
        elseif node.type == MOI.Nonlinear.NODE_MOI_VARIABLE
            VariableRef(model, MOI.VariableIndex(node.index))
        elseif node.type == MOI.Nonlinear.NODE_PARAMETER
            NonlinearParameter(model, node.index)
        elseif node.type == MOI.Nonlinear.NODE_SUBEXPRESSION
            NonlinearExpression(model, node.index)
        elseif node.type == MOI.Nonlinear.NODE_VALUE
            expr.values[node.index]
        else
            # node.type == MOI.Nonlinear.NODE_COMPARISON
            # node.type == MOI.Nonlinear.NODE_LOGIC
            error("Unsupported node")
        end
    end
    return parsed[1]
end

# MutableArithmetics.jl

# These converts are used in the {add,sub}mul definition for AbstractJuMPScalar.

Base.convert(::Type{NonlinearExpr}, x::AbstractVariableRef) = x

function Base.convert(::Type{NonlinearExpr}, x::GenericAffExpr)
    args = Any[]
    for (variable, coef) in x.terms
        if isone(coef)
            push!(args, variable)
        elseif !iszero(coef)
            push!(args, NonlinearExpr(:*, coef, variable))
        end
    end
    if !iszero(x.constant) || isempty(args)
        push!(args, x.constant)
    end
    if length(args) == 1
        return args[1]
    end
    return NonlinearExpr(:+, args)
end

function Base.convert(::Type{NonlinearExpr}, x::GenericQuadExpr)
    args = Any[]
    for (variable, coef) in x.aff.terms
        if isone(coef)
            push!(args, variable)
        elseif !iszero(coef)
            push!(args, NonlinearExpr(:*, coef, variable))
        end
    end
    for (pair, coef) in x.terms
        if isone(coef)
            push!(args, NonlinearExpr(:*, pair.a, pair.b))
        elseif !iszero(coef)
            push!(args, NonlinearExpr(:*, coef, pair.a, pair.b))
        end
    end
    if !iszero(x.aff.constant) || isempty(args)
        push!(args, x.aff.constant)
    end
    if length(args) == 1
        return args[1]
    end
    return NonlinearExpr(:+, args)
end

function _MA.promote_operation(
    ::Union{typeof(+),typeof(-),typeof(*)},
    ::Type{NonlinearExpr},
    ::Type{<:AbstractJuMPScalar},
)
    return NonlinearExpr
end

function _MA.promote_operation(
    ::Union{typeof(+),typeof(-),typeof(*)},
    ::Type{<:AbstractJuMPScalar},
    ::Type{NonlinearExpr},
)
    return NonlinearExpr
end

function _MA.promote_operation(
    ::Union{typeof(+),typeof(-),typeof(*)},
    ::Type{NonlinearExpr},
    ::Type{NonlinearExpr},
)
    return NonlinearExpr
end

function moi_function(x::Vector{<:AbstractJuMPScalar})
    return MOI.VectorNonlinearFunction{Float64}(
        Any[moi_function(xi) for xi in x],
    )
end

function moi_function_type(::Type{<:Vector{<:AbstractJuMPScalar}})
    return MOI.VectorNonlinearFunction{Float64}
end

function jump_function(model::Model, f::MOI.VectorNonlinearFunction)
    return AbstractJuMPScalar[jump_function(model, arg) for arg in f.args]
end

function jump_function_type(::Model, ::Type{<:MOI.VectorNonlinearFunction})
    return Vector{AbstractJuMPScalar}
end

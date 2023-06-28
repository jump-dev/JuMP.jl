#  Copyright 2017, Iain Dunning, Joey Huchette, Miles Lubin, and contributors
#  This Source Code Form is subject to the terms of the Mozilla Public
#  License, v. 2.0. If a copy of the MPL was not distributed with this
#  file, You can obtain one at https://mozilla.org/MPL/2.0/.

"""
    GenericNonlinearExpr{V}(head::Symbol, args::Vector{Any})
    GenericNonlinearExpr{V}(head::Symbol, args::Any...)

The scalar-valued nonlinear function `head(args...)`, represented as a symbolic
expression tree, with the call operator `head` and ordered arguments in `args`.

`V` is the type of [`AbstractVariableRef`](@ref) present in the expression, and
is used to help dispatch JuMP extensions.

## `head`

The `head::Symbol` must be an operator supported by the model.

The default list of supported univariate operators is given by:

 * [`MOI.Nonlinear.DEFAULT_UNIVARIATE_OPERATORS`](@ref)

and the default list of supported multivariate operators is given by:

 * [`MOI.Nonlinear.DEFAULT_MULTIVARIATE_OPERATORS`](@ref)

Additional operators can be registered by setting a [`MOI.UserDefinedFunction`](@ref)
attribute.

See the full list of operators supported by a [`MOI.ModelLike`](@ref) by
querying [`MOI.ListOfSupportedNonlinearOperators`](@ref).

## `args`

The vector `args` contains the arguments to the nonlinear function. If the
operator is univariate, it must contain one element. Otherwise, it may contain
multiple elements.

Given a subtype of [`AbstractVariableRef`](@ref), `V`, for `GenericNonlinearExpr{V}`,
each element must be one of the following:

 * A constant value of type `<:Number`
 * A `V`
 * A [`GenericAffExpr{C,V}`](@ref)
 * A [`GenericQuadExpr{C,V}`](@ref)
 * A [`GenericNonlinearExpr{V}`](@ref)

## Unsupported operators

If the optimizer does not support `head`, an [`MOI.UnsupportedNonlinearOperator`](@ref)
error will be thrown.

There is no guarantee about when this error will be thrown; it may be thrown
when the function is first added to the model, or it may be thrown when
[`optimize!`](@ref) is called.

## Example

To represent the function ``f(x) = sin(x)^2``, do:

```jldoctest
julia> model = Model();

julia> @variable(model, x)
x

julia> f = sin(x)^2
sin(x) ^ 2.0

julia> f = GenericNonlinearExpr{VariableRef}(
           :^,
           GenericNonlinearExpr{VariableRef}(:sin, x),
           2.0,
       )
sin(x) ^ 2.0
```
"""
struct GenericNonlinearExpr{V<:AbstractVariableRef} <: AbstractJuMPScalar
    head::Symbol
    args::Vector{Any}

    function GenericNonlinearExpr(head::Symbol, args::Vector{Any})
        index = findfirst(Base.Fix2(isa, AbstractJuMPScalar), args)
        if index === nothing
            error(
                "Unable to create a nonlinear expression because it did not " *
                "contain any JuMP scalars. head = $head, args = $args.",
            )
        end
        return new{variable_ref_type(args[index])}(head, args)
    end

    function GenericNonlinearExpr{V}(
        head::Symbol,
        args::Vector{Any},
    ) where {V<:AbstractVariableRef}
        return new{V}(head, args)
    end
end

"""
    NonlinearExpr

Alias for `GenericNonlinearExpr{VariableRef}`, the specific
[`GenericNonlinearExpr`](@ref) used by JuMP.
"""
const NonlinearExpr = GenericNonlinearExpr{VariableRef}

variable_ref_type(::GenericNonlinearExpr{V}) where {V} = V

# We include this method so that we can refactor the internal representation of
# GenericNonlinearExpr without having to rewrite the method overloads.
function GenericNonlinearExpr{V}(
    head::Symbol,
    args...,
) where {V<:AbstractVariableRef}
    return GenericNonlinearExpr{V}(head, Any[args...])
end

Base.length(x::GenericNonlinearExpr) = length(x.args)
Base.getindex(x::GenericNonlinearExpr, i::Int) = x.args[i]

const _PREFIX_OPERATORS =
    (:+, :-, :*, :/, :^, :||, :&&, :>, :<, :(<=), :(>=), :(==))

_needs_parentheses(::Union{Number,AbstractVariableRef}) = false
_needs_parentheses(::Any) = true
function _needs_parentheses(x::GenericNonlinearExpr)
    return x.head in _PREFIX_OPERATORS && length(x) > 1
end

function function_string(::MIME"text/plain", x::GenericNonlinearExpr)
    io, stack = IOBuffer(), Any[x]
    while !isempty(stack)
        arg = pop!(stack)
        if arg isa GenericNonlinearExpr
            if arg.head in _PREFIX_OPERATORS && length(arg) > 1
                if _needs_parentheses(arg[1])
                    print(io, "(")
                end
                if _needs_parentheses(arg.args[end])
                    push!(stack, ")")
                end
                for i in length(arg):-1:2
                    push!(stack, arg[i])
                    if _needs_parentheses(arg.args[i])
                        push!(stack, "(")
                    end
                    push!(stack, " $(arg.head) ")
                    if _needs_parentheses(arg.args[i-1])
                        push!(stack, ")")
                    end
                end
                push!(stack, arg[1])
            else
                print(io, arg.head, "(")
                push!(stack, ")")
                for i in length(arg):-1:2
                    push!(stack, arg[i])
                    push!(stack, ", ")
                end
                push!(stack, arg[1])
            end
        else
            print(io, arg)
        end
    end
    seekstart(io)
    return read(io, String)
end

function function_string(::MIME"text/latex", x::GenericNonlinearExpr)
    io, stack = IOBuffer(), Any[x]
    while !isempty(stack)
        arg = pop!(stack)
        if arg isa GenericNonlinearExpr
            if arg.head in _PREFIX_OPERATORS && length(arg) > 1
                print(io, "\\left({")
                push!(stack, "}\\right)")
                for i in length(arg):-1:2
                    push!(stack, arg[i])
                    push!(stack, "} $(arg.head) {")
                end
                push!(stack, arg[1])
            else
                print(io, "\\textsf{", arg.head, "}\\left({")
                push!(stack, "}\\right)")
                for i in length(arg):-1:2
                    push!(stack, arg[i])
                    push!(stack, "}, {")
                end
                push!(stack, arg[1])
            end
        else
            print(io, arg)
        end
    end
    seekstart(io)
    return read(io, String)
end

_isequal(x, y) = x == y
_isequal(x::T, y::T) where {T<:AbstractJuMPScalar} = isequal_canonical(x, y)

function isequal_canonical(x::GenericNonlinearExpr, y::GenericNonlinearExpr)
    return x.head == y.head &&
           length(x) == length(y) &&
           all(i -> _isequal(x[i], y[i]), 1:length(x))
end

function MOI.Nonlinear.parse_expression(
    data::MOI.Nonlinear.Model,
    expr::MOI.Nonlinear.Expression,
    x::GenericNonlinearExpr,
    parent::Int,
)
    stack = Tuple{Int,Any}[(parent, x)]
    while !isempty(stack)
        parent_node, arg = pop!(stack)
        if arg isa GenericNonlinearExpr
            _parse_without_recursion_inner(stack, data, expr, arg, parent_node)
        else
            # We can use recursion here, because GenericNonlinearExpr only occur in
            # other GenericNonlinearExpr.
            MOI.Nonlinear.parse_expression(data, expr, arg, parent_node)
        end
    end
    return
end

function _get_node_type(data, x)
    id = get(data.operators.univariate_operator_to_id, x.head, nothing)
    if length(x) == 1 && id !== nothing
        return id, MOI.Nonlinear.NODE_CALL_UNIVARIATE
    end
    id = get(data.operators.multivariate_operator_to_id, x.head, nothing)
    if id !== nothing
        return id, MOI.Nonlinear.NODE_CALL_MULTIVARIATE
    end
    id = get(data.operators.comparison_operator_to_id, x.head, nothing)
    if id !== nothing
        return id, MOI.Nonlinear.NODE_COMPARISON
    end
    id = get(data.operators.logic_operator_to_id, x.head, nothing)
    if id !== nothing
        return id, MOI.Nonlinear.NODE_LOGIC
    end
    return throw(MOI.UnsupportedNonlinearOperator(x.head))
end

function _parse_without_recursion_inner(stack, data, expr, x, parent)
    id, node_type = _get_node_type(data, x)
    push!(expr.nodes, MOI.Nonlinear.Node(node_type, id, parent))
    parent = length(expr.nodes)
    for i in length(x):-1:1  # Args need to be pushed onto the stack in reverse
        push!(stack, (parent, x[i]))
    end
    return
end

# Method definitions

function Base.zero(::Type{GenericNonlinearExpr{V}}) where {V}
    return GenericNonlinearExpr{V}(:+, 0.0)
end

function Base.one(::Type{GenericNonlinearExpr{V}}) where {V}
    return GenericNonlinearExpr{V}(:+, 1.0)
end

# Univariate operators

for f in MOI.Nonlinear.DEFAULT_UNIVARIATE_OPERATORS
    op = Meta.quot(f)
    if f == :+
        continue  # We don't need this.
    elseif f == :-
        @eval function Base.:-(x::GenericNonlinearExpr{V}) where {V}
            return GenericNonlinearExpr{V}(:-, x)
        end
    elseif isdefined(Base, f)
        @eval function Base.$(f)(x::AbstractJuMPScalar)
            return GenericNonlinearExpr{variable_ref_type(x)}($op, x)
        end
    elseif isdefined(MOI.Nonlinear, :SpecialFunctions)
        # The operator is defined in some other package.
        SF = MOI.Nonlinear.SpecialFunctions
        if isdefined(SF, f)
            @eval function $(SF).$(f)(x::AbstractJuMPScalar)
                return GenericNonlinearExpr{variable_ref_type(x)}($op, x)
            end
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
            return GenericNonlinearExpr{variable_ref_type(x)}($op, x, rhs)
        end
        function Base.$(f)(x::_Constant, y::AbstractJuMPScalar)
            lhs = convert(Float64, _constant_to_number(x))
            return GenericNonlinearExpr{variable_ref_type(y)}($op, lhs, y)
        end
        function Base.$(f)(x::AbstractJuMPScalar, y::AbstractJuMPScalar)
            return GenericNonlinearExpr{variable_ref_type(x)}($op, x, y)
        end
    end
end

function _MA.operate!!(
    ::typeof(_MA.add_mul),
    x::GenericNonlinearExpr,
    y::AbstractJuMPScalar,
)
    if x.head == :+
        push!(x.args, y)
        return x
    end
    return +(x, y)
end

"""
    flatten(expr::GenericNonlinearExpr)

Flatten a nonlinear expression by lifting nested `+` and `*` nodes into a single
n-ary operation.

## Motivation

Nonlinear expressions created using operator overloading can be deeply nested
and unbalanced. For example, `prod(x for i in 1:4)` creates
`*(x, *(x, *(x, x)))` instead of the more preferable `*(x, x, x, x)`.

## Example

```jldoctest
julia> model = Model();

julia> @variable(model, x)
x

julia> y = prod(x for i in 1:4)
((x²) * x) * x

julia> flatten(y)
(x²) * x * x

julia> flatten(sin(y))
sin((x²) * x * x)
```
"""
function flatten(expr::GenericNonlinearExpr{V}) where {V}
    root = GenericNonlinearExpr{V}(expr.head, Any[])
    nodes_to_visit = Any[(root, arg) for arg in reverse(expr.args)]
    while !isempty(nodes_to_visit)
        parent, arg = pop!(nodes_to_visit)
        if !(arg isa GenericNonlinearExpr)
            # Not a nonlinear expression, so can use recursion.
            push!(parent.args, flatten(arg))
        elseif parent.head in (:+, :*) && arg.head == parent.head
            # A special case: the arg can be lifted to an n-ary argument of the
            # parent.
            for n in reverse(arg.args)
                push!(nodes_to_visit, (parent, n))
            end
        else
            # The default case for nonlinear expressions. Put the args on the
            # stack, so that we may walk them later.
            for n in reverse(arg.args)
                push!(nodes_to_visit, (arg, n))
            end
            empty!(arg.args)
            push!(parent.args, arg)
        end
    end
    return root
end

flatten(expr) = expr

function _ifelse(a::AbstractJuMPScalar, x, y)
    return GenericNonlinearExpr{variable_ref_type(a)}(:ifelse, Any[a, x, y])
end

for (f, op) in (
    :_and => :&&,
    :_or => :||,
    :_less_than => :(<),
    :_greater_than => :(>),
    :_less_equal => :(<=),
    :_greater_equal => :(>=),
    :_equal_to => :(==),
)
    op = Meta.quot(op)
    @eval begin
        function $(f)(x::AbstractJuMPScalar, y)
            return GenericNonlinearExpr{variable_ref_type(x)}($op, x, y)
        end
        function $(f)(x, y::AbstractJuMPScalar)
            return GenericNonlinearExpr{variable_ref_type(y)}($op, x, y)
        end
        function $(f)(x::AbstractJuMPScalar, y::AbstractJuMPScalar)
            return GenericNonlinearExpr{variable_ref_type(x)}($op, x, y)
        end
    end
end

# JuMP interop

function owner_model(expr::GenericNonlinearExpr)
    for arg in expr.args
        if !(arg isa AbstractJuMPScalar)
            continue
        end
        model = owner_model(arg)
        if model !== nothing
            return model
        end
    end
    return nothing
end

function check_belongs_to_model(
    expr::GenericNonlinearExpr,
    model::AbstractModel,
)
    for arg in expr.args
        if arg isa AbstractJuMPScalar
            check_belongs_to_model(arg, model)
        end
    end
    return
end

function moi_function(f::GenericNonlinearExpr)
    ret = MOI.ScalarNonlinearFunction(f.head, Any[])
    stack = Tuple{MOI.ScalarNonlinearFunction,Any}[]
    for arg in reverse(f.args)
        push!(stack, (ret, arg))
    end
    while !isempty(stack)
        parent, arg = pop!(stack)
        if arg isa GenericNonlinearExpr
            new_ret = MOI.ScalarNonlinearFunction(arg.head, Any[])
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

function jump_function(model::GenericModel, f::MOI.ScalarNonlinearFunction)
    V = variable_ref_type(typeof(model))
    ret = GenericNonlinearExpr{V}(f.head, Any[])
    stack = Tuple{GenericNonlinearExpr,Any}[]
    for arg in reverse(f.args)
        push!(stack, (ret, arg))
    end
    while !isempty(stack)
        parent, arg = pop!(stack)
        if arg isa MOI.ScalarNonlinearFunction
            new_ret = GenericNonlinearExpr{V}(arg.head, Any[])
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

function jump_function_type(
    model::GenericModel,
    ::Type{<:MOI.ScalarNonlinearFunction},
)
    return GenericNonlinearExpr{variable_ref_type(typeof(model))}
end

moi_function_type(::Type{<:GenericNonlinearExpr}) = MOI.ScalarNonlinearFunction

function constraint_object(c::NonlinearConstraintRef)
    nlp = nonlinear_model(c.model)
    data = nlp.constraints[index(c)]
    return ScalarConstraint(jump_function(c.model, data.expression), data.set)
end

function jump_function(model::GenericModel, expr::MOI.Nonlinear.Expression)
    V = variable_ref_type(typeof(model))
    nlp = nonlinear_model(model)
    parsed = Vector{Any}(undef, length(expr.nodes))
    adj = MOI.Nonlinear.adjacency_matrix(expr.nodes)
    rowvals = SparseArrays.rowvals(adj)
    for i in length(expr.nodes):-1:1
        node = expr.nodes[i]
        parsed[i] = if node.type == MOI.Nonlinear.NODE_CALL_UNIVARIATE
            GenericNonlinearExpr{V}(
                nlp.operators.univariate_operators[node.index],
                parsed[rowvals[SparseArrays.nzrange(adj, i)[1]]],
            )
        elseif node.type == MOI.Nonlinear.NODE_CALL_MULTIVARIATE
            GenericNonlinearExpr{V}(
                nlp.operators.multivariate_operators[node.index],
                Any[parsed[rowvals[j]] for j in SparseArrays.nzrange(adj, i)],
            )
        elseif node.type == MOI.Nonlinear.NODE_MOI_VARIABLE
            V(model, MOI.VariableIndex(node.index))
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

function value(f::Function, expr::GenericNonlinearExpr)
    return _evaluate_expr(MOI.Nonlinear.OperatorRegistry(), f, expr)
end

function value(a::GenericNonlinearExpr; result::Int = 1)
    return value(a) do x
        return value(x; result = result)
    end
end

function _evaluate_expr(
    ::MOI.Nonlinear.OperatorRegistry,
    f::Function,
    expr::AbstractJuMPScalar,
)
    return value(f, expr)
end

function _evaluate_expr(
    ::MOI.Nonlinear.OperatorRegistry,
    ::Function,
    expr::Real,
)
    return convert(Float64, expr)
end

function _evaluate_user_defined_function(
    registry,
    f,
    expr::GenericNonlinearExpr,
)
    model = owner_model(expr)
    op, nargs = expr.head, length(expr.args)
    udf = MOI.get(model, MOI.UserDefinedFunction(op, nargs))
    if udf === nothing
        return error(
            "Unable to evaluate nonlinear operator $op because it is not " *
            "registered",
        )
    end
    args = [_evaluate_expr(registry, f, arg) for arg in expr.args]
    return first(udf)(args...)
end

function _evaluate_expr(
    registry::MOI.Nonlinear.OperatorRegistry,
    f::Function,
    expr::GenericNonlinearExpr,
)
    op = expr.head
    # TODO(odow): uses private function
    if !MOI.Nonlinear._is_registered(registry, op, length(expr.args))
        return _evaluate_user_defined_function(registry, f, expr)
    end
    if length(expr.args) == 1 && haskey(registry.univariate_operator_to_id, op)
        arg = _evaluate_expr(registry, f, expr.args[1])
        return MOI.Nonlinear.eval_univariate_function(registry, op, arg)
    elseif haskey(registry.multivariate_operator_to_id, op)
        args = [_evaluate_expr(registry, f, arg) for arg in expr.args]
        return MOI.Nonlinear.eval_multivariate_function(registry, op, args)
    elseif haskey(registry.logic_operator_to_id, op)
        @assert length(expr.args) == 2
        x = _evaluate_expr(registry, f, expr.args[1])
        y = _evaluate_expr(registry, f, expr.args[2])
        return MOI.Nonlinear.eval_logic_function(registry, op, x, y)
    else
        @assert haskey(registry.comparison_operator_to_id, op)
        @assert length(expr.args) == 2
        x = _evaluate_expr(registry, f, expr.args[1])
        y = _evaluate_expr(registry, f, expr.args[2])
        return MOI.Nonlinear.eval_comparison_function(registry, op, x, y)
    end
end

# MutableArithmetics.jl

# These converts are used in the {add,sub}mul definition for AbstractJuMPScalar.

Base.convert(::Type{<:GenericNonlinearExpr}, x::AbstractVariableRef) = x

function Base.convert(
    ::Type{<:GenericNonlinearExpr},
    x::GenericAffExpr{C,V},
) where {C,V}
    args = Any[]
    for (variable, coef) in x.terms
        if isone(coef)
            push!(args, variable)
        elseif !iszero(coef)
            push!(args, GenericNonlinearExpr{V}(:*, coef, variable))
        end
    end
    if !iszero(x.constant) || isempty(args)
        push!(args, x.constant)
    end
    if length(args) == 1
        return args[1]
    end
    return GenericNonlinearExpr{V}(:+, args)
end

function Base.convert(
    ::Type{<:GenericNonlinearExpr},
    x::GenericQuadExpr{C,V},
) where {C,V}
    args = Any[]
    for (variable, coef) in x.aff.terms
        if isone(coef)
            push!(args, variable)
        elseif !iszero(coef)
            push!(args, GenericNonlinearExpr{V}(:*, coef, variable))
        end
    end
    for (pair, coef) in x.terms
        if isone(coef)
            push!(args, GenericNonlinearExpr{V}(:*, pair.a, pair.b))
        elseif !iszero(coef)
            push!(args, GenericNonlinearExpr{V}(:*, coef, pair.a, pair.b))
        end
    end
    if !iszero(x.aff.constant) || isempty(args)
        push!(args, x.aff.constant)
    end
    if length(args) == 1
        return args[1]
    end
    return GenericNonlinearExpr{V}(:+, args)
end

function _MA.promote_operation(
    ::Union{typeof(+),typeof(-),typeof(*)},
    ::Type{GenericNonlinearExpr{V}},
    ::Type{<:AbstractJuMPScalar},
) where {V<:AbstractVariableRef}
    return GenericNonlinearExpr{V}
end

function _MA.promote_operation(
    ::Union{typeof(+),typeof(-),typeof(*)},
    ::Type{<:AbstractJuMPScalar},
    ::Type{GenericNonlinearExpr{V}},
) where {V<:AbstractVariableRef}
    return GenericNonlinearExpr{V}
end

function _MA.promote_operation(
    ::Union{typeof(+),typeof(-),typeof(*)},
    ::Type{GenericNonlinearExpr{V}},
    ::Type{GenericNonlinearExpr{V}},
) where {V<:AbstractVariableRef}
    return GenericNonlinearExpr{V}
end

function _MA.promote_operation(
    ::Union{typeof(+),typeof(-),typeof(*)},
    ::Type{GenericNonlinearExpr{U}},
    ::Type{GenericNonlinearExpr{V}},
) where {U<:AbstractVariableRef,V<:AbstractVariableRef}
    return error(
        "Unable to promote two different types of nonlinear expression",
    )
end

"""
    UserDefinedFunction(head::Symbol)

A struct representing a user-defined function named `head`. This function must
have already been added to the model using [`add_user_defined_function`](@ref)
or [`@register`](@ref).

## Example

```jldoctest
julia> model = Model();

julia> @variable(model, x)
x

julia> f(x::Float64) = x^2
f (generic function with 1 method)

julia> ∇f(x::Float64) = 2 * x
∇f (generic function with 1 method)

julia> ∇²f(x::Float64) = 2.0
∇²f (generic function with 1 method)

julia> add_user_defined_function(model, :foo, 1, f, ∇f, ∇²f)
UserDefinedFunction(:foo)

julia> bar = UserDefinedFunction(:foo)
UserDefinedFunction(:foo)

julia> @objective(model, Min, bar(x))
foo(x)
```
"""
struct UserDefinedFunction
    head::Symbol
end

function (f::UserDefinedFunction)(args...)
    return GenericNonlinearExpr(f.head, Any[a for a in args])
end

"""
    add_user_defined_function(
        model::Model,
        op::Symbol,
        dim::Int,
        f::Function,
        [∇f::Function,]
        [∇²f::Function,]
    )

Add a user-definend function with `dim` input arguments to `model` and associate
it with the operator `op`.

The function `f` evaluates the function. The optional function `∇f` evaluates
the first derivative, and the optional function `∇²f` evaluates the second
derivative. `∇²f` may be provided only if `∇f` is also provided.

## Example

```jldoctest
julia> model = Model();

julia> @variable(model, x)
x

julia> f(x::Float64) = x^2
f (generic function with 1 method)

julia> ∇f(x::Float64) = 2 * x
∇f (generic function with 1 method)

julia> ∇²f(x::Float64) = 2.0
∇²f (generic function with 1 method)

julia> foo = add_user_defined_function(model, :foo, 1, f, ∇f, ∇²f)
UserDefinedFunction(:foo)

julia> @objective(model, Min, foo(x))
foo(x)
```
"""
function add_user_defined_function(
    model::GenericModel,
    op::Symbol,
    dim::Int,
    args::Vararg{Function,N},
) where {N}
    if !(1 <= N <= 3)
        error(
            "Unable to register user-defined function $op: invalid number of " *
            "functions provided. Got $N, but expected 1 (if function only), " *
            "2 (if function and gradient), or 3 (if function, gradient, and " *
            "hesssian provided)",
        )
    end
    # TODO(odow): we could add other checks here, but we won't for now because
    # down-stream solvers in MOI can add their own checks, and any solver using
    # MOI.Nonlinear will automatically check for autodiff and common mistakes
    # and throw a nice informative error.
    MOI.set(model, MOI.UserDefinedFunction(op, dim), args)
    return UserDefinedFunction(op)
end

"""
    @register(model, operator, dim, args...)

Register a user-defined function in `model`, and create a new variable
[`UserDefinedFunction`](@ref) called `operator` in the current scope.

## Example

```jldoctest
julia> model = Model();

julia> @variable(model, x)
x

julia> f(x::Float64) = x^2
f (generic function with 1 method)

julia> ∇f(x::Float64) = 2 * x
∇f (generic function with 1 method)

julia> ∇²f(x::Float64) = 2.0
∇²f (generic function with 1 method)

julia> @register(model, foo, 1, f, ∇f, ∇²f)
UserDefinedFunction(:foo)

julia> @objective(model, Min, foo(x))
foo(x)
```

## Non-macro version

This macro is provided as helpful syntax that matches the style of the rest of
the JuMP macros. However, you may also create user-defined functions outside the
macros using [`add_user_defined_function`](@ref). For example:

```julia
julia> model = Model();

julia> @register(model, f, 1, x -> x^2)
UserDefinedFunction(:f)
```
is equivalent to
```julia
julia> model = Model();

julia> f = add_user_defined_function(model, :f, 1, x -> x^2)
UserDefinedFunction(:f)
```
"""
macro register(model, op, args...)
    rhs = Expr(
        :call,
        add_user_defined_function,
        esc(model),
        Meta.quot(op),
        esc.(args)...,
    )
    return Expr(:(=), esc(op), rhs)
end

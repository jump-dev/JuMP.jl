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

Additional operators can be add using [`@operator`](@ref).

See the full list of operators supported by a [`MOI.ModelLike`](@ref) by
querying the [`MOI.ListOfSupportedNonlinearOperators`](@ref) attribute.

## `args`

The vector `args` contains the arguments to the nonlinear function. If the
operator is univariate, it must contain one element. Otherwise, it may contain
multiple elements.

Given a subtype of [`AbstractVariableRef`](@ref), `V`, for `GenericNonlinearExpr{V}`,
each element must be one of the following:

 * A constant value of type `<:Real`
 * A `V`
 * A [`GenericAffExpr{T,V}`](@ref)
 * A [`GenericQuadExpr{T,V}`](@ref)
 * A [`GenericNonlinearExpr{V}`](@ref)

where `T<:Real` and `T == value_type(V)`.

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

    function GenericNonlinearExpr{V}(
        head::Symbol,
        args::Vararg{Any},
    ) where {V<:AbstractVariableRef}
        for arg in args
            _throw_if_not_real(arg)
            _throw_if_legacy_error(arg)
        end
        return new{V}(head, Any[a for a in args])
    end

    function GenericNonlinearExpr{V}(
        head::Symbol,
        args::Vector{Any},
    ) where {V<:AbstractVariableRef}
        for arg in args
            _throw_if_not_real(arg)
            _throw_if_legacy_error(arg)
        end
        return new{V}(head, args)
    end
end

variable_ref_type(::Type{GenericNonlinearExpr}, ::Any) = nothing

function variable_ref_type(::Type{GenericNonlinearExpr}, x::AbstractJuMPScalar)
    return variable_ref_type(x)
end

function _has_variable_ref_type(a)
    return variable_ref_type(GenericNonlinearExpr, a) !== nothing
end

function _variable_ref_type(head, args)
    if (i = findfirst(_has_variable_ref_type, args)) !== nothing
        V = variable_ref_type(GenericNonlinearExpr, args[i])
        return V::Type{<:AbstractVariableRef}
    end
    return error(
        "Unable to create a nonlinear expression because it did not contain " *
        "any JuMP scalars. head = `:$head`, args = `$args`.",
    )
end

function GenericNonlinearExpr(head::Symbol, args::Vector{Any})
    return GenericNonlinearExpr{_variable_ref_type(head, args)}(head, args)
end

function GenericNonlinearExpr(head::Symbol, args::Vararg{Any,N}) where {N}
    return GenericNonlinearExpr{_variable_ref_type(head, args)}(head, args...)
end

"""
    NonlinearExpr

Alias for `GenericNonlinearExpr{VariableRef}`, the specific
[`GenericNonlinearExpr`](@ref) used by JuMP.
"""
const NonlinearExpr = GenericNonlinearExpr{VariableRef}

variable_ref_type(::Type{GenericNonlinearExpr{V}}) where {V} = V

const _PREFIX_OPERATORS =
    (:+, :-, :*, :/, :^, :||, :&&, :>, :<, :(<=), :(>=), :(==))

_needs_parentheses(::Union{Number,AbstractVariableRef}) = false
_needs_parentheses(::Any) = true
function _needs_parentheses(x::GenericNonlinearExpr)
    return x.head in _PREFIX_OPERATORS && length(x.args) > 1
end

_parens(::MIME) = "(", ")", "", "", ""
_parens(::MIME"text/latex") = "\\left(", "\\right)", "{", "}", "\\textsf"

"""
    op_string(mime::MIME, x::GenericNonlinearExpr, ::Val{op}) where {op}

Return the string that should be printed for the operator `op` when
[`function_string`](@ref) is called with `mime` and `x`.

## Example

```jldoctest
julia> model = Model();

julia> @variable(model, x[1:2], Bin);

julia> f = @expression(model, x[1] || x[2]);

julia> op_string(MIME("text/plain"), f, Val(:||))
"||"
```
"""
op_string(::MIME, ::GenericNonlinearExpr, ::Val{op}) where {op} = string(op)
op_string(::MIME"text/latex", ::GenericNonlinearExpr, ::Val{:&&}) = "\\wedge"
op_string(::MIME"text/latex", ::GenericNonlinearExpr, ::Val{:||}) = "\\vee"
op_string(::MIME"text/latex", ::GenericNonlinearExpr, ::Val{:<=}) = "\\le"
op_string(::MIME"text/latex", ::GenericNonlinearExpr, ::Val{:>=}) = "\\ge"
op_string(::MIME"text/latex", ::GenericNonlinearExpr, ::Val{:(==)}) = "="

function function_string(mime::MIME, x::GenericNonlinearExpr)
    p_left, p_right, p_open, p_close, p_textsf = _parens(mime)
    io, stack = IOBuffer(), Any[x]
    while !isempty(stack)
        arg = pop!(stack)
        if arg isa GenericNonlinearExpr
            op = op_string(mime, arg, Val(arg.head))
            if arg.head in _PREFIX_OPERATORS && length(arg.args) > 1
                print(io, p_open)
                push!(stack, p_close)
                l = ceil(_TERM_LIMIT_FOR_PRINTING[] / 2)
                r = floor(_TERM_LIMIT_FOR_PRINTING[] / 2)
                skip_indices = (1+l):(length(arg.args)-r)
                for i in length(arg.args):-1:1
                    if i in skip_indices
                        if i == skip_indices[end]
                            push!(
                                stack,
                                _terms_omitted(mime, length(skip_indices)),
                            )
                            push!(stack, " $op $p_open")
                        end
                        continue
                    elseif _needs_parentheses(arg.args[i])
                        push!(stack, p_right)
                        push!(stack, arg.args[i])
                        push!(stack, p_left)
                    else
                        push!(stack, arg.args[i])
                    end
                    if i > 1
                        push!(stack, "$p_close $op $p_open")
                    end
                end
            else
                print(io, p_textsf, p_open, op, p_close, p_left, p_open)
                push!(stack, p_close * p_right)
                for i in length(arg.args):-1:2
                    push!(stack, arg.args[i])
                    push!(stack, "$p_close, $p_open")
                end
                if length(arg.args) >= 1
                    push!(stack, arg.args[1])
                end
            end
        elseif arg isa AbstractJuMPScalar
            print(io, function_string(mime, arg))
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
           length(x.args) == length(y.args) &&
           all(i -> _isequal(x.args[i], y.args[i]), 1:length(x.args))
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
            # We can use recursion here, because GenericNonlinearExpr only occur
            # in other GenericNonlinearExpr.
            MOI.Nonlinear.parse_expression(data, expr, arg, parent_node)
        end
    end
    return
end

function _get_node_type(data, x::GenericNonlinearExpr)
    id = get(data.operators.univariate_operator_to_id, x.head, nothing)
    if length(x.args) == 1 && id !== nothing
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

function _parse_without_recursion_inner(
    stack,
    data,
    expr,
    x::GenericNonlinearExpr,
    parent,
)
    id, node_type = _get_node_type(data, x)
    push!(expr.nodes, MOI.Nonlinear.Node(node_type, id, parent))
    parent = length(expr.nodes)
    # Args need to be pushed onto the stack in reverse
    for i in length(x.args):-1:1
        push!(stack, (parent, x.args[i]))
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

Base.conj(x::GenericNonlinearExpr) = x
Base.real(x::GenericNonlinearExpr) = x
Base.imag(x::GenericNonlinearExpr) = zero(x)
Base.abs2(x::GenericNonlinearExpr) = x^2
Base.isreal(::GenericNonlinearExpr) = true

# Univariate operators

_is_real(::Any) = false
_is_real(::Real) = true
_is_real(::AbstractVariableRef) = true
_is_real(::GenericAffExpr{<:Real}) = true
_is_real(::GenericQuadExpr{<:Real}) = true
_is_real(::GenericNonlinearExpr) = true
_is_real(::NonlinearExpression) = true
_is_real(::NonlinearParameter) = true

function _throw_if_not_real(x)
    if !_is_real(x)
        error(
            "Cannot build `GenericNonlinearExpr` because a term is " *
            "complex-valued: `($x)::$(typeof(x))`",
        )
    end
    return
end

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
            _throw_if_not_real(x)
            return GenericNonlinearExpr{variable_ref_type(x)}($op, x)
        end
    elseif isdefined(MOI.Nonlinear, :SpecialFunctions)
        # The operator is defined in some other package.
        SF = MOI.Nonlinear.SpecialFunctions
        if isdefined(SF, f)
            @eval function $(SF).$(f)(x::AbstractJuMPScalar)
                _throw_if_not_real(x)
                return GenericNonlinearExpr{variable_ref_type(x)}($op, x)
            end
        end
    end
end

# Multivariate operators

# The multivariate operators in MOI are +, -, *, ^, /, ifelse, atan, min, max
#
# However, ifelse is a builtin, so we can't add methods to it.

# We need only very generic fallbacks for these, because all other cases are
# caught with more specific methods.
for f in (:+, :-, :*, :^, :/, :atan, :min, :max)
    op = Meta.quot(f)
    @eval begin
        function Base.$(f)(x::AbstractJuMPScalar, y::_Constant)
            _throw_if_not_real(x)
            _throw_if_not_real(y)
            rhs = convert(Float64, _constant_to_number(y))
            return GenericNonlinearExpr{variable_ref_type(x)}($op, x, rhs)
        end
        function Base.$(f)(x::_Constant, y::AbstractJuMPScalar)
            _throw_if_not_real(x)
            _throw_if_not_real(y)
            lhs = convert(Float64, _constant_to_number(x))
            return GenericNonlinearExpr{variable_ref_type(y)}($op, lhs, y)
        end
        function Base.$(f)(x::AbstractJuMPScalar, y::AbstractJuMPScalar)
            _throw_if_not_real(x)
            _throw_if_not_real(y)
            return GenericNonlinearExpr{variable_ref_type(x)}($op, x, y)
        end
    end
end

# Base has unary methods `min(x::Real) = x` and `max(x::Real) = x`, so I guess
# we need to replicate them.
Base.min(x::AbstractJuMPScalar) = x
Base.max(x::AbstractJuMPScalar) = x

function _MA.operate!!(
    ::typeof(_MA.add_mul),
    x::GenericNonlinearExpr,
    args::Vararg{Any,N},
) where {N}
    _throw_if_not_real(x)
    if any(isequal(_MA.Zero()), args)
        return x
    end
    # It may seem like we should do this performance optimization, but it is NOT
    # safe. See JuMP#3825. The issue is that even though we're calling operate!!
    # `x` is not mutable, because it may, amoungst other things, be aliased in
    # one of the args
    # elseif x.head == :+
    #     push!(x.args, *(args...))
    #     return x
    # end
    return +(x, *(args...))
end

function _MA.operate!!(
    ::typeof(_MA.add_mul),
    ::GenericNonlinearExpr,
    x::AbstractArray,
)
    return _throw_operator_error(_MA.add_mul, x)
end

"""
    flatten!(expr::GenericNonlinearExpr)

Flatten a nonlinear expression in-place by lifting nested `+` and `*` nodes into
a single n-ary operation.

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

julia> flatten!(y)
(x²) * x * x

julia> flatten!(sin(prod(x for i in 1:4)))
sin((x²) * x * x)
```
"""
function flatten!(expr::GenericNonlinearExpr{V}) where {V}
    if !any(Base.Fix1(_needs_flatten, expr), expr.args)
        return expr
    end
    stack = Tuple{GenericNonlinearExpr{V},Int,GenericNonlinearExpr{V}}[]
    for i in length(expr.args):-1:1
        if _needs_flatten(expr, expr.args[i])
            push!(stack, (expr, i, expr.args[i]))
        end
    end
    while !isempty(stack)
        parent, i, arg = pop!(stack)
        if parent.head in (:+, :*) && arg.head == parent.head
            n = length(parent.args)
            resize!(parent.args, n + length(arg.args) - 1)
            for j in length(arg.args):-1:1
                parent_index = j == 1 ? i : n + j - 1
                if _needs_flatten(parent, arg.args[j])
                    push!(stack, (parent, parent_index, arg.args[j]))
                else
                    parent.args[parent_index] = arg.args[j]
                end
            end
        else
            parent.args[i] = arg
            for j in length(arg.args):-1:1
                if _needs_flatten(arg, arg.args[j])
                    push!(stack, (arg, j, arg.args[j]))
                end
            end
        end
    end
    return expr
end

flatten!(expr) = expr

_is_expr(::Any, ::Any) = false
_is_expr(x::GenericNonlinearExpr, op::Symbol) = x.head == op

_needs_flatten(::GenericNonlinearExpr, ::Any) = false

function _needs_flatten(parent::GenericNonlinearExpr, arg::GenericNonlinearExpr)
    if _is_expr(parent, :+)
        return _is_expr(arg, :+)
    elseif _is_expr(parent, :*)
        return _is_expr(arg, :*)
    else
        # Heuristic: we decide to flatten if `parent` is not a + or * operator,
        # but if one level down there are + or * nodes. This let's us flatten
        # sin(+(x, +(y, z)) => sin(+(x, y, z)), but not a more complicated
        # expression like log(sin(+(x, +(y, z))).
        #
        # If you have a benchmark that requires modifying this code, consider
        # instead adding `flatten!(::Any; force::Bool)` that would allow the
        # user to override this decision and flatten the entire tree.
        return any(Base.Fix2(_is_expr, :+), arg.args) ||
               any(Base.Fix2(_is_expr, :*), arg.args)
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

moi_function(x::Number) = x

function moi_function(f::GenericNonlinearExpr{V}) where {V}
    ret = MOI.ScalarNonlinearFunction(f.head, similar(f.args))
    stack = Tuple{MOI.ScalarNonlinearFunction,Int,GenericNonlinearExpr{V}}[]
    for i in length(f.args):-1:1
        if f.args[i] isa GenericNonlinearExpr{V}
            push!(stack, (ret, i, f.args[i]))
        else
            ret.args[i] = moi_function(f.args[i])
        end
    end
    while !isempty(stack)
        parent, i, arg = pop!(stack)
        child = MOI.ScalarNonlinearFunction(arg.head, similar(arg.args))
        parent.args[i] = child
        for j in length(arg.args):-1:1
            if arg.args[j] isa GenericNonlinearExpr{V}
                push!(stack, (child, j, arg.args[j]))
            else
                child.args[j] = moi_function(arg.args[j])
            end
        end
    end
    return ret
end

jump_function(::GenericModel{T}, x::Number) where {T} = convert(T, x)

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
    nlp = nonlinear_model(c.model)::MOI.Nonlinear.Model
    data = nlp.constraints[index(c)]
    return ScalarConstraint(jump_function(c.model, data.expression), data.set)
end

function jump_function(model::GenericModel, expr::MOI.Nonlinear.Expression)
    V = variable_ref_type(typeof(model))
    nlp = nonlinear_model(model)::MOI.Nonlinear.Model
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
            "Unable to evaluate nonlinear operator $op because it was " *
            "not added as an operator.",
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

# MutableArithmetics.jl and promotion

function Base.promote_rule(
    ::Type{GenericNonlinearExpr{V}},
    ::Type{V},
) where {V<:AbstractVariableRef}
    return GenericNonlinearExpr{V}
end

function Base.promote_rule(
    ::Type{GenericNonlinearExpr{V}},
    ::Type{<:Number},
) where {V<:AbstractVariableRef}
    return GenericNonlinearExpr{V}
end

function Base.promote_rule(
    ::Type{GenericNonlinearExpr{V}},
    ::Type{<:Union{GenericAffExpr{C,V},GenericQuadExpr{C,V}}},
) where {C,V<:AbstractVariableRef}
    return GenericNonlinearExpr{V}
end

function Base.convert(::Type{GenericNonlinearExpr{V}}, x::V) where {V}
    return GenericNonlinearExpr{V}(:+, Any[x])
end

function Base.convert(::Type{GenericNonlinearExpr{V}}, x::Number) where {V}
    return GenericNonlinearExpr{V}(:+, Any[x])
end

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
    if length(args) == 1 && args[1] isa GenericNonlinearExpr{V}
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
    if length(args) == 1 && args[1] isa GenericNonlinearExpr{V}
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
    NonlinearOperator(func::Function, head::Symbol)

A callable struct (functor) representing a function named `head`.

When called with [`AbstractJuMPScalar`](@ref)s, the struct returns a
[`GenericNonlinearExpr`](@ref).

When called with non-JuMP types, the struct returns the evaluation of
`func(args...)`.

Unless `head` is special-cased by the optimizer, the operator must have already
been added to the model using [`add_nonlinear_operator`](@ref) or
[`@operator`](@ref).

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

julia> @operator(model, op_f, 1, f, ∇f, ∇²f)
NonlinearOperator(f, :op_f)

julia> bar = NonlinearOperator(f, :op_f)
NonlinearOperator(f, :op_f)

julia> @objective(model, Min, bar(x))
op_f(x)

julia> bar(2.0)
4.0
```
"""
struct NonlinearOperator{F}
    func::F
    head::Symbol
end

# Make it so that we don't print the complicated type parameter
function Base.show(io::IO, f::NonlinearOperator)
    return print(io, "NonlinearOperator($(f.func), :$(f.head))")
end

function (f::NonlinearOperator)(args::Vararg{Any,N}) where {N}
    types = variable_ref_type.(GenericNonlinearExpr, args)
    if (i = findfirst(!isnothing, types)) !== nothing
        return GenericNonlinearExpr{types[i]}(f.head, args...)
    end
    return f.func(args...)
end

"""
    add_nonlinear_operator(
        model::Model,
        dim::Int,
        f::Function,
        [∇f::Function,]
        [∇²f::Function];
        [name::Symbol = Symbol(f),]
    )

Add a new nonlinear operator with `dim` input arguments to `model` and
associate it with the name `name`.

The function `f` evaluates the operator and must return a scalar.

The optional function `∇f` evaluates the first derivative, and the optional
function `∇²f` evaluates the second derivative.

`∇²f` may be provided only if `∇f` is also provided.

## Univariate syntax

If `dim == 1`, then the method signatures of each function must be:

 * `f(::T)::T where {T<:Real}`
 * `∇f(::T)::T where {T<:Real}`
 * `∇²f(::T)::T where {T<:Real}`

## Multivariate syntax

If `dim > 1`, then the method signatures of each function must be:

 * `f(x::T...)::T where {T<:Real}`
 * `∇f(g::AbstractVector{T}, x::T...)::Nothing where {T<:Real}`
 * `∇²f(H::AbstractMatrix{T}, x::T...)::Nothing where {T<:Real}`

Where the gradient vector `g` and Hessian matrix `H` are filled in-place. For
the Hessian, you must fill in the non-zero lower-triangular entries only.
Setting an off-diagonal upper-triangular element may error.

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

julia> op_f = add_nonlinear_operator(model, 1, f, ∇f, ∇²f)
NonlinearOperator(f, :f)

julia> @objective(model, Min, op_f(x))
f(x)

julia> op_f(2.0)
4.0
```
"""
function add_nonlinear_operator(
    model::GenericModel,
    dim::Int,
    f::Function,
    args::Vararg{Function,N};
    name::Symbol = Symbol(f),
) where {N}
    nargs = 1 + N
    if !(1 <= nargs <= 3)
        error(
            "Unable to add operator $name: invalid number of functions " *
            "provided. Got $nargs, but expected 1 (if function only), 2 (if " *
            "function and gradient), or 3 (if function, gradient, and " *
            "hesssian provided)",
        )
    end
    # TODO(odow): we could add other checks here, but we won't for now because
    # down-stream solvers in MOI can add their own checks, and any solver using
    # MOI.Nonlinear will automatically check for autodiff and common mistakes
    # and throw a nice informative error.
    MOI.set(model, MOI.UserDefinedFunction(name, dim), tuple(f, args...))
    return NonlinearOperator(f, name)
end

function _catch_redefinition_constant_error(op::Symbol, f::Function, args...)
    if op == Symbol(f)
        error("""
        Unable to add the nonlinear operator `:$op` with the same name as
        an existing function.

        For example, this code will error:
        ```julia
        model = Model()
        f(x) = x^2
        @operator(model, f, 1, f)
        ```
        because it is equivalent to:
        ```julia
        model = Model()
        f(x) = x^2
        f = add_nonlinear_operator(model, 1, f; name = :f)
        ```

        To fix, use a unique name, like `op_$op`:
        ```julia
        model = Model()
        f(x) = x^2
        @operator(model, op_f, 1, f)
        @expression(model, op_f(x))
        ```
        """)
    end
    return
end

"""
    @operator(model, operator, dim, f[, ∇f[, ∇²f]])

Add the nonlinear operator `operator` in `model` with `dim` arguments, and
create a new [`NonlinearOperator`](@ref) object called `operator` in the current
scope.

The function `f` evaluates the operator and must return a scalar.

The optional function `∇f` evaluates the first derivative, and the optional
function `∇²f` evaluates the second derivative.

`∇²f` may be provided only if `∇f` is also provided.

## Univariate syntax

If `dim == 1`, then the method signatures of each function must be:

 * `f(::T)::T where {T<:Real}`
 * `∇f(::T)::T where {T<:Real}`
 * `∇²f(::T)::T where {T<:Real}`

## Multivariate syntax

If `dim > 1`, then the method signatures of each function must be:

 * `f(x::T...)::T where {T<:Real}`
 * `∇f(g::AbstractVector{T}, x::T...)::Nothing where {T<:Real}`
 * `∇²f(H::AbstractMatrix{T}, x::T...)::Nothing where {T<:Real}`

Where the gradient vector `g` and Hessian matrix `H` are filled in-place. For
the Hessian, you must fill in the non-zero lower-triangular entries only.
Setting an off-diagonal upper-triangular element may error.

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

julia> @operator(model, op_f, 1, f, ∇f, ∇²f)
NonlinearOperator(f, :op_f)

julia> @objective(model, Min, op_f(x))
op_f(x)

julia> op_f(2.0)
4.0

julia> model[:op_f]
NonlinearOperator(f, :op_f)

julia> model[:op_f](x)
op_f(x)
```

## Non-macro version

This macro is provided as helpful syntax that matches the style of the rest of
the JuMP macros. However, you may also add operators outside the macro
using [`add_nonlinear_operator`](@ref). For example:

```jldoctest
julia> model = Model();

julia> f(x) = x^2
f (generic function with 1 method)

julia> @operator(model, op_f, 1, f)
NonlinearOperator(f, :op_f)
```
is equivalent to
```jldoctest
julia> model = Model();

julia> f(x) = x^2
f (generic function with 1 method)

julia> op_f = model[:op_f] = add_nonlinear_operator(model, 1, f; name = :op_f)
NonlinearOperator(f, :op_f)
```
"""
macro operator(model, op, dim, f, args...)
    return _finalize_macro(
        esc(model),
        quote
            _catch_redefinition_constant_error($(Meta.quot(op)), $(esc(f)))
            add_nonlinear_operator(
                $(esc(model)),
                $(esc(dim)),
                $(esc(f)),
                $(esc.(args)...);
                name = $(Meta.quot(op)),
            )
        end,
        __source__;
        register_name = op,
        wrap_let = true,
    )
end

function jump_function_type(
    ::GenericModel{T},
    ::Type{MOI.VectorNonlinearFunction},
) where {T}
    return Vector{GenericNonlinearExpr{GenericVariableRef{T}}}
end

function jump_function(
    model::GenericModel{T},
    f::MOI.VectorNonlinearFunction,
) where {T}
    return GenericNonlinearExpr{GenericVariableRef{T}}[
        jump_function(model, fi) for fi in MOI.Utilities.eachscalar(f)
    ]
end

function moi_function_type(::Type{<:AbstractVector{<:GenericNonlinearExpr}})
    return MOI.VectorNonlinearFunction
end

function moi_function(f::AbstractVector{<:GenericNonlinearExpr})
    return MOI.VectorNonlinearFunction(f)
end

function MOI.VectorNonlinearFunction(f::Vector{<:AbstractJuMPScalar})
    return MOI.VectorNonlinearFunction(map(moi_function, f))
end

# LinearAlgebra overloads to throw nicer error messages. These may be changed to
# return expressions in the future.

function LinearAlgebra.det(::AbstractMatrix{<:AbstractJuMPScalar})
    return throw(MOI.UnsupportedNonlinearOperator(:det))
end

function LinearAlgebra.logdet(::AbstractMatrix{<:AbstractJuMPScalar})
    return throw(MOI.UnsupportedNonlinearOperator(:logdet))
end

function LinearAlgebra.norm(::AbstractArray{<:AbstractJuMPScalar}, ::Real)
    return throw(MOI.UnsupportedNonlinearOperator(:norm))
end

function LinearAlgebra.nullspace(::AbstractVector{<:AbstractJuMPScalar})
    return throw(MOI.UnsupportedNonlinearOperator(:nullspace))
end

function LinearAlgebra.nullspace(::AbstractMatrix{<:AbstractJuMPScalar})
    return throw(MOI.UnsupportedNonlinearOperator(:nullspace))
end

function LinearAlgebra.qr(::AbstractMatrix{<:AbstractJuMPScalar})
    return throw(MOI.UnsupportedNonlinearOperator(:qr))
end

function add_to_expression!(f::GenericNonlinearExpr, args...)
    return error(
        """
        `add_to_expression!` is not supported for expressions of type
        `$(typeof(f))` because they cannot be modified in-place.
        Instead of `add_to_expression!(expr, args..)`, use one of the following:
        ```julia
        expr += *(args...)
        # or
        import MutableArithmetics as MA
        expr = MA.add_mul!!(expr, args...)
        ```
        """,
    )
end

#  Copyright 2017, Iain Dunning, Joey Huchette, Miles Lubin, and contributors
#  This Source Code Form is subject to the terms of the Mozilla Public
#  License, v. 2.0. If a copy of the MPL was not distributed with this
#  file, You can obtain one at https://mozilla.org/MPL/2.0/.

function _create_binary_switch(ids, exprs)
    if length(exprs) <= 3
        out = Expr(:if, Expr(:call, :(==), :id, ids[1]), exprs[1])
        if length(exprs) > 1
            push!(out.args, _create_binary_switch(ids[2:end], exprs[2:end]))
        end
        return out
    else
        mid = length(exprs) >>> 1
        return Expr(
            :if,
            Expr(:call, :(<=), :id, ids[mid]),
            _create_binary_switch(ids[1:mid], exprs[1:mid]),
            _create_binary_switch(ids[mid+1:end], exprs[mid+1:end]),
        )
    end
end

let exprs = map(SYMBOLIC_UNIVARIATE_EXPRESSIONS) do arg
        return :(return $(arg[1])(x), $(arg[2]))
    end
    expr = _create_binary_switch(1:length(exprs), exprs)
    @eval @inline function _eval_univariate(id, x::T) where {T}
        $(expr)
        return error("Invalid operator_id")
    end
end

let exprs = map(SYMBOLIC_UNIVARIATE_EXPRESSIONS) do arg
        if arg === :(nothing)  # f''(x) isn't defined
            :(error("Invalid operator_id"))
        else
            :(return $(arg[3]))
        end
    end
    expr = _create_binary_switch(1:length(exprs), exprs)
    @eval @inline function _eval_univariate_2nd_deriv(id, x::T) where {T}
        $(expr)
        return error("Invalid operator_id")
    end
end

struct UnivariateOperator{F,F′,F′′}
    f::F
    f′::F′
    f′′::F′′
end

struct MultivariateOperator{F,F′}
    N::Int
    f::F
    ∇f::F′
    function MultivariateOperator{N}(f::Function, ∇f::Function) where {N}
        return new{typeof(f),typeof(∇f)}(N, f, ∇f)
    end
end

"""
    DEFAULT_UNIVARIATE_OPERATORS

The list of univariate operators that are supported by default.
"""
const DEFAULT_UNIVARIATE_OPERATORS = first.(SYMBOLIC_UNIVARIATE_EXPRESSIONS)

"""
    DEFAULT_MULTIVARIATE_OPERATORS

The list of multivariate operators that are supported by default.
"""
const DEFAULT_MULTIVARIATE_OPERATORS = [:+, :-, :*, :^, :/, :ifelse]

struct OperatorRegistry
    # NODE_CALL_UNIVARIATE
    univariate_operators::Vector{Symbol}
    univariate_operator_to_id::Dict{Symbol,Int}
    univariate_user_operator_start::Int
    registered_univariate_operators::Vector{UnivariateOperator}
    # NODE_CALL_MULTIVARIATE
    multivariate_operators::Vector{Symbol}
    multivariate_operator_to_id::Dict{Symbol,Int}
    multivariate_user_operator_start::Int
    registered_multivariate_operators::Vector{MultivariateOperator}
    # NODE_LOGIC
    logic_operators::Vector{Symbol}
    logic_operator_to_id::Dict{Symbol,Int}
    # NODE_COMPARISON
    comparison_operators::Vector{Symbol}
    comparison_operator_to_id::Dict{Symbol,Int}
    function OperatorRegistry()
        univariate_operators = copy(DEFAULT_UNIVARIATE_OPERATORS)
        multivariate_operators = copy(DEFAULT_MULTIVARIATE_OPERATORS)
        logic_operators = [:&&, :||]
        comparison_operators = [:<=, :(==), :>=, :<, :>]
        return new(
            # NODE_CALL_UNIVARIATE
            univariate_operators,
            Dict{Symbol,Int}(
                op => i for (i, op) in enumerate(univariate_operators)
            ),
            length(univariate_operators),
            UnivariateOperator[],
            # NODE_CALL
            multivariate_operators,
            Dict{Symbol,Int}(
                op => i for (i, op) in enumerate(multivariate_operators)
            ),
            length(multivariate_operators),
            MultivariateOperator[],
            # NODE_LOGIC
            logic_operators,
            Dict{Symbol,Int}(op => i for (i, op) in enumerate(logic_operators)),
            # NODE_COMPARISON
            comparison_operators,
            Dict{Symbol,Int}(
                op => i for (i, op) in enumerate(comparison_operators)
            ),
        )
    end
end

const _FORWARD_DIFF_METHOD_ERROR_HELPER = raw"""
Common reasons for this include:

 * the function assumes `Float64` will be passed as input, it must work for any
   generic `Real` type.
 * the function allocates temporary storage using `zeros(3)` or similar. This
   defaults to `Float64`, so use `zeros(T, 3)` instead.

As an example, instead of:
```julia
function my_function(x::Float64...)
    y = zeros(length(x))
    for i in 1:length(x)
        y[i] = x[i]^2
    end
    return sum(y)
end
```
use:
```julia
function my_function(x::T...) where {T<:Real}
    y = zeros(T, length(x))
    for i in 1:length(x)
        y[i] = x[i]^2
    end
    return sum(y)
end
```

Review the stacktrace below for more information, but it can often be hard to
understand why and where your function is failing. You can also debug this
outside of JuMP as follows:
```julia
import ForwardDiff

# If the input dimension is 1
x = 1.0
my_function(a) = a^2
ForwardDiff.derivative(my_function, x)

# If the input dimension is more than 1
x = [1.0, 2.0]
my_function(a, b) = a^2 + b^2
ForwardDiff.gradient(x -> my_function(x...), x)
```
"""

_intercept_ForwardDiff_MethodError(err, ::Symbol) = rethrow(err)

function _intercept_ForwardDiff_MethodError(::MethodError, op::Symbol)
    return error(
        "JuMP's autodiff of the user-defined function $(op) failed with a " *
        "MethodError.\n\n$(_FORWARD_DIFF_METHOD_ERROR_HELPER)",
    )
end

function _checked_derivative(f::F, op::Symbol) where {F}
    return function (x)
        try
            return ForwardDiff.derivative(f, x)
        catch err
            _intercept_ForwardDiff_MethodError(err, op)
        end
    end
end

"""
    _validate_register_assumptions(
        f::Function,
        name::Symbol,
        dimension::Integer,
    )

A function that attempts to check if `f` is suitable for registration via
[`register`](@ref) and throws an informative error if it is not.

Because we don't know the domain of `f`, this function may encounter false
negatives. But it should catch the majority of cases in which users supply
non-differentiable functions that rely on `::Float64` assumptions.
"""
function _validate_register_assumptions(
    f::Function,
    name::Symbol,
    dimension::Integer,
)
    # Assumption 1: check that `f` can be called with `Float64` arguments.
    y = 0.0
    try
        if dimension == 1
            y = f(0.0)
        else
            y = f(zeros(dimension))
        end
    catch
        # We hit some other error, perhaps we called a function like log(0).
        # Ignore for now, and hope that a useful error is shown to the user
        # during the solve.
    end
    if !(y isa Real)
        error(
            "Expected return type of `Float64` from the user-defined " *
            "function :$(name), but got `$(typeof(y))`.",
        )
    end
    # Assumption 2: check that `f` can be differentiated using `ForwardDiff`.
    try
        if dimension == 1
            ForwardDiff.derivative(f, 0.0)
        else
            ForwardDiff.gradient(x -> f(x...), zeros(dimension))
        end
    catch err
        if err isa MethodError
            error(
                "Unable to register the function :$name because it does not " *
                "support differentiation via ForwardDiff.\n\n" *
                _FORWARD_DIFF_METHOD_ERROR_HELPER,
            )
        end
        # We hit some other error, perhaps we called a function like log(0).
        # Ignore for now, and hope that a useful error is shown to the user
        # during the solve.
    end
    return
end

function UnivariateOperator(op::Symbol, f::Function)
    _validate_register_assumptions(f, op, 1)
    f′ = _checked_derivative(f, op)
    f′′ = _checked_derivative(f′, op)
    return UnivariateOperator(f, f′, f′′)
end

function UnivariateOperator(op::Symbol, f::Function, f′::Function)
    _validate_register_assumptions(f′, op, 1)
    f′′ = _checked_derivative(f′, op)
    return UnivariateOperator(f, f′, f′′)
end

function UnivariateOperator(::Symbol, f::Function, f′::Function, f′′::Function)
    return UnivariateOperator(f, f′, f′′)
end

function MultivariateOperator{N}(op::Symbol, f::Function) where {N}
    _validate_register_assumptions(f, op, N)
    g = x -> f(x...)
    ∇f = function (ret, x)
        try
            ForwardDiff.gradient!(ret, g, x)
        catch err
            _intercept_ForwardDiff_MethodError(err, op)
        end
        return
    end
    return MultivariateOperator{N}(g, ∇f)
end

function MultivariateOperator{N}(::Symbol, f::Function, ∇f::Function) where {N}
    return MultivariateOperator{N}(x -> f(x...), (g, x) -> ∇f(g, x...))
end

function register_operator(
    registry::OperatorRegistry,
    op::Symbol,
    nargs::Int,
    f::Function...,
)
    if nargs == 1
        if haskey(registry.univariate_operator_to_id, op)
            error("Operator $op is already registered.")
        elseif haskey(registry.multivariate_operator_to_id, op)
            error("Operator $op is already registered.")
        end
        operator = UnivariateOperator(op, f...)
        push!(registry.univariate_operators, op)
        push!(registry.registered_univariate_operators, operator)
        registry.univariate_operator_to_id[op] =
            length(registry.univariate_operators)
    else
        if haskey(registry.multivariate_operator_to_id, op)
            error("Operator $op is already registered.")
        elseif haskey(registry.univariate_operator_to_id, op)
            error("Operator $op is already registered.")
        end
        operator = MultivariateOperator{nargs}(op, f...)
        push!(registry.multivariate_operators, op)
        push!(registry.registered_multivariate_operators, operator)
        registry.multivariate_operator_to_id[op] =
            length(registry.multivariate_operators)
    end
    return
end

function _is_registered(registry::OperatorRegistry, op::Symbol, nargs::Int)
    if op in (:<=, :>=, :(==), :<, :>, :&&, :||)
        return true
    end
    if nargs == 1 && haskey(registry.univariate_operator_to_id, op)
        return true
    end
    return haskey(registry.multivariate_operator_to_id, op)
end

function _warn_auto_register(op::Symbol, nargs::Int)
    @warn("""Function $op automatically registered with $nargs arguments.

    Calling the function with a different number of arguments will result in an
    error.

    While you can safely ignore this warning, we recommend that you manually
    register the function as follows:
    ```Julia
    model = Model()
    register(model, :$op, $nargs, $op; autodiff = true)
    ```""")
    return
end

function register_if_needed(
    registry::OperatorRegistry,
    op::Symbol,
    nargs::Int,
    f::Function,
)
    if !_is_registered(registry, op, nargs)
        register_operator(registry, op, nargs, f)
        _warn_auto_register(op, nargs)
    end
    return
end

function assert_registered(registry::OperatorRegistry, op::Symbol, nargs::Int)
    if !_is_registered(registry, op, nargs)
        msg = """
        Unrecognized function \"$(op)\" used in nonlinear expression.

        You must register it as a user-defined function before building
        the model. For example, replacing `N` with the appropriate number
        of arguments, do:
        ```julia
        model = Model()
        register(model, :$(op), N, $(op), autodiff=true)
        # ... variables and constraints ...
        ```
        """
        error(msg)
    end
    return
end

check_return_type(::Type{T}, ret::T) where {T} = nothing

function check_return_type(::Type{T}, ret) where {T}
    return error(
        "Expected return type of $T from a user-defined function, but got " *
        "$(typeof(ret)).",
    )
end

function eval_univariate_function(
    registry::OperatorRegistry,
    op::Symbol,
    x::T,
) where {T}
    id = registry.univariate_operator_to_id[op]
    if id <= registry.univariate_user_operator_start
        f, _ = _eval_univariate(id, x)
        return f::T
    end
    offset = id - registry.univariate_user_operator_start
    operator = registry.registered_univariate_operators[offset]
    ret = operator.f(x)
    check_return_type(T, ret)
    return ret::T
end

function eval_univariate_gradient(
    registry::OperatorRegistry,
    op::Symbol,
    x::T,
) where {T}
    id = registry.univariate_operator_to_id[op]
    if id <= registry.univariate_user_operator_start
        _, f′ = _eval_univariate(id, x)
        return f′::T
    end
    offset = id - registry.univariate_user_operator_start
    operator = registry.registered_univariate_operators[offset]
    ret = operator.f′(x)
    check_return_type(T, ret)
    return ret::T
end

function eval_univariate_hessian(
    registry::OperatorRegistry,
    op::Symbol,
    x::T,
) where {T}
    id = registry.univariate_operator_to_id[op]
    if id <= registry.univariate_user_operator_start
        return _eval_univariate_2nd_deriv(id, x)
    end
    offset = id - registry.univariate_user_operator_start
    operator = registry.registered_univariate_operators[offset]
    ret = operator.f′′(x)
    check_return_type(T, ret)
    return ret::T
end

function eval_multivariate_function(
    registry::OperatorRegistry,
    op::Symbol,
    x::AbstractVector{T},
)::T where {T}
    if op == :+
        return +(x...)
    elseif op == :-
        return -(x...)
    elseif op == :*
        return *(x...)
    elseif op == :^
        @assert length(x) == 2
        return x[1]^x[2]
    elseif op == :/
        @assert length(x) == 2
        return x[1] / x[2]
    elseif op == :ifelse
        @assert length(x) == 3
        return ifelse(Bool(x[1]), x[2], x[3])
    end
    id = registry.multivariate_operator_to_id[op]
    offset = id - registry.multivariate_user_operator_start
    operator = registry.registered_multivariate_operators[offset]
    @assert length(x) == operator.N
    ret = operator.f(x)
    check_return_type(T, ret)
    return ret::T
end

function eval_multivariate_gradient(
    registry::OperatorRegistry,
    op::Symbol,
    g::AbstractVector{T},
    x::AbstractVector{T},
) where {T}
    if op == :+
        fill!(g, one(T))
    elseif op == :-
        fill!(g, -one(T))
        g[1] = one(T)
    elseif op == :*
        total = *(x...)
        if total == zero(T)
            for i in 1:length(x)
                g[i] = prod(x[j] for j in 1:length(x) if i != j)
            end
        else
            for i in 1:length(x)
                g[i] = total / x[i]
            end
        end
    elseif op == :^
        @assert length(x) == 2
        if x[2] == 1
            g[1] = one(T)
        elseif x[2] == 2
            g[1] = T(2) * x[1]
        else
            g[1] = x[2] * x[1]^(x[2] - one(T))
        end
        # TODO(odow): fix me. We should use NaNMath.jl here
        if x[1] < 0
            g[2] = T(NaN)
        else
            g[2] = x[1]^x[2] * log(x[1])
        end
    elseif op == :/
        @assert length(x) == 2
        g[1] = one(T) / x[2]
        g[2] = -x[1] / x[2]^2
    elseif op == :ifelse
        @assert length(x) == 3
        g[1] = zero(T)  # It doesn't matter what this is.
        g[2] = x[1] == one(T)
        g[3] = x[1] == zero(T)
    else
        id = registry.multivariate_operator_to_id[op]
        offset = id - registry.multivariate_user_operator_start
        operator = registry.registered_multivariate_operators[offset]
        @assert length(x) == operator.N
        operator.∇f(g, x)
    end
    return
end

function eval_multivariate_hessian(
    ::OperatorRegistry,
    ::Symbol,
    ::AbstractMatrix{T},
    ::AbstractVector{T},
) where {T}
    return error("Not implemented")
end

# These are not extendable!
function eval_logic_function(op::Symbol, lhs::T, rhs::T)::Bool where {T}
    if op == :&&
        return lhs && rhs
    else
        @assert op == :||
        return lhs || rhs
    end
end

# These are not extendable!
function eval_comparison_function(op::Symbol, lhs::T, rhs::T)::Bool where {T}
    if op == :<=
        return lhs <= rhs
    elseif op == :>=
        return lhs >= rhs
    elseif op == :(==)
        return lhs == rhs
    elseif op == :<
        return lhs < rhs
    else
        @assert op == :>
        return lhs > rhs
    end
end

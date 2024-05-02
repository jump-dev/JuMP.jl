#  Copyright 2017, Iain Dunning, Joey Huchette, Miles Lubin, and contributors
#  This Source Code Form is subject to the terms of the Mozilla Public
#  License, v. 2.0. If a copy of the MPL was not distributed with this
#  file, You can obtain one at https://mozilla.org/MPL/2.0/.

const _op_add = NonlinearOperator(+, :+)
const _op_sub = NonlinearOperator(-, :-)
const _op_mul = NonlinearOperator(*, :*)
const _op_div = NonlinearOperator(/, :/)
const _op_pow = NonlinearOperator(^, :^)

"""
    @force_nonlinear(expr)

Change the parsing of `expr` to construct [`GenericNonlinearExpr`](@ref) instead
of [`GenericAffExpr`](@ref) or [`GenericQuadExpr`](@ref).

This macro works by walking `expr` and substituting all calls to `+`, `-`, `*`,
`/`, and `^` in favor of ones that construct [`GenericNonlinearExpr`](@ref).

This macro will error if the resulting expression does not produce a
[`GenericNonlinearExpr`](@ref) because, for example, it is used on an expression
that does not use the basic arithmetic operators.

## When to use this macro

In most cases, you should not use this macro.

Use this macro only if the intended output type is a [`GenericNonlinearExpr`](@ref)
and the regular macro calls destroy problem structure, or in rare cases, if the
regular macro calls introduce a large amount of intermediate variables, for
example, because they promote types to a common quadratic expression.

## Example

### Use-case one: preserve problem structure.

```jldoctest
julia> model = Model();

julia> @variable(model, x);

julia> @expression(model, (x - 0.1)^2)
x² - 0.2 x + 0.010000000000000002

julia> @expression(model, @force_nonlinear((x - 0.1)^2))
(x - 0.1) ^ 2

julia> (x - 0.1)^2
x² - 0.2 x + 0.010000000000000002

julia> @force_nonlinear((x - 0.1)^2)
(x - 0.1) ^ 2
```

### Use-case two: reduce allocations

In this example, we know that `x * 2.0 * (1 + x) * x` is going to construct a
nonlinear expression.

However, the default parsing first constructs:

 * the [`GenericAffExpr`](@ref) `a = x * 2.0`,
 * another [`GenericAffExpr`](@ref) `b = 1 + x`
 * the [`GenericQuadExpr`](@ref) `c = a * b`
 * a [`GenericNonlinearExpr`](@ref) `*(c, x)`

In contrast, the modified parsing constructs:

 * the [`GenericNonlinearExpr`](@ref) `a = GenericNonlinearExpr(:+, 1, x)`
 * the [`GenericNonlinearExpr`](@ref) `GenericNonlinearExpr(:*, x, 2.0, a, x)`

This results in significantly fewer allocations.

```jldoctest
julia> model = Model();

julia> @variable(model, x);

julia> @expression(model, x * 2.0 * (1 + x) * x)
(2 x² + 2 x) * x

julia> @expression(model, @force_nonlinear(x * 2.0 * (1 + x) * x))
x * 2.0 * (1 + x) * x

julia> @allocated @expression(model, x * 2.0 * (1 + x) * x)
3200

julia> @allocated @expression(model, @force_nonlinear(x * 2.0 * (1 + x) * x))
640
```
"""
macro force_nonlinear(expr)
    error_fn = Containers.build_error_fn(:force_nonlinear, (expr,), __source__)
    ret = MacroTools.postwalk(expr) do x
        if Meta.isexpr(x, :call)
            if x.args[1] == :+
                return Expr(:call, _op_add, x.args[2:end]...)
            elseif x.args[1] == :-
                return Expr(:call, _op_sub, x.args[2:end]...)
            elseif x.args[1] == :*
                return Expr(:call, _op_mul, x.args[2:end]...)
            elseif x.args[1] == :/
                return Expr(:call, _op_div, x.args[2:end]...)
            elseif x.args[1] == :^
                return Expr(:call, _op_pow, x.args[2:end]...)
            end
        end
        return x
    end
    return Expr(:call, _force_nonlinear, error_fn, esc(ret))
end

_force_nonlinear(::F, ret::GenericNonlinearExpr) where {F} = ret

function _force_nonlinear(error_fn::F, ret::Any) where {F}
    return error_fn(
        "expression did not produce a `GenericNonlinearExpr`. Got a " *
        "`$(typeof(ret))`: $(ret)",
    )
end

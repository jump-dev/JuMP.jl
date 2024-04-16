#  Copyright 2017, Iain Dunning, Joey Huchette, Miles Lubin, and contributors
#  This Source Code Form is subject to the terms of the Mozilla Public
#  License, v. 2.0. If a copy of the MPL was not distributed with this
#  file, You can obtain one at https://mozilla.org/MPL/2.0/.

const _op_add = NonlinearOperator(+, :+)
const _op_sub = NonlinearOperator(-, :-)
const _op_mul = NonlinearOperator(*, :*)
const _op_div = NonlinearOperator(/, :/)

"""
    @nl(expr)

Change the parsing of `expr` to construct [`GenericNonlinearExpr`](@ref) instead
of [`GenericAffExpr`](@ref) or [`GenericQuadExpr`](@ref).

This macro works by walking `expr` and substituting all calls to `+`, `-`, `*`
and `/` in favor of ones that construct [`GenericNonlinearExpr`](@ref).

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

julia> @expression(model, @nl((x - 0.1)^2))
(x - 0.1) ^ 2.0

julia> (x - 0.1)^2
x² - 0.2 x + 0.010000000000000002

julia> @nl((x - 0.1)^2)
(x - 0.1) ^ 2.0
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

julia> @allocated @expression(model, x * 2.0 * (1 + x) * x)
3200

julia> @allocated @expression(model, @nl(x * 2.0 * (1 + x) * x))
640
```
"""
macro nl(expr)
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
            end
        end
        return x
    end
    return esc(ret)
end

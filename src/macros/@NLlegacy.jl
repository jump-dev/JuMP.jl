#  Copyright 2017, Iain Dunning, Joey Huchette, Miles Lubin, and contributors
#  This Source Code Form is subject to the terms of the Mozilla Public
#  License, v. 2.0. If a copy of the MPL was not distributed with this
#  file, You can obtain one at https://mozilla.org/MPL/2.0/.

"""
    _parse_nonlinear_expression(model::GenericModel, x::Expr)

JuMP needs to build Nonlinear expression objects in macro scope. This has two
main challenges:

 1. We need to evaluate local variables into the expressions. This is reasonably
    easy, anywhere we see a symbol that is not a function call, replace it by
    `esc(x)`.

 2. We need to identify un-registered user-defined functions so that we can
    attempt to automatically register them if their symbolic name exists in the
    scope. I (@odow) originally introduced the auto-registration in
    https://github.com/jump-dev/JuMP.jl/pull/2537 to fix a common pain-point in
    JuMP, but after working through this I believe it was a mistake. It's a lot
    of hassle! One problem is that the design of Nonlinear has moved the
    expression parsing from macro-expansion time to runtime. I think this is a
    big win for readability of the system, but it means we loose access to the
    caller's local scope. My solution to maintain backwards compatibility is to
    check that every function call is registered before parsing the expression.
"""
function _parse_nonlinear_expression(model, x)
    code = quote
        _init_NLP($model)
    end
    operators = Set{Tuple{Symbol,Int}}()
    _assert_constant_comparison(code, x)
    y = _parse_nonlinear_expression_inner(code, x, operators)
    user_defined_operators = filter(operators) do (op, i)
        if op in (:<=, :>=, :(==), :<, :>, :&&, :||)
            return false
        elseif i == 1 && op in MOI.Nonlinear.DEFAULT_UNIVARIATE_OPERATORS
            return false
        elseif i > 1 && op in MOI.Nonlinear.DEFAULT_MULTIVARIATE_OPERATORS
            return false
        end
        return true
    end
    if length(user_defined_operators) > 0
        op_var = gensym()
        push!(code.args, :($op_var = $(model).nlp_model.operators))
        for (op, i) in collect(user_defined_operators)
            push!(code.args, _auto_register_expression(op_var, op, i))
        end
    end
    return code, y
end

# JuMP special-cases the error for constant LHS and RHS terms in a comparison
# expression. In JuMP v1.1 and earlier, it evaluated the outer terms in the
# current scope, and checked to see if they were Real. To keep the same behavior
# we do the same here.
function _assert_constant_comparison(code::Expr, expr::Expr)
    if Meta.isexpr(expr, :comparison)
        lhs, rhs = gensym(), gensym()
        push!(code.args, esc(:($lhs = $(expr.args[1]))))
        push!(code.args, esc(:($rhs = $(expr.args[5]))))
        expr.args[1], expr.args[5] = lhs, rhs
    end
    return
end

_assert_constant_comparison(::Expr, ::Any) = nothing

function _auto_register_expression(op_var, op, i)
    q_op = Meta.quot(op)
    return quote
        try
            MOI.Nonlinear.register_operator_if_needed(
                $op_var,
                $q_op,
                $i,
                $(esc(op)),
            )
        catch
        end
        MOI.Nonlinear.assert_registered($op_var, $q_op, $i)
    end
end

function _normalize_unicode(x::Symbol)
    if x == :≤
        return :<=
    elseif x == :≥
        return :>=
    else
        return x
    end
end

function _parse_nonlinear_expression_inner(::Any, x::Symbol, ::Any)
    x = _normalize_unicode(x)
    if x in (:<=, :>=, :(==), :<, :>, :&&, :||)
        return Meta.quot(x)
    end
    return esc(x)
end

# Numbers and other literal constants.
_parse_nonlinear_expression_inner(::Any, x, ::Any) = x

function _is_generator(x)
    return Meta.isexpr(x, :call) &&
           length(x.args) >= 2 &&
           (
               Meta.isexpr(x.args[end], :generator) ||
               Meta.isexpr(x.args[end], :flatten)
           )
end

function _parse_nonlinear_expression_inner(code, x::Expr, operators)
    if Meta.isexpr(x, :block)
        error(
            "`begin...end` blocks are not supported in nonlinear macros. The " *
            "nonlinear expression must be a single statement.",
        )
    end
    if Meta.isexpr(x, :ref)
        return esc(x)
    elseif Meta.isexpr(x, :.)
        return esc(x)
    elseif _is_generator(x)
        return _parse_generator_expression(code, x, operators)
    elseif Meta.isexpr(x, Symbol("'"))
        # Treat the adjoint operator as a special case, because people often
        # use linear algebra in macros.
        return esc(x)
    end
    y = gensym()
    y_expr = :($y = Expr($(Meta.quot(x.head))))
    offset = 1
    if Meta.isexpr(x, :call)
        if !(x.args[1] isa Symbol)
            error(
                "Unsupported function $(x.args[1]). All function calls must " *
                "be `Symbol`s.",
            )
        end
        op = _normalize_unicode(x.args[1])
        push!(operators, (op, length(x.args) - 1))
        push!(y_expr.args[2].args, Meta.quot(op))
        offset += 1
    end
    for i in offset:length(x.args)
        arg = _parse_nonlinear_expression_inner(code, x.args[i], operators)
        push!(y_expr.args[2].args, arg)
    end
    push!(code.args, y_expr)
    return y
end

_is_sum(s::Symbol) = (s == :sum) || (s == :∑) || (s == :Σ)
_is_prod(s::Symbol) = (s == :prod) || (s == :∏)

function _parse_generator_expression(code, x, operators)
    y = gensym()
    y_expr, default = if _is_sum(x.args[1])
        :($y = Expr(:call, :+)), 0
    elseif _is_prod(x.args[1])
        :($y = Expr(:call, :*)), 1
    elseif x.args[1] == :maximum
        :($y = Expr(:call, :max)), nothing
    elseif x.args[1] == :minimum
        :($y = Expr(:call, :min)), nothing
    else
        error("Unsupported generator `:$(x.args[1])`")
    end
    body = x.args[2]
    has_init = false
    # foo(generator; init = value)
    if Meta.isexpr(x.args[2], :parameters)
        for kw in x.args[2].args
            if !Meta.isexpr(kw, :kw) || kw.args[1] != :init
                error("Unsupported nonlinear expression: $x")
            end
            if kw.args[2] != default
                default = esc(kw.args[2])
                push!(y_expr.args[2].args, default)
                has_init = true
            end
        end
        body = x.args[3]
    end
    # foo(generator, init = value)
    if Meta.isexpr(x.args[2], :generator, 3)
        kw = x.args[2].args[3]
        if Meta.isexpr(kw, :(=), 2) && kw.args[1] == :init
            pop!(x.args[2].args)
            if kw.args[2] != default
                default = esc(kw.args[2])
                push!(y_expr.args[2].args, default)
                has_init = true
            end
        end
    end
    block = _MA.rewrite_generator(
        body,
        t -> begin
            new_code = quote end
            arg = _parse_nonlinear_expression_inner(new_code, t, operators)
            push!(new_code.args, :(push!($y.args, $arg)))
            new_code
        end,
    )
    # Special case that was handled by JuMP in the past.
    error_string = "reducing over an empty collection in `$(x.args[1])` is not allowed"
    push!(code.args, quote
        $y_expr
        $block
        if length($y.args) == $(has_init ? 2 : 1)
            if $default === nothing
                throw(ArgumentError($error_string))
            else
                $y = $default
            end
        end
    end)
    return y
end

###
### @NLobjective(s)
###

"""
    @NLobjective(model, sense, expression)

Add a nonlinear objective to `model` with optimization sense `sense`.
`sense` must be `Max` or `Min`.

!!! compat
    This macro is part of the legacy nonlinear interface. Consider using the
    new nonlinear interface documented in [Nonlinear Modeling](@ref). In most
    cases, you can replace `@NLobjective` with [`@objective`](@ref).

## Example

```jldoctest
julia> model = Model();

julia> @variable(model, x)
x

julia> @NLobjective(model, Max, 2x + 1 + sin(x))

julia> print(model)
Max 2.0 * x + 1.0 + sin(x)
Subject to
```
"""
macro NLobjective(model, sense, x)
    error_fn =
        Containers.build_error_fn(:NLobjective, (model, sense, x), __source__)
    sense_expr = _parse_moi_sense(error_fn, sense)
    esc_model = esc(model)
    parsing_code, expr = _parse_nonlinear_expression(esc_model, x)
    code = quote
        $parsing_code
        set_nonlinear_objective($esc_model, $sense_expr, $expr)
    end
    return _finalize_macro(esc_model, code, __source__)
end

###
### @NLconstraint(s)
###

"""
    @NLconstraint(model::GenericModel, expr)

Add a constraint described by the nonlinear expression `expr`. See also
[`@constraint`](@ref).

!!! compat
    This macro is part of the legacy nonlinear interface. Consider using the
    new nonlinear interface documented in [Nonlinear Modeling](@ref). In most
    cases, you can replace `@NLconstraint` with [`@constraint`](@ref).

## Example

```jldoctest
julia> model = Model();

julia> @variable(model, x)
x

julia> @NLconstraint(model, sin(x) <= 1)
sin(x) - 1.0 ≤ 0

julia> @NLconstraint(model, [i = 1:3], sin(i * x) <= 1 / i)
3-element Vector{NonlinearConstraintRef{ScalarShape}}:
 (sin(1.0 * x) - 1.0 / 1.0) - 0.0 ≤ 0
 (sin(2.0 * x) - 1.0 / 2.0) - 0.0 ≤ 0
 (sin(3.0 * x) - 1.0 / 3.0) - 0.0 ≤ 0
```
"""
macro NLconstraint(m, x, args...)
    error_fn =
        Containers.build_error_fn(:NLconstraint, (m, x, args...), __source__)
    esc_m = esc(m)
    if Meta.isexpr(x, :block)
        error_fn("Invalid syntax. Did you mean to use `@NLconstraints`?")
    end
    # Two formats:
    # - @NLconstraint(m, a*x <= 5)
    # - @NLconstraint(m, myref[a=1:5], sin(x^a) <= 5)
    extra, kw_args, requested_container = Containers._extract_kw_args(args)
    if length(extra) > 1 || length(kw_args) > 0
        error_fn("too many arguments.")
    end
    # Canonicalize the arguments
    c = length(extra) == 1 ? x : nothing
    con = length(extra) == 1 ? extra[1] : x
    # Strategy: build up the code for non-macro add_constraint, and if needed
    # we will wrap in loops to assign to the ConstraintRefs
    name, idxvars, indices =
        Containers.parse_ref_sets(error_fn, c; invalid_index_variables = [m])
    parsing_code, expr = _parse_nonlinear_expression(esc_m, con)
    code = quote
        $parsing_code
        add_nonlinear_constraint($esc_m, $expr)
    end
    looped =
        Containers.container_code(idxvars, indices, code, requested_container)
    creation_code = quote
        _init_NLP($esc_m)
        $looped
    end
    return _finalize_macro(
        esc_m,
        creation_code,
        __source__;
        register_name = name,
    )
end

"""
    @NLconstraints(model, args...)

Adds multiple nonlinear constraints to model at once, in the same fashion as
the [`@NLconstraint`](@ref) macro.

The model must be the first argument, and multiple constraints can be added on
multiple lines wrapped in a `begin ... end` block.

The macro returns a tuple containing the constraints that were defined.

!!! compat
    This macro is part of the legacy nonlinear interface. Consider using the
    new nonlinear interface documented in [Nonlinear Modeling](@ref). In most
    cases, you can replace `@NLconstraints` with [`@constraints`](@ref).

## Example

```jldoctest
julia> model = Model();

julia> @variable(model, x);

julia> @variable(model, y);

julia> @variable(model, t);

julia> @variable(model, z[1:2]);

julia> a = [4, 5];

julia> @NLconstraints(model, begin
           t >= sqrt(x^2 + y^2)
           [i = 1:2], z[i] <= log(a[i])
       end)
((t - sqrt(x ^ 2.0 + y ^ 2.0)) - 0.0 ≥ 0, NonlinearConstraintRef{ScalarShape}[(z[1] - log(4.0)) - 0.0 ≤ 0, (z[2] - log(5.0)) - 0.0 ≤ 0])
```
"""
macro NLconstraints(model, block)
    return _plural_macro_code(model, block, Symbol("@NLconstraint"))
end

###
### @NLexpression(s)
###

"""
    @NLexpression(args...)

Efficiently build a nonlinear expression which can then be inserted in other
nonlinear constraints and the objective. See also [`@expression`].

!!! compat
    This macro is part of the legacy nonlinear interface. Consider using the
    new nonlinear interface documented in [Nonlinear Modeling](@ref). In most
    cases, you can replace `@NLexpression` with [`@expression`](@ref).

## Example

```jldoctest api_nlexpression
julia> model = Model();

julia> @variable(model, x)
x

julia> @variable(model, y)
y

julia> @NLexpression(model, my_expr, sin(x)^2 + cos(x^2))
subexpression[1]: sin(x) ^ 2.0 + cos(x ^ 2.0)

julia> @NLconstraint(model, my_expr + y >= 5)
(subexpression[1] + y) - 5.0 ≥ 0

julia> @NLobjective(model, Min, my_expr)
```

Indexing over sets and anonymous expressions are also supported:
```jldoctest api_nlexpression
julia> @NLexpression(model, my_expr_1[i=1:3], sin(i * x))
3-element Vector{NonlinearExpression}:
 subexpression[2]: sin(1.0 * x)
 subexpression[3]: sin(2.0 * x)
 subexpression[4]: sin(3.0 * x)

julia> my_expr_2 = @NLexpression(model, log(1 + sum(exp(my_expr_1[i]) for i in 1:2)))
subexpression[5]: log(1.0 + (exp(subexpression[2]) + exp(subexpression[3])))
```
"""
macro NLexpression(args...)
    error_fn = Containers.build_error_fn(:NLexpression, args, __source__)
    args, kw_args, requested_container = Containers._extract_kw_args(args)
    if length(args) <= 1
        error_fn(
            "To few arguments ($(length(args))); must pass the model and nonlinear expression as arguments.",
        )
    elseif length(args) == 2
        m, x = args
        c = nothing
    elseif length(args) == 3
        m, c, x = args
    end
    if Meta.isexpr(args[2], :block)
        error_fn("Invalid syntax. Did you mean to use `@NLexpressions`?")
    end
    if length(args) > 3 || length(kw_args) > 0
        error_fn("To many arguments ($(length(args))).")
    end
    name, idxvars, indices = Containers.parse_ref_sets(error_fn, c)
    name, idxvars, indices = Containers.parse_ref_sets(
        error_fn,
        c;
        invalid_index_variables = [args[1]],
    )
    esc_m = esc(m)
    parsing_code, expr = _parse_nonlinear_expression(esc_m, x)
    code = quote
        $parsing_code
        add_nonlinear_expression($esc_m, $expr)
    end
    creation_code =
        Containers.container_code(idxvars, indices, code, requested_container)
    return _finalize_macro(
        esc_m,
        creation_code,
        __source__;
        register_name = name,
    )
end

"""
    @NLexpressions(model, args...)

Adds multiple nonlinear expressions to model at once, in the same fashion as the
[`@NLexpression`](@ref) macro.

The model must be the first argument, and multiple expressions can be added on
multiple lines wrapped in a `begin ... end` block.

The macro returns a tuple containing the expressions that were defined.

!!! compat
    This macro is part of the legacy nonlinear interface. Consider using the
    new nonlinear interface documented in [Nonlinear Modeling](@ref). In most
    cases, you can replace `@NLexpressions` with [`@expressions`](@ref).

## Example

```jldoctest
julia> model = Model();

julia> @variable(model, x);

julia> @variable(model, y);

julia> @variable(model, z[1:2]);

julia> a = [4, 5];

julia> @NLexpressions(model, begin
           my_expr, sqrt(x^2 + y^2)
           my_expr_1[i = 1:2], log(a[i]) - z[i]
       end)
(subexpression[1]: sqrt(x ^ 2.0 + y ^ 2.0), NonlinearExpression[subexpression[2]: log(4.0) - z[1], subexpression[3]: log(5.0) - z[2]])
```
"""
macro NLexpressions(model, block)
    return _plural_macro_code(model, block, Symbol("@NLexpression"))
end

###
### @NLparameter(s)
###

"""
    @NLparameter(model, param == value)

Create and return a nonlinear parameter `param` attached to the model `model`
with initial value set to `value`. Nonlinear parameters may be used only in
nonlinear expressions.

## Example

```jldoctest
julia> model = Model();

julia> @NLparameter(model, x == 10)
x == 10.0

julia> value(x)
10.0
```

    @NLparameter(model, value = param_value)

Create and return an anonymous nonlinear parameter `param` attached to the model
`model` with initial value set to `param_value`. Nonlinear parameters may be
used only in nonlinear expressions.

## Example

```jldoctest
julia> model = Model();

julia> x = @NLparameter(model, value = 10)
parameter[1] == 10.0

julia> value(x)
10.0
```

    @NLparameter(model, param_collection[...] == value_expr)

Create and return a collection of nonlinear parameters `param_collection`
attached to the model `model` with initial value set to `value_expr` (may
depend on index sets).
Uses the same syntax for specifying index sets as [`@variable`](@ref).

## Example

```jldoctest
julia> model = Model();

julia> @NLparameter(model, y[i = 1:3] == 2 * i)
3-element Vector{NonlinearParameter}:
 parameter[1] == 2.0
 parameter[2] == 4.0
 parameter[3] == 6.0

julia> value(y[2])
4.0
```

    @NLparameter(model, [...] == value_expr)

Create and return an anonymous collection of nonlinear parameters attached to
the model `model` with initial value set to `value_expr` (may depend on index
sets). Uses the same syntax for specifying index sets as [`@variable`](@ref).

!!! compat
    This macro is part of the legacy nonlinear interface. Consider using the
    new nonlinear interface documented in [Nonlinear Modeling](@ref). In most
    cases, you can replace a call like `@NLparameter(model, p == value)` with
    `@variable(model, p in Parameter(value))`.

## Example

```jldoctest
julia> model = Model();

julia> y = @NLparameter(model, [i = 1:3] == 2 * i)
3-element Vector{NonlinearParameter}:
 parameter[1] == 2.0
 parameter[2] == 4.0
 parameter[3] == 6.0

julia> value(y[2])
4.0
```
"""
macro NLparameter(model, args...)
    esc_m = esc(model)
    error_fn =
        Containers.build_error_fn(:NLparameter, (model, args...), __source__)
    pos_args, kw_args, requested_container = Containers._extract_kw_args(args)
    value = missing
    for arg in kw_args
        if arg.args[1] == :value
            value = arg.args[2]
        end
    end
    kw_args = filter(kw -> kw.args[1] != :value, kw_args)
    if !ismissing(value) && length(pos_args) > 0
        error_fn(
            "Invalid syntax: no positional args allowed for anonymous " *
            "parameters.",
        )
    elseif length(pos_args) > 1
        error_fn("Invalid syntax: too many positional arguments.")
    elseif length(kw_args) > 0
        error_fn("Invalid syntax: unsupported keyword arguments.")
    elseif ismissing(value) && Meta.isexpr(pos_args[1], :block)
        error_fn("Invalid syntax: did you mean to use `@NLparameters`?")
    elseif ismissing(value)
        ex = pos_args[1]
        if !Meta.isexpr(ex, :call) ||
           length(ex.args) != 3 ||
           ex.args[1] != :(==)
            error_fn(
                "Invalid syntax: expected syntax of form `param == value`.",
            )
        end
    end
    param = nothing
    if ismissing(value)
        param, value = pos_args[1].args[2], pos_args[1].args[3]
    end
    name, index_vars, index_values = Containers.parse_ref_sets(
        error_fn,
        param;
        invalid_index_variables = [model],
    )
    code = quote
        if !isa($(esc(value)), Number)
            $(esc(error_fn))("Parameter value is not a number.")
        end
        add_nonlinear_parameter($esc_m, $(esc(value)))
    end
    creation_code = Containers.container_code(
        index_vars,
        index_values,
        code,
        requested_container,
    )
    return _finalize_macro(
        esc_m,
        creation_code,
        __source__;
        register_name = name,
    )
end

"""
     @NLparameters(model, args...)

Create and return multiple nonlinear parameters attached to model `model`, in
the same fashion as [`@NLparameter`](@ref) macro.

The model must be the first argument, and multiple parameters can be added on
multiple lines wrapped in a `begin ... end` block. Distinct parameters need to
be placed on separate lines as in the following example.

The macro returns a tuple containing the parameters that were defined.

!!! compat
    This macro is part of the legacy nonlinear interface. Consider using the
    new nonlinear interface documented in [Nonlinear Modeling](@ref). In most
    cases, you can replace a call like
    ```julia
    @NLparameters(model, begin
        p == value
    end)
    ```
    with
    ```julia
    @variables(model, begin
        p in Parameter(value)
    end)
    ```

## Example

```jldoctest
julia> model = Model();

julia> @NLparameters(model, begin
           x == 10
           b == 156
       end);

julia> value(x)
10.0
```
"""
macro NLparameters(model, block)
    return _plural_macro_code(model, block, Symbol("@NLparameter"))
end

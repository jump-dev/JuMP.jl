#  Copyright 2017, Iain Dunning, Joey Huchette, Miles Lubin, and contributors
#  This Source Code Form is subject to the terms of the Mozilla Public
#  License, v. 2.0. If a copy of the MPL was not distributed with this
#  file, You can obtain one at https://mozilla.org/MPL/2.0/.

"""
    @expression(args...)

Efficiently builds a linear or quadratic expression but does not add to model
immediately. Instead, returns the expression which can then be inserted in other
constraints.

## Example

```jldoctest expression_docstring
julia> model = Model();

julia> @variable(model, x[1:5]);

julia> @variable(model, y);

julia> @variable(model, z);

julia> @expression(model, shared, sum(i * x[i] for i in 1:5))
x[1] + 2 x[2] + 3 x[3] + 4 x[4] + 5 x[5]

julia> @constraint(model, shared + y >= 5)
x[1] + 2 x[2] + 3 x[3] + 4 x[4] + 5 x[5] + y ≥ 5

julia> @constraint(model, shared + z <= 10)
x[1] + 2 x[2] + 3 x[3] + 4 x[4] + 5 x[5] + z ≤ 10
```

The `ref` accepts index sets in the same way as `@variable`, and those indices
can be used in the construction of the expressions:

```jldoctest expression_docstring
julia> @expression(model, expr[i = 1:3], i * sum(x[j] for j in 1:3))
3-element Vector{AffExpr}:
 x[1] + x[2] + x[3]
 2 x[1] + 2 x[2] + 2 x[3]
 3 x[1] + 3 x[2] + 3 x[3]
```

Anonymous syntax is also supported:

```jldoctest expression_docstring
julia> expr = @expression(model, [i in 1:3], i * sum(x[j] for j in 1:3))
3-element Vector{AffExpr}:
 x[1] + x[2] + x[3]
 2 x[1] + 2 x[2] + 2 x[3]
 3 x[1] + 3 x[2] + 3 x[3]
```
"""
macro expression(args...)
    error_fn(str...) = _macro_error(:expression, args, __source__, str...)
    args, kw_args, requested_container = Containers._extract_kw_args(args)
    if length(args) == 3
        m = esc(args[1])
        c = args[2]
        x = args[3]
    elseif length(args) == 2
        m = esc(args[1])
        c = gensym()
        x = args[2]
    else
        error_fn("needs at least two arguments.")
    end
    length(kw_args) == 0 || error_fn("unrecognized keyword argument")
    if Meta.isexpr(args[2], :block)
        error_fn("Invalid syntax. Did you mean to use `@expressions`?")
    end
    anonvar =
        Meta.isexpr(c, :vect) || Meta.isexpr(c, :vcat) || length(args) == 2
    variable = gensym()

    idxvars, indices = Containers.build_ref_sets(error_fn, c)
    if args[1] in idxvars
        error_fn(
            "Index $(args[1]) is the same symbol as the model. Use a " *
            "different name for the index.",
        )
    end
    expr_var, build_code = _rewrite_expression(x)
    code = quote
        $build_code
        # Don't leak a `_MA.Zero` if the expression is an empty summation, or
        # other structure that returns `_MA.Zero()`.
        _replace_zero($m, $expr_var)
    end
    code =
        Containers.container_code(idxvars, indices, code, requested_container)
    # Wrap the entire code block in a let statement to make the model act as
    # a type stable local variable.
    code = _wrap_let(m, code)
    # don't do anything with the model, but check that it's valid anyway
    if anonvar
        macro_code = code
    else
        macro_code = _macro_assign_and_return(
            code,
            variable,
            Containers._get_name(c);
            model_for_registering = m,
        )
    end
    return _finalize_macro(m, macro_code, __source__)
end

"""
    @expressions(model, args...)

Adds multiple expressions to model at once, in the same fashion as the
[`@expression`](@ref) macro.

The model must be the first argument, and multiple expressions can be added on
multiple lines wrapped in a `begin ... end` block.

The macro returns a tuple containing the expressions that were defined.

## Example

```jldoctest
julia> model = Model();

julia> @variable(model, x);

julia> @variable(model, y);

julia> @variable(model, z[1:2]);

julia> a = [4, 5];

julia> @expressions(model, begin
           my_expr, x^2 + y^2
           my_expr_1[i = 1:2], a[i] - z[i]
       end)
(x² + y², AffExpr[-z[1] + 4, -z[2] + 5])
```
"""
macro expressions(model, block)
    return _plural_macro_code(model, block, Symbol("@expression"))
end

#  Copyright 2017, Iain Dunning, Joey Huchette, Miles Lubin, and contributors
#  This Source Code Form is subject to the terms of the Mozilla Public
#  License, v. 2.0. If a copy of the MPL was not distributed with this
#  file, You can obtain one at https://mozilla.org/MPL/2.0/.

"""
    @expression(model::GenericModel, expression)
    @expression(model::GenericModel, [index_sets...], expression)
    @expression(model::GenericModel, name, expression)
    @expression(model::GenericModel, name[index_sets...], expression)

Efficiently builds and returns an expression.

The `name` argument is optional. If index sets are passed, a container is built
and the expression may depend on the indices of the index sets.

## Keyword arguments

 * `container = :Auto`: force the container type by passing `container = Array`,
   `container = DenseAxisArray`, `container = SparseAxisArray`, or any another
   container type which is supported by a JuMP extension.

## Example

```jldoctest
julia> model = Model();

julia> @variable(model, x[1:5]);

julia> @expression(model, shared, sum(i * x[i] for i in 1:5))
x[1] + 2 x[2] + 3 x[3] + 4 x[4] + 5 x[5]

julia> shared
x[1] + 2 x[2] + 3 x[3] + 4 x[4] + 5 x[5]
```

In the same way as [`@variable`](@ref), the second argument may define index
sets, and those indices can be used in the construction of the expressions:

```jldoctest
julia> model = Model();

julia> @variable(model, x[1:3]);

julia> @expression(model, expr[i = 1:3], i * sum(x[j] for j in 1:3))
3-element Vector{AffExpr}:
 x[1] + x[2] + x[3]
 2 x[1] + 2 x[2] + 2 x[3]
 3 x[1] + 3 x[2] + 3 x[3]
```

Anonymous syntax is also supported:

```jldoctest
julia> model = Model();

julia> @variable(model, x[1:3]);

julia> expr = @expression(model, [i in 1:3], i * sum(x[j] for j in 1:3))
3-element Vector{AffExpr}:
 x[1] + x[2] + x[3]
 2 x[1] + 2 x[2] + 2 x[3]
 3 x[1] + 3 x[2] + 3 x[3]
```
"""
macro expression(input_args...)
    error_fn = Containers.build_error_fn(:expression, input_args, __source__)
    args, kwargs = Containers.parse_macro_arguments(
        error_fn,
        input_args;
        num_positional_args = 2:3,
        valid_kwargs = [:container],
    )
    if Meta.isexpr(args[2], :block)
        error_fn("Invalid syntax. Did you mean to use `@expressions`?")
    end
    name_expr = length(args) == 3 ? args[2] : nothing
    name, index_vars, indices = Containers.parse_ref_sets(
        error_fn,
        name_expr;
        invalid_index_variables = [args[1]],
    )
    model = esc(args[1])
    expr, build_code = _rewrite_expression(args[end])
    code = quote
        $build_code
        # Don't leak a `_MA.Zero` if the expression is an empty summation, or
        # other structure that returns `_MA.Zero()`.
        _replace_zero($model, $expr)
    end
    return _finalize_macro(
        model,
        Containers.container_code(index_vars, indices, code, kwargs),
        __source__;
        register_name = name,
        wrap_let = true,
    )
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

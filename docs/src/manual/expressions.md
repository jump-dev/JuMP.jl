```@meta
DocTestSetup = quote
    using JuMP
end
```

# Expressions

JuMP has three types of expressions: affine, quadratic, and nonlinear. These
expressions can be inserted into constraints or into the objective. This is
particularly useful if an expression is used in multiple places in the model.

## Affine expressions

There are four ways of constructing an affine expression in JuMP: with the
[`@expression`](@ref) macro, with operator overloading, with the `AffExpr`
constructor, and with [`add_to_expression!`](@ref).

### Macros

The recommended way to create an affine expression is via the
[`@expression`](@ref) macro.

```jldoctest affine_macro
julia> model = Model();

julia> @variable(model, x)
x

julia> @variable(model, y)
y

julia> ex = @expression(model, 2x + y - 1)
2 x + y - 1
```

This expression can be used in the objective or added to a constraint. For
example:
```jldoctest affine_macro
julia> @objective(model, Min, 2 * ex - 1)
4 x + 2 y - 3

julia> objective_function(model)
4 x + 2 y - 3
```

Just like variables and constraints, named expressions can also be created. For
example
```jldoctest
julia> model = Model();

julia> @variable(model, x[i = 1:3]);

julia> @expression(model, expr[i = 1:3], i * sum(x[j] for j in i:3));

julia> expr
3-element Vector{AffExpr}:
 x[1] + x[2] + x[3]
 2 x[2] + 2 x[3]
 3 x[3]
```

!!! tip
    You can read more about containers in the [Containers](@ref) section.

### Operator overloading

Expressions can also be created without macros. However, note that in some
cases, this can be much slower that constructing an expression using macros.

```jldoctest
julia> model = Model();

julia> @variable(model, x)
x

julia> @variable(model, y)
y

julia> ex = 2x + y - 1
2 x + y - 1
```

### Constructors

A third way to create an affine expression is by the `AffExpr` constructor. The
first argument is the constant term, and the remaining arguments are
variable-coefficient pairs.

```jldoctest
julia> model = Model();

julia> @variable(model, x)
x

julia> @variable(model, y)
y

julia> ex = AffExpr(-1.0, x => 2.0, y => 1.0)
2 x + y - 1
```

### `add_to_expression!`

The fourth way to create an affine expression is by using
[`add_to_expression!`](@ref). Compared to the operator overloading method,
this approach is faster because it avoids constructing temporary objects.
The [`@expression`](@ref) macro uses [`add_to_expression!`](@ref)
behind-the-scenes.
```jldoctest
julia> model = Model();

julia> @variable(model, x)
x

julia> @variable(model, y)
y

julia> ex = AffExpr(-1.0)
-1

julia> add_to_expression!(ex, 2.0, x)
2 x - 1

julia> add_to_expression!(ex, 1.0, y)
2 x + y - 1
```

!!! warning
    Read the section [Initializing arrays](@ref) for some cases to be careful
    about when using [`add_to_expression!`](@ref).

### Removing zero terms

Use [`drop_zeros!`](@ref) to remove terms from an affine expression with a `0`
coefficient.

```jldoctest
julia> model = Model();

julia> @variable(model, x)
x

julia> @expression(model, ex, x + 1 - x)
0 x + 1

julia> drop_zeros!(ex)

julia> ex
1
```

### Coefficients

Use [`coefficient`](@ref) to return the coefficient associated with a variable
in an affine expression.

```jldoctest
julia> model = Model();

julia> @variable(model, x)
x

julia> @variable(model, y)
y

julia> @expression(model, ex, 2x + 1)
2 x + 1

julia> coefficient(ex, x)
2.0

julia> coefficient(ex, y)
0.0
```

## Quadratic expressions

Like affine expressions, there are four ways of constructing a quadratic
expression in JuMP: macros, operator overloading, constructors, and
[`add_to_expression!`](@ref).

### Macros

The [`@expression`](@ref) macro can be used to create quadratic expressions by
including quadratic terms.

```jldoctest
julia> model = Model();

julia> @variable(model, x)
x

julia> @variable(model, y)
y

julia> ex = @expression(model, x^2 + 2 * x * y + y^2 + x + y - 1)
x² + 2 x*y + y² + x + y - 1
```

### Operator overloading

Operator overloading can also be used to create quadratic expressions. The same
performance warning (discussed in the affine expression section) applies.

```jldoctest
julia> model = Model();

julia> @variable(model, x)
x

julia> @variable(model, y)
y

julia> ex = x^2 + 2 * x * y + y^2 + x + y - 1
x² + 2 x*y + y² + x + y - 1
```

### Constructors

Quadratic expressions can also be created using the `QuadExpr` constructor. The
first argument is an affine expression, and the remaining arguments are pairs,
where the first term is a `JuMP.UnorderedPair` and the second term is the
coefficient.

```jldoctest
julia> model = Model();

julia> @variable(model, x)
x

julia> @variable(model, y)
y

julia> aff_expr = AffExpr(-1.0, x => 1.0, y => 1.0)
x + y - 1

julia> quad_expr = QuadExpr(
           aff_expr,
           UnorderedPair(x, x) => 1.0,
           UnorderedPair(x, y) => 2.0,
           UnorderedPair(y, y) => 1.0,
       )
x² + 2 x*y + y² + x + y - 1
```

### `add_to_expression!`

Finally, [`add_to_expression!`](@ref) can also be used to add quadratic terms.

```jldoctest
julia> model = Model();

julia> @variable(model, x)
x

julia> @variable(model, y)
y

julia> ex = QuadExpr(x + y - 1.0)
x + y - 1

julia> add_to_expression!(ex, 1.0, x, x)
x² + x + y - 1

julia> add_to_expression!(ex, 2.0, x, y)
x² + 2 x*y + x + y - 1

julia> add_to_expression!(ex, 1.0, y, y)
x² + 2 x*y + y² + x + y - 1
```

!!! warning
    Read the section [Initializing arrays](@ref) for some cases to be careful
    about when using [`add_to_expression!`](@ref).

### Removing zero terms

Use [`drop_zeros!`](@ref) to remove terms from a quadratic expression with a `0`
coefficient.

```jldoctest
julia> model = Model();

julia> @variable(model, x)
x

julia> @expression(model, ex, x^2 + x + 1 - x^2)
0 x² + x + 1

julia> drop_zeros!(ex)

julia> ex
x + 1
```

### Coefficients

Use [`coefficient`](@ref) to return the coefficient associated with a pair of variables
in a quadratic expression.

```jldoctest
julia> model = Model();

julia> @variable(model, x)
x

julia> @variable(model, y)
y

julia> @expression(model, ex, 2*x*y + 3*x)
2 x*y + 3 x

julia> coefficient(ex, x, y)
2.0

julia> coefficient(ex, x, x)
0.0

julia> coefficient(ex, y, x)
2.0

julia> coefficient(ex, x)
3.0
```

## Nonlinear expressions

Nonlinear expressions can be constructed only using the [`@NLexpression`](@ref)
macro and can be used only in [`@NLobjective`](@ref), [`@NLconstraint`](@ref),
and other [`@NLexpression`](@ref)s. For more details, see the [Nonlinear
Modeling](@ref) section.

## Initializing arrays

JuMP implements `zero(AffExpr)` and `one(AffExpr)` to support various functions
in `LinearAlgebra` (for example, accessing the off-diagonal of a `Diagonal`
matrix).
```jldoctest
julia> zero(AffExpr)
0

julia> one(AffExpr)
1
```

However, this can result in a subtle bug if you call
[`add_to_expression!`](@ref) or the [MutableArithmetics API](https://github.com/jump-dev/MutableArithmetics.jl)
on an element created by `zeros` or `ones`:
```jldoctest
julia> x = zeros(AffExpr, 2)
2-element Vector{AffExpr}:
 0
 0

julia> add_to_expression!(x[1], 1.1)
1.1

julia> x
2-element Vector{AffExpr}:
 1.1
 1.1
```

Notice how we modified `x[1]`, but we also changed `x[2]`!

This happened because `zeros(AffExpr, 2)` calls `zero(AffExpr)` once to obtain a
zero element, and then creates an appropriately sized array filled with the same
element.

This also happens with broadcasting calls containing a conversion of `0` or `1`:
```jldoctest
julia> x = Vector{AffExpr}(undef, 2)
2-element Vector{AffExpr}:
 #undef
 #undef

julia> x .= 0
2-element Vector{AffExpr}:
 0
 0

julia> add_to_expression!(x[1], 1.1)
1.1

julia> x
2-element Vector{AffExpr}:
 1.1
 1.1
```

The recommended way to create an array of empty expressions is as follows:
```jldoctest
julia> x = Vector{AffExpr}(undef, 2)
2-element Vector{AffExpr}:
 #undef
 #undef

julia> for i in eachindex(x)
           x[i] = AffExpr(0.0)
       end

julia> add_to_expression!(x[1], 1.1)
1.1

julia> x
2-element Vector{AffExpr}:
 1.1
 0
```

Alternatively, use non-mutating operation to avoid updating `x[1]` in-place:
```jldoctest
julia> x = zeros(AffExpr, 2)
2-element Vector{AffExpr}:
 0
 0

julia> x[1] += 1.1
1.1

julia> x
2-element Vector{AffExpr}:
 1.1
 0
```
Note that for large expressions this will be slower due to the allocation of
additional temporary objects.

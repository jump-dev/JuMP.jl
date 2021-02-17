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
model = Model()
@variable(model, x)
@variable(model, y)
ex = @expression(model, 2x + y - 1)

# output

2 x + y - 1
```

This expression can be used in the objective or added to a constraint. For
example:
```jldoctest affine_macro
@objective(model, Min, 2 * ex - 1)
objective_function(model)

# output

4 x + 2 y - 3
```

Just like variables and constraints, named expressions can also be created. For
example
```jldoctest
model = Model()
@variable(model, x[i = 1:3])
@expression(model, expr[i = 1:3], i * sum(x[j] for j in i:3))
expr

# output

3-element Array{GenericAffExpr{Float64,VariableRef},1}:
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
model = Model()
@variable(model, x)
@variable(model, y)
ex = 2x + y - 1

# output

2 x + y - 1
```

### Constructors

A third way to create an affine expression is by the `AffExpr` constructor. The
first argument is the constant term, and the remaining arguments are
variable-coefficient pairs.

```jldoctest
model = Model()
@variable(model, x)
@variable(model, y)
ex = AffExpr(-1.0, x => 2.0, y => 1.0)

# output

2 x + y - 1
```

### `add_to_expression!`

The fourth way to create an affine expression is by using
[`add_to_expression!`](@ref). Compared to the operator overloading method,
this approach is faster because it avoids constructing temporary objects.
The [`@expression`](@ref) macro uses [`add_to_expression!`](@ref)
behind-the-scenes.
```jldoctest
model = Model()
@variable(model, x)
@variable(model, y)
ex = AffExpr(-1.0)
add_to_expression!(ex, 2.0, x)
add_to_expression!(ex, 1.0, y)

# output

2 x + y - 1
```

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
model = Model()
@variable(model, x)
@variable(model, y)
ex = @expression(model, x^2 + 2 * x * y + y^2 + x + y - 1)

# output

x² + 2 y*x + y² + x + y - 1
```

### Operator overloading

Operator overloading can also be used to create quadratic expressions. The same
performance warning (discussed in the affine expression section) applies.

```jldoctest
model = Model()
@variable(model, x)
@variable(model, y)
ex = x^2 + 2 * x * y + y^2 + x + y - 1

# output

x² + 2 x*y + y² + x + y - 1
```

### Constructors

Quadratic expressions can also be created using the `QuadExpr` constructor. The
first argument is an affine expression, and the remaining arguments are pairs,
where the first term is a `JuMP.UnorderedPair` and the second term is the
coefficient.

```jldoctest
model = Model()
@variable(model, x)
@variable(model, y)
aff_expr = AffExpr(-1.0, x => 1.0, y => 1.0)
quad_expr = QuadExpr(aff_expr, UnorderedPair(x, x) => 1.0,
                     UnorderedPair(x, y) => 2.0, UnorderedPair(y, y) => 1.0)

# output

x² + 2 x*y + y² + x + y - 1
```

### `add_to_expression!`

Finally, [`add_to_expression!`](@ref) can also be used to add quadratic terms.

```jldoctest
model = Model()
@variable(model, x)
@variable(model, y)
ex = QuadExpr(x + y - 1.0)
add_to_expression!(ex, 1.0, x, x)
add_to_expression!(ex, 2.0, x, y)
add_to_expression!(ex, 1.0, y, y)

# output

x² + 2 x*y + y² + x + y - 1
```

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
and other [`@NLexpression`](@ref)s. Moreover, quadratic and affine expressions
cannot be used in the nonlinear macros. For more details, see the [Nonlinear
Modeling](@ref) section.

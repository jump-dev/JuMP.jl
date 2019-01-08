```@meta
DocTestSetup = quote
    using JuMP
end
```

# Expressions

JuMP has three types of expressions: affine, quadratic, and nonlinear.

## Affine expressions

An affine expression is

There are four ways of constructing an affine expression in JuMP: using macros,
operator overloading, constructors, and via `JuMP.add_to_expression!`.

### Macros

```jldoctest
model = Model()
@variable(model, x)
@variable(model, y)
ex = @expression(model, 2x + y - 1)

# output

2 x + y - 1
```

### Operator overloading

```jldoctest
model = Model()
@variable(model, x)
@variable(model, y)
ex = 2x + y - 1

# output

2 x + y - 1
```

### Constructors

```jldoctest
model = Model()
@variable(model, x)
@variable(model, y)
ex = AffExpr(-1.0, x => 2.0, y => 1.0)

# output

2 x + y - 1
```

### `add_to_expression!`

```jldoctest
model = Model()
@variable(model, x)
@variable(model, y)
ex = JuMP.AffExpr(-1.0)
JuMP.add_to_expression!(ex, 2.0, x)
JuMP.add_to_expression!(ex, 1.0, y)

# output

2 x + y - 1
```

## Quadratic expressions

A quadratic expression is

Like affine expressions, there are four ways of constructing a quadratic
expression in JuMP: using macros, operator overloading, constructors, and via
`JuMP.add_to_expression!`.

### Macros

```jldoctest
model = Model()
@variable(model, x)
@variable(model, y)
ex = @expression(model, x^2 + 2 * x * y + y^2 + x + y - 1)

# output

x² + 2 x*y + y² + x + y - 1
```

### Operator overloading

```jldoctest
model = Model()
@variable(model, x)
@variable(model, y)
ex = x^2 + 2 * x * y + y^2 + x + y - 1

# output

x² + 2 x*y + y² + x + y - 1
```

### Constructors

```jldoctest
model = Model()
@variable(model, x)
@variable(model, y)
aff_expr = AffExpr(-1.0, x => 1.0, y => 1.0)
quad_expr = JuMP.QuadExpr(aff_expr,
    JuMP.UnorderedPair(x, x) => 1.0,
    JuMP.UnorderedPair(x, y) => 2.0,
    JuMP.UnorderedPair(y, y) => 1.0
)

# output

x² + 2 x*y + y² + x + y - 1
```

### `add_to_expression!`

```jldoctest
model = Model()
@variable(model, x)
@variable(model, y)
ex = JuMP.QuadExpr(x + y - 1.0)
JuMP.add_to_expression!(ex, 1.0, x, x)
JuMP.add_to_expression!(ex, 2.0, x, y)
JuMP.add_to_expression!(ex, 1.0, y, y)

# output

x² + 2 x*y + y² + x + y - 1
```

## Nonlinear expressions

Nonlinear expressions can only be constructed using the [`@NLexpression`](@ref)
macro. For more details, see the [Nonlinear Modelling](@ref) section.

## Reference

```@docs
@expression
@NLexpression
JuMP.add_to_expression!
```

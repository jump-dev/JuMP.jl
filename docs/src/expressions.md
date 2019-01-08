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

There are four ways of constructing an affine expression in JuMP: using the
[`@expression`](@ref) macro, by operator overloading, using the `AffExpr`
constructor, and via [`JuMP.add_to_expression!`](@ref).

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

This expression can be used in the objective, or added to a constraint. For
example:
```jldoctest affine_macro
@objective(model, Min, 2 * ex - 1)
JuMP.objective_function(model)

# output

4 x + 2 y - 3
```

Just like variables and constraints, named expressions can also be created. For
example
```jldoctest
model = Model()
@variable(model, x[i = 1:3])
@expression(model, expr[i=1:3], i * sum(x[j] for j in i:3))
expr

# output

3-element Array{JuMP.GenericAffExpr{Float64,VariableRef},1}:
 x[1] + x[2] + x[3]
 2 x[2] + 2 x[3]
 3 x[3]
```

### Operator overloading

Expressions can also be created outside the macro. However, not that this is
much slower that constructing an expression using the macro. This should only be
used if the performance is not critical.

```jldoctest
model = Model()
@variable(model, x)
@variable(model, y)
ex = 2x + y - 1

# output

2 x + y - 1
```

### Constructors

A third was to create an affine expression is by the `AffExpr` constructor. The
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

### `JuMP.add_to_expression!`

The fourth way to create an affine expression is by using
[`JuMP.add_to_expression!`](@ref). Compared to the operator overloading method,
this approach is fast, and is what the [`@expression`](@ref) macro implements
behind-the-scenes.
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

Like affine expressions, there are four ways of constructing a quadratic
expression in JuMP: using macros, operator overloading, constructors, and via
[`JuMP.add_to_expression!`](@ref).

### Macros

The [`@expression`](@ref) macro can be used to create quadratic expressions by
including quadratic terms.

```jldoctest
model = Model()
@variable(model, x)
@variable(model, y)
ex = @expression(model, x^2 + 2 * x * y + y^2 + x + y - 1)

# output

x² + 2 x*y + y² + x + y - 1
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
quad_expr = JuMP.QuadExpr(aff_expr,
    JuMP.UnorderedPair(x, x) => 1.0,
    JuMP.UnorderedPair(x, y) => 2.0,
    JuMP.UnorderedPair(y, y) => 1.0
)

# output

x² + 2 x*y + y² + x + y - 1
```

### `JuMP.add_to_expression!`

Finally, [`JuMP.add_to_expression!`](@ref) can also be used to add quadratic
terms.

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
macro. For more details, see the [Nonlinear Modeling](@ref) section.

## Reference

```@docs
@expression
JuMP.add_to_expression!
```

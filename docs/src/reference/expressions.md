# [Expressions](@id ExpressionAPI)

More information can be found in the [Expressions](@ref) section of the manual.


## Macros

```@docs
@expression
@expressions
```

## Affine expressions

```@docs
GenericAffExpr
AffExpr
linear_terms
```

## Quadratic expressions


```@docs
GenericQuadExpr
QuadExpr
UnorderedPair
quad_terms
```

## Utilities and modifications

```@docs
constant
coefficient
isequal_canonical
add_to_expression!
drop_zeros!
map_coefficients
map_coefficients_inplace!
```

## JuMP-to-MOI converters

```@docs
variable_ref_type
jump_function
jump_function_type
moi_function
moi_function_type
```

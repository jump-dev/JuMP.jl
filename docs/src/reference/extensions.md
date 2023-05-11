# [Extensions](@id ExtensionAPI)

More information can be found in the [Extensions](@ref extensions_manual)
section of the manual.

## Define a new set

```@docs
AbstractVectorSet
AbstractScalarSet
```

## Extend `@variable`

```@docs
ScalarVariable
VariableInfo
add_variable
build_variable
```

## Extend `@constraint`

```@docs
build_constraint
add_constraint
model_convert
AbstractShape
shape
reshape_vector
reshape_set
dual_shape
ScalarShape
VectorShape
SquareMatrixShape
SymmetricMatrixShape
operator_to_set
parse_constraint
parse_constraint_head
parse_constraint_call
```

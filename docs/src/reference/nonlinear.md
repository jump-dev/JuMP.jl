# [Nonlinear Modeling](@id NonlinearAPI)

More information can be found in the [Nonlinear](@ref nonlinear_developers)
section of the manual.

```@docs
Nonlinear.NonlinearData
```

## [Expressions](@id nonlinear_api_expressions)

```@docs
Nonlinear.ExpressionIndex
Nonlinear.add_expression
```

## [Parameters](@id nonlinear_api_parameters)

```@docs
Nonlinear.ParameterIndex
Nonlinear.add_parameter
Nonlinear.set_parameter
```

## [Objectives](@id nonlinear_api_objectives)

```@docs
Nonlinear.set_objective
```

## [Constraints](@id nonlinear_api_constraints)

```@docs
Nonlinear.ConstraintIndex
Nonlinear.add_constraint
Nonlinear.delete
```

## [User-defined operators](@id nonlinear_api_operators)

```@docs
Nonlinear.register_operator
```

## Automatic-differentiation backends

```@docs
Nonlinear.AbstractAutomaticDifferentiation
Nonlinear.Default
Nonlinear.set_differentiation_backend
```

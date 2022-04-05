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
Nonlinear.OperatorRegistry
Nonlinear.DEFAULT_UNIVARIATE_OPERATORS
Nonlinear.DEFAULT_MULTIVARIATE_OPERATORS
Nonlinear.register_operator
Nonlinear.register_operator_if_needed
Nonlinear.assert_registered
Nonlinear.check_return_type
Nonlinear.eval_univariate_function
Nonlinear.eval_univariate_gradient
Nonlinear.eval_univariate_hessian
Nonlinear.eval_multivariate_function
Nonlinear.eval_multivariate_gradient
Nonlinear.eval_logic_function
Nonlinear.eval_comparison_function
```

## Automatic-differentiation backends

```@docs
Nonlinear.AbstractAutomaticDifferentiation
Nonlinear.Default
Nonlinear.SparseReverseMode
Nonlinear.set_differentiation_backend
```

## Data-structure

```@docs
Nonlinear.Node
Nonlinear.NodeType
Nonlinear.NonlinearExpression
Nonlinear.adjacency_matrix
```

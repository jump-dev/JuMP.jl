# [Nonlinear Modeling](@id NonlinearAPI)

More information can be found in the [Nonlinear Modeling](@ref) section of the
manual.

## [Constraints](@id ref_nl_constraints)

```@docs
@NLconstraint
@NLconstraints
NonlinearConstraintIndex
num_nonlinear_constraints
add_nonlinear_constraint
all_nonlinear_constraints
nonlinear_dual_start_value
set_nonlinear_dual_start_value
```

## [Expressions](@id ref_nl_expressions)

```@docs
@NLexpression
@NLexpressions
NonlinearExpression
add_nonlinear_expression
```

## [Objectives](@id ref_nl_objectives)

```@docs
@NLobjective
set_nonlinear_objective
```

## [Parameters](@id ref_nl_parameters)

```@docs
@NLparameter
@NLparameters
NonlinearParameter
value(::JuMP.NonlinearParameter)
set_value(::JuMP.NonlinearParameter, ::Number)
```

## User-defined functions

```@docs
register
```

## Derivatives

```@docs
NLPEvaluator
```

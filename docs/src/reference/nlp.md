# [Nonlinear Modeling](@id NonlinearAPI)

More information can be found in the [Nonlinear Modeling](@ref) section of the
manual.

## [Constraints](@id ref_nl_constraints)

```@docs
@NLconstraint
@NLconstraints
NonlinearConstraintIndex
num_nl_constraints
add_NL_constraint
all_nl_constraints
nl_dual_start_value
set_nl_dual_start_value
```

## [Expressions](@id ref_nl_expressions)

```@docs
@NLexpression
@NLexpressions
NonlinearExpression
add_NL_expression
```

## [Objectives](@id ref_nl_objectives)

```@docs
@NLobjective
set_NL_objective
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

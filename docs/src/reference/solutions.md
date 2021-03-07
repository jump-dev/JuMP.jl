# [Solutions](@id SolutionAPI)

More information can be found in the [Solutions](@ref) section of the manual.

## Basic utilities

```@docs
JuMP.optimize!
NoOptimizer
OptimizeNotCalled
solution_summary
```

## Termination status

```@docs
termination_status
MOI.TerminationStatusCode
raw_status
result_count
```

## Primal solutions

```@docs
primal_status
has_values
value
MOI.ResultStatusCode
```

## Dual solutions

```@docs
dual_status
has_duals
dual
shadow_price
reduced_cost
```

## Basic attributes

```@docs
objective_value
objective_bound
dual_objective_value
solve_time
relative_gap
simplex_iterations
barrier_iterations
node_count
```

## [Conflicts](@id ref_conflicts)

```@docs
JuMP.compute_conflict!
MOI.compute_conflict!
MOI.ConflictStatus
MOI.ConflictStatusCode
MOI.ConstraintConflictStatus
MOI.ConflictParticipationStatusCode
JuMP.copy_conflict
```

## Sensitivity

```@docs
lp_sensitivity_report
SensitivityReport
lp_objective_perturbation_range
lp_rhs_perturbation_range
```

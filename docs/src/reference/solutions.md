# [Solutions](@id SolutionAPI)

More information can be found in the [Querying Solutions](@ref) section of the
manual.

```@docs
JuMP.optimize!
NoOptimizer

termination_status
MOI.TerminationStatusCode

raw_status
primal_status
MOI.ResultStatusCode

result_count

has_values
value

dual_status
has_duals
dual
shadow_price
reduced_cost

objective_bound
objective_value
dual_objective_value

solve_time
relative_gap
simplex_iterations
barrier_iterations
node_count

OptimizeNotCalled
MOI.optimize!

JuMP.compute_conflict!
MOI.compute_conflict!
MOI.ConflictStatus
MOI.ConflictStatusCode
MOI.ConstraintConflictStatus
MOI.ConflictParticipationStatusCode

lp_objective_perturbation_range
lp_rhs_perturbation_range
lp_sensitivity_report
SensitivityReport
```

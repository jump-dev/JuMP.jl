Objective
=========

This page describes macros and functions related to linear and quadratic
objective functions only, unless otherwise indicated. For nonlinear objective
functions, see [Nonlinear Modeling](@ref).

Use the [`@objective`](@ref) macro to set linear and quadratic objective
functions in a JuMP model. The functions [`set_objective_sense`](@ref) and
[`set_objective_function`](@ref) provide an equivalent lower-level interface.

To update a term in the objective function, see
[`set_objective_coefficient`](@ref). 

To query the objective function from a model, see [`objective_sense`](@ref),
[`objective_function`](@ref), and [`objective_function_type`](@ref).

To query the optimal objective value or best known bound after a solve, see
[`objective_value`](@ref) and [`objective_bound`](@ref). These two functions
apply to nonlinear objectives also. The optimal value of the
[dual](@ref constraint_duality) objective can be obtained via
[`dual_objective_value`](@ref).


## Reference

```@docs
@objective
JuMP.set_objective_sense
JuMP.set_objective_function
JuMP.set_objective_coefficient

JuMP.objective_sense
JuMP.objective_function
JuMP.objective_function_type

JuMP.objective_bound
JuMP.objective_value
JuMP.dual_objective_value
```

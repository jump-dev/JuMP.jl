Objective
=========

TODO: Describe how the objective is represented (link to MOI docs)

Objective functions
-------------------

TODO: Describe how JuMP expressions relate to MOI functions. How to set, query,
and modify an objective function.

Setting the objective function and objective sense:
```@docs
@objective
JuMP.set_objective_sense
JuMP.set_objective_function
```

Querying the objective function and objective sense:
```@docs
JuMP.objective_sense
JuMP.objective_function
JuMP.objective_function_type
```

Querying the objective value and bound:
```@docs
JuMP.objective_bound
JuMP.objective_value
```

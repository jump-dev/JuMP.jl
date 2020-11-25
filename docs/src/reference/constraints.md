# Constraint API

```@docs
@constraint
@SDconstraint

name(::JuMP.ConstraintRef{Model, <:JuMP.MOI.ConstraintIndex})
set_name(::JuMP.ConstraintRef{Model, <:JuMP.MOI.ConstraintIndex}, ::String)
constraint_by_name

normalized_coefficient
set_normalized_coefficient

normalized_rhs
set_normalized_rhs

add_to_function_constant

is_valid
JuMP.delete

list_of_constraint_types

all_constraints

num_constraints

constraint_object

set_dual_start_value
dual_start_value

ConstraintRef

SecondOrderCone
RotatedSecondOrderCone
PSDCone

AbstractConstraint
ScalarConstraint
VectorConstraint

index(::ConstraintRef)
optimizer_index(::ConstraintRef{Model})
```

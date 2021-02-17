# [Constraints](@id ConstraintAPI)

More information can be found in the [Constraints](@ref) section of the manual.

```@docs
@constraint
@constraints
@SDconstraint
@SDconstraints

name(::ConstraintRef{Model,<:JuMP._MOICON})
set_name(::ConstraintRef{Model,<:JuMP._MOICON}, ::String)
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

set_dual_start_value
dual_start_value

ConstraintRef

SecondOrderCone
RotatedSecondOrderCone
PSDCone
SOS1
SOS2
SkewSymmetricMatrixSpace
SkewSymmetricMatrixShape

AbstractConstraint
ScalarConstraint
VectorConstraint

index(::ConstraintRef)
optimizer_index(::ConstraintRef{Model})

constraint_object
moi_set

function_string
constraints_string
in_set_string
show_constraints_summary

ConstraintNotOwned
```

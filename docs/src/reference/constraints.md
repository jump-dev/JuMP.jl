# [Constraints](@id ConstraintAPI)

More information can be found in the [Constraints](@ref jump_constraints)
section of the manual.

## Macros

```@docs
@constraint
@constraints
@SDconstraint
@SDconstraints
ConstraintRef
AbstractConstraint
ScalarConstraint
VectorConstraint
```

## Names

```@docs
name(::ConstraintRef{Model,<:JuMP._MOICON})
set_name(::ConstraintRef{Model,<:JuMP._MOICON}, ::String)
constraint_by_name
```

## Modification

```@docs
normalized_coefficient
set_normalized_coefficient

normalized_rhs
set_normalized_rhs

add_to_function_constant
```

## Deletion

```@docs
JuMP.delete
is_valid
ConstraintNotOwned
```

## Query constraints

```@docs
list_of_constraint_types
all_constraints
num_constraints
index(::ConstraintRef)
optimizer_index(::ConstraintRef{Model})
constraint_object
```

## Start values

```@docs
set_dual_start_value
dual_start_value
```

## Special sets

```@docs
SecondOrderCone
RotatedSecondOrderCone
PSDCone
SOS1
SOS2
SkewSymmetricMatrixSpace
SkewSymmetricMatrixShape
moi_set
```

## Printing

```@docs
function_string
constraints_string
in_set_string
show_constraints_summary
```

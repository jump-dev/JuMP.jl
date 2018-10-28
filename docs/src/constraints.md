Constraints
===========

DRAFT: Describe how constraints are represented (link to MOI docs). Constraints
are very similar to variables in (1) how names work (2) how attributes work, and
(3) the macro syntax for constructing them. They're a bit different because
they're parameterized by function-set type. Describe constraints vs.
`ConstraintRefs`. Describe `JuMP.constraint_object`. How to delete constraints.
How to modify constraints by setting attributes and `MOI.modifyconstraint!`.
Describe semidefinite constraints and symmetry handling. Refer to NLP docs for
nonlinear constraints.

```@docs
@constraint
@SDconstraint
```

## Sets

As mentioned in the documentation of the [`@constraint`](@ref) and
[`@SDconstraint`](@ref) macros, the following sets can be used to create
constraints in addition to [any MOI set](http://www.juliaopt.org/MathOptInterface.jl/stable/apireference.html#Sets-1).

```@docs
SecondOrderCone
RotatedSecondOrderCone
PSDCone
```

## Constraint modifications

`JuMP.set_coefficient(constraint, variable, value)`

```@docs
JuMP.set_coefficient
```

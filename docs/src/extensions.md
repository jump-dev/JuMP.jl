```@meta
CurrentModule = JuMP
```

Extending JuMP
==============

TODO: How to extend JuMP: discussion on different ways to build on top of JuMP.
How to extend JuMP's macros and how to avoid doing this.

## Extending MOI

TODO: Create new MOI function/sets, how to use it in JuMP, create new bridge

## Extending JuMP macros

In order to provide a convenient syntax for the user to create variables,
constraints or set the objective of a JuMP extension, it might be required to
use macros similar to [`@variable`](@ref), [`@constraint`](@ref) and
[`@objective`](@ref).
It is recommended to first check whether it is possible to extend one of these
three macros before creating a new one so as to leverage all their features and
provide a more consistent interface to the user.

## Extending the `@variable` macro

TODO: parse/build/add

## Extending the `@constraint` macro

The [`@constraint`](@ref) macro always calls the same three functions:
* `parse_constraint`: is called at parsing time, it parses the constraint
  expression and returns a [`build_constraint`](@ref) call expression;
* [`build_constraint`](@ref): given the functions and sets involved in the
  constraints, it returns a `AbstractConstraint`;
* [`add_constraint`](@ref): given the model, the `AbstractConstraint`
  constructedi in [`build_constraint`](@ref) and the constraint name, it stores
  them in the model and returns a `ConstraintRef`.

Adding methods to these functions is the recommended way to extend the
[`@constraint`](@ref) macro.

### Adding `build_constraint` methods

There is typically two choices when creating a [`build_constraint`](@ref)
method, either return an `AbstractConstraint` already supported by the
model, i.e. `ScalarConstraint` or `VectorConstraint`, or a custom
`AbstractConstraint` with a corresponding [`add_constraint`](@ref) method (see
[Adding `add_constraint` methods](@ref)).

```@docs
build_constraint
```

#### Shapes

Shapes allow vector constraints, which are represented as flat vectors in MOI,
to retain a matrix shape at the JuMP level. There is a `shape` field in
`VectorConstraint` that can be set in [`build_constraint`](@ref) and that is
used to reshape the result computed in [`value`](@ref) and [`dual`](@ref).

```@docs
AbstractShape
shape
reshape
dual_shape
ScalarShape
VectorShape
SquareMatrixShape
SymmetricMatrixShape
```

### Adding `add_constraint` methods

TODO

```@docs
add_constraint
```

### Adding `parse_constraint` methods

TODO

## Extending the [`@objective`](@ref) macro

TODO

## Defining new JuMP models

TODO

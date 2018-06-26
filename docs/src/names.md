Names
=====

There a two different naming aspects that need to be distinguished when
creating variables/contraints (resp. a container of variables/constraints):

* The name of the local variable created (if any) holding the reference (resp.
  the container of references) which corresponds to the name that can be used
  to retrieve it using `m[:name]`.
* The name of the variable/constraint (resp. each variable/constraint in the
  container) used for printing. This corresponds to the
  `MOI.VariableName`/`MOI.ConstraintName` attribute.

When creating a variable using the syntax `@variable(m; kwargs...)`, creating
a constraint using the syntax `@constraint(m, expr)` or when creating a
container with the syntax `[...]` in a macro, we say that the variable or
constraint is *anonymous*.
For anonymous variables/constraints, no local variable is created holding the
reference or container of references and it is not stored in the model, i.e. it
is not possible to retrieve it using `m[:name]`.

Otherwise, when it is not anonymous, the name used both for the local variable
created and the key for retrieving the reference or container of references in
the model are determined from the macro expression. For instance, when creating
a container with the syntax `name[...]` or when creating a constraint with
`@constraint(m, name, expr)`, the name used is `name`.

The name of the variable/constraint used for printing is based on the base name
which is specified by the `basename` keyword argument. When the `basename`
keyword argument is not specified, the name depends on whether the variable is
anonymous:

* if the variable/constraint is anonymous, then the
  `MOI.VariableName`/`MOI.ConstraintName` attribute is not set and the name
  used for printing is `noname`,
* otherwise, the base name is set to the name used for the local variable
  created.

The name of the variable/constraint set to the
`MOI.VariableName`/`MOI.ConstraintName` attribute and used for printing is then
`basename` for single variable/constraint and `basename[i1,i2,...,in]` for the
reference at indices `i1`, `i2`, ..., `in` in a container.

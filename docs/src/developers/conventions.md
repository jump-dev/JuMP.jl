# Conventions for interfacing between JuMP and MathOptInterface

The purpose of this guide is to document the conventions that we have developed for
interfacing between JuMP and MathOptInterface.

## Attributes

JuMP provides [`get_attribute`](@ref) and [`set_attribute`](@ref) as thin shims
for [`MOI.get`](@ref) and [`MOI.set`](@ref). However, there are two common cases
where the thin shims are not sufficient:

 * when the value of the attribute is an [`AbstractJuMPScalar`](@ref) that needs
   to be mapped to and from [`MOI.AbstractFunction`](@ref)
 * when the value of the attribute depends on the shape of the constraint, for
   example, when getting [`MOI.ConstraintDualStart`](@ref) of matrix-valued
   constraints.

In these two cases, the convention is to keep [`get_attribute`](@ref) as a thin
shim that does not modify the attribute value, and to develop new functions that
modify or reshape the value as appropriate.

As an example, JuMP provides [`dual_start_value`](@ref) and [`set_dual_start_value`](@ref)
to get and set the [`MOI.ConstraintDualStart`](@ref) in the original matrix
shape, while [`get_attribute`](@ref) and [`set_attribute`](@ref) can be used to
get and set the value in the vectorized shape:

```@repl
using JuMP
model = Model();
@variable(model, x[1:2, 1:2], PSD);
c = VariableInSetRef(x);
set_dual_start_value(c, [1 0; 0 1])
dual_start_value(c)
get_attribute(c, MOI.ConstraintDualStart())
set_attribute(c, MOI.ConstraintDualStart(), [2.0, -1.0, 1.0])
dual_start_value(c)
```

## `unset_` methods

There are a variety of attributes in JuMP and MOI that can be "set" and "unset."
For example, there is [`MOI.Silent`](@ref), and the corresponding
[`set_silent`](@ref) and [`unset_silent`](@ref).

Note how [`set_silent`](@ref) and [`unset_silent`](@ref) take a single argument
(the model), where `set_silent(model)` corresponds to `MOI.set(model, MOI.Silent(), true)`
and `unset_silent(model)` corresponds to `MOI.set(model, MOI.Silent(), false)`.
We could have instead implemented a single method `set_silent(model, flag::Bool)`
that corresponded to `MOI.set(model, MOI.Silent(), flag)`. Another example is
[`unset_time_limit_sec`](@ref), which is equivalent to
`set_time_limit_sec(model, nothing)`.

We have come to regard the `unset_` design as a mistake, because it leads to a
proliferation of unique function names instead of leveraging Julia's strength
for multiple dispatch.

The existing `unset_` names are retained for backwards compatibility, but, going
forward, provide a single `set_` method and document what value type should be
provided to restore the model to the default setting. Thus, we have
[`set_string_names_on_creation`](@ref), but no corresponding
`unset_string_names_on_creation`.

# [Variables](@id VariableAPI)

More information can be found in the [Variables](@ref jump_variables) section of
the manual.

## Macros

```@docs
@variable
@variables
```

## Basic utilities

```@docs
VariableRef
num_variables
all_variables
owner_model
index(::GenericVariableRef)
optimizer_index(::GenericVariableRef)
check_belongs_to_model
VariableNotOwned
VariableConstrainedOnCreation
VariablesConstrainedOnCreation
ComplexPlane
ComplexVariable
```

## Names

```@docs
name(::JuMP.VariableRef)
set_name(::JuMP.VariableRef, ::String)
variable_by_name
```

## Start values

```@docs
has_start_value
set_start_value
start_value
```

## Lower bounds

```@docs
has_lower_bound
lower_bound
set_lower_bound
delete_lower_bound
LowerBoundRef
```

## Upper bounds

```@docs
has_upper_bound
upper_bound
set_upper_bound
delete_upper_bound
UpperBoundRef
```

## Fixed bounds

```@docs
is_fixed
fix_value
fix
unfix
FixRef
```

## Integer variables

```@docs
is_integer
set_integer
unset_integer
IntegerRef
```

## Binary variables

```@docs
is_binary
set_binary
unset_binary
BinaryRef
```

## Integrality utilities

```@docs
relax_integrality
fix_discrete_variables
```

## Extensions

```@docs
AbstractVariable
AbstractVariableRef
parse_one_operator_variable
```

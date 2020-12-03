# Variable API

More information can be found in the [Variables](@ref) section of the manual.

```@docs
@variable
@variables

name(::JuMP.VariableRef)
set_name(::JuMP.VariableRef, ::String)
variable_by_name

set_start_value
start_value

has_lower_bound
lower_bound
set_lower_bound
delete_lower_bound
LowerBoundRef

has_upper_bound
upper_bound
set_upper_bound
delete_upper_bound
UpperBoundRef

is_fixed
fix_value
fix
unfix
FixRef

is_integer
set_integer
unset_integer
IntegerRef

is_binary
set_binary
unset_binary
BinaryRef

relax_integrality

all_variables
owner_model
num_variables
VariableRef
index(::VariableRef)
optimizer_index(::VariableRef)

add_variable
check_belongs_to_model
build_variable
parse_one_operator_variable

AbstractVariable
AbstractVariableRef
VariableNotOwned
VariableConstrainedOnCreation
VariablesConstrainedOnCreation
```

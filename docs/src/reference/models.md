# Model API

More information can be found in the [Interacting with solvers](@ref) section of
the manual.

```@docs
Model

direct_model

mode

set_optimizer
optimizer_with_attributes

object_dictionary
unregister

set_silent
unset_silent

set_time_limit_sec
unset_time_limit_sec
time_limit_sec

solver_name

backend

Base.empty!(::Model)

get_optimizer_attribute
set_optimizer_attribute
set_optimizer_attributes

bridge_constraints
print_bridge_graph

write_to_file
Base.write(::IO, ::Model; ::MOI.FileFormats.FileFormat)

read_from_file
Base.read(::IO, ::Type{Model}; ::MOI.FileFormats.FileFormat)

ReferenceMap
copy_model
copy_extension_data
Base.copy(::AbstractModel)

operator_warn
error_if_direct_mode
```

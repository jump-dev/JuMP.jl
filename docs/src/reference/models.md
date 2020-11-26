# Model API

More information can be found in the [Interacting with solvers](@ref) section of
the manual.

```@docs
Model()
Model(::Any)

direct_model

set_optimizer
optimizer_with_attributes

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

write_to_file
Base.write(::IO, ::Model; ::MOI.FileFormats.FileFormat)

read_from_file
Base.read(::IO, ::Type{Model}; ::MOI.FileFormats.FileFormat)
```

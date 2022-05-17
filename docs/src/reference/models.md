# [Models](@id ModelAPI)

More information can be found in the [Models](@ref jump_models) section of
the manual.

## Constructors

```@docs
Model
direct_model
```

## Enums

```@docs
ModelMode
AUTOMATIC
MANUAL
DIRECT
```

## Basic functions

```@docs
backend
unsafe_backend
name(::AbstractModel)
solver_name
Base.empty!(::Model)
Base.isempty(::Model)
mode
object_dictionary
unregister
latex_formulation
set_string_names
```

## Working with attributes

```@docs
set_optimizer
optimizer_with_attributes
get_optimizer_attribute
set_optimizer_attribute
set_optimizer_attributes
set_silent
unset_silent
set_time_limit_sec
unset_time_limit_sec
time_limit_sec
```

## Copying

```@docs
ReferenceMap
copy_model
copy_extension_data
Base.copy(::AbstractModel)
```
## I/O

```@docs
write_to_file
Base.write(::IO, ::Model; ::MOI.FileFormats.FileFormat)
read_from_file
Base.read(::IO, ::Type{Model}; ::MOI.FileFormats.FileFormat)
```

## Caching Optimizer

```@docs
MOIU.reset_optimizer(::JuMP.Model)
MOIU.drop_optimizer(::JuMP.Model)
MOIU.attach_optimizer(::JuMP.Model)
```

## Bridge tools

```@docs
bridge_constraints
print_bridge_graph
```

## Extension tools

```@docs
AbstractModel
operator_warn
error_if_direct_mode
```

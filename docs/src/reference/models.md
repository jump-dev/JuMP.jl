# [Models](@id ModelAPI)

More information can be found in the [Models](@ref jump_models) section of
the manual.

## Constructors

```@docs
GenericModel
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
Base.empty!(::GenericModel)
Base.isempty(::GenericModel)
mode
object_dictionary
unregister
latex_formulation
set_string_names_on_creation
value_type
```

## Working with attributes

```@docs
set_optimizer
optimizer_with_attributes
get_attribute
set_attribute
set_attributes
set_silent
unset_silent
set_time_limit_sec
unset_time_limit_sec
time_limit_sec
set_start_values
get_optimizer_attribute
set_optimizer_attribute
set_optimizer_attributes
```

## Copying

```@docs
GenericReferenceMap
copy_model
copy_extension_data
Base.copy(::AbstractModel)
```
## I/O

```@docs
write_to_file
Base.write(::IO, ::GenericModel; ::MOI.FileFormats.FileFormat)
read_from_file
Base.read(::IO, ::Type{GenericModel{T}}; ::MOI.FileFormats.FileFormat) where {T}
```

## Caching Optimizer

```@docs
MOIU.reset_optimizer(::JuMP.GenericModel)
MOIU.drop_optimizer(::JuMP.GenericModel)
MOIU.attach_optimizer(::JuMP.GenericModel)
```

## Bridge tools

```@docs
add_bridge
remove_bridge
bridge_constraints
print_active_bridges
print_bridge_graph
```

## Extension tools

```@docs
AbstractModel
operator_warn
error_if_direct_mode
set_optimize_hook
```

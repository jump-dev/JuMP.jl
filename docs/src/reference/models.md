# [Models](@id ModelAPI)

```@docs
AUTOMATIC
MANUAL
DIRECT
Base.empty!(::Model)
Base.isempty(::Model)
Base.copy(::AbstractModel)
Base.write(::IO, ::Model; ::MOI.FileFormats.FileFormat)
Base.read(::IO, ::Type{Model}; ::MOI.FileFormats.FileFormat)
MOIU.reset_optimizer(::JuMP.Model)
MOIU.drop_optimizer(::JuMP.Model)
MOIU.attach_optimizer(::JuMP.Model)
```


Interacting with solvers
========================


```@docs
with_optimizer
```

The solvers can be set with
```@docs
JuMP.setoptimizer
```

```@docs
Model
```

## Direct mode

For advanced users, the model can be created without using a caching nor a
bridge optimizer using the [`JuMP.direct_model`](@ref) function:
```@docs
JuMP.direct_model
```

TODO: Describe the connection between JuMP and solvers. Automatic vs. Manual
mode. CachingOptimizer. How to set/change solvers. How to set parameters (solver
specific and generic). Status codes. Accessing the result.
How to accurately measure the solve time.

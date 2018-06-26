Containers
==========

Containers can be created using the `generatecontainer` function
```@docs
JuMP.generatecontainer
```

Containers in macro
-------------------

In the [`@variable`](@ref) (resp. [`@constraint`](@ref)) macro, containers of
variables (resp. constraints) can be created the following syntax

* `name[idxs1,idxs2,...,idxn]` creating an `n`-dimensional container of name
  `name`; or
* `[idxs1,idxs2,...,idxn]` creating an *anonymous* (see [Names](@ref))
  `n`-dimensional container.

Each expression `idxsi` can either be

* of the form `idxset` specifying that the `i`th index set of the container
  is `idxset`; or
* of the form `idxname=idxset` specifying that the `i`th index set of the
  container is `idxset` and allowing values used in the macro expression and
  keyword arguments to be expressions depending on the `idxname`.

The macro then creates the container using the [`JuMP.generatecontainer`](@ref) function
with the following arguments:

1. `VariableRef` for the [`@variable`](@ref) macro and `ConstraintRef` for the
   [`@constraint`](@ref) macro.
2. The index variables and arbitrary symbols for dimensions for which no
   variable index is specified.
3. The index sets specified.
4. The value of the `keyword` argument if given or `:Auto`.

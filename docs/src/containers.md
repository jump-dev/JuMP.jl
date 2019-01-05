```@meta
DocTestSetup = quote
    using JuMP
end
```

Containers
==========

JuMP provides a specialized container similar to
[`AxisArrays`](https://github.com/JuliaArrays/AxisArrays.jl) that enables
indexing with non-integer indices. Normally these are created automatically
by JuMP's macros. The following constructors can be used to create them
manually.

```@docs
JuMP.Containers.DenseAxisArray
```

Containers in macros
--------------------

The `generate_container` function encodes the logic for how containers are
constructed in JuMP's macros.
```@docs
JuMP.Containers.generate_container
```

In the [`@variable`](@ref) (resp. [`@constraint`](@ref)) macro, containers of
variables (resp. constraints) can be created the following syntax

* `name[index_set_1, index_set_2, ..., index_set_n]` creating an `n`-dimensional
  container of name `name`; or
* `[index_set_1, index_set_2, ..., index_set_n]` creating an *anonymous*
  `n`-dimensional container.

Each expression `index_set_i` can either be

* of the form `index_set` specifying that the `i`th index set of the container
  is `index_set`; or
* of the form `index_name=index_set` specifying that the `i`th index set of the
  container is `index_set` and allowing values used in the macro expression and
  keyword arguments to be expressions depending on the `index_name`.

The macro then creates the container using the
[`JuMP.Containers.generate_container`](@ref) function with the following
arguments:

1. `VariableRef` for the [`@variable`](@ref) macro and `ConstraintRef` for the
   [`@constraint`](@ref) macro.
2. The index variables and arbitrary symbols for dimensions for which no
   variable index is specified.
3. The index sets specified.
4. The value of the `keyword` argument if given or `:Auto`.

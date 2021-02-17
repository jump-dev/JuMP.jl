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
JuMP.Containers.SparseAxisArray
```

Containers in macros
--------------------

The `container` function encodes the logic for how containers are
constructed in JuMP's macros. The `@container` macro is available to create
containers independently of any JuMP model.
```@docs
JuMP.Containers.container
JuMP.Containers.default_container
JuMP.Containers.VectorizedProductIterator
JuMP.Containers.NestedIterator
JuMP.Containers.@container
```

In the [`@variable`](@ref) (resp. [`@constraint`](@ref)) macro, containers of
variables (resp. constraints) can be created with the following syntax:

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
[`JuMP.Containers.container`](@ref) function with the following
arguments:

1. A function taking as argument the value of the indices and returning the
   value to be stored in the container, e.g. a variable for the
   [`@variable`](@ref) macro and a constraint for the [`@constraint`](@ref)
   macro.
2. An iterator over the indices of the container.
4. The value of the `container` keyword argument if given.

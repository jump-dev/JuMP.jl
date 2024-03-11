# DimensionalData.jl

[DimensionalData.jl](https://github.com/rafaqz/DimensionalData.jl) provides
tools and abstractions for working with rectangular arrays that have named
dimensions.

!!! compat
    Using the DimensionalData extension with JuMP requires Julia v1.9 or later.

The DimensionalData extension in JuMP lets you construct a `DimensionalData.DimArray`
as an alternative to [`Containers.DenseAxisArray`](@ref) in the JuMP macros.

## License

DimensionalData.jl is licensed under the [MIT license](https://github.com/rafaqz/DimensionalData.jl/blob/main/LICENSE).

## Installation

Install DimensionalData using `Pkg.add`:

```julia
import Pkg
Pkg.add("DimensionalData")
```

## Use with JuMP

Activate the extension by loading both JuMP and DimensionalData:

```jldoctest ext_dimensional_data
julia> using JuMP, DimensionalData
```

Then, pass `container = DimensionalData.DimArray` in the [`@variable`](@ref),
[`@constraint`](@ref), or [`@expression`](@ref) macros:
```jldoctest ext_dimensional_data
julia> model = Model();

julia> @variable(
           model,
           x[i = 2:4, j = ["a", "b"]] >= i,
           container = DimensionalData.DimArray,
       )
╭─────────────────────────────╮
│ 3×2 DimArray{VariableRef,2} │
├─────────────────────────────┴─────────────────── dims ┐
  ↓ i Sampled{Int64} 2:4 ForwardOrdered Regular Points,
  → j Categorical{String} ["a", "b"] ForwardOrdered
└───────────────────────────────────────────────────────┘
 ↓ →  "a"     "b"
 2    x[2,a]  x[2,b]
 3    x[3,a]  x[3,b]
 4    x[4,a]  x[4,b]
```

Here `x` is a `DimensionalData.Dim` array object, so indexing uses the
DimensionalData syntax:
```jldoctest ext_dimensional_data
julia> x[At(2), At("a")]
x[2,a]

julia> x[2, 2]
x[3,b]
```

You can use `container = DimensionalData.DimArray` in the [`@expression`](@ref)
macro:
```jldoctest ext_dimensional_data
julia> @expression(
           model,
           expr[j = ["a", "b"]],
           sum(x[At(i), At(j)] for i in 2:4),
           container = DimensionalData.DimArray,
       )
╭───────────────────────────────╮
│ 2-element DimArray{AffExpr,1} │
├───────────────────────────────┴───────────── dims ┐
  ↓ j Categorical{String} ["a", "b"] ForwardOrdered
└───────────────────────────────────────────────────┘
 "a"  x[2,a] + x[3,a] + x[4,a]
 "b"  x[2,b] + x[3,b] + x[4,b]
```
and in [`@constraint`](@ref):
```jldoctest ext_dimensional_data
julia> @constraint(
           model,
           [j = ["a", "b"]],
           expr[At(j)] <= 1,
           container = DimensionalData.DimArray,
       )
╭──────────────────────────────────────────────────────────────────────────────╮
│ 2-element DimArray{ConstraintRef{Model, MathOptInterface.ConstraintIndex{MathOptInterface.ScalarAffineFunction{Float64}, MathOptInterface.LessThan{Float64}}, ScalarShape},1} │
├──────────────────────────────────────────────────────────────────────── dims ┤
  ↓ j Categorical{String} ["a", "b"] ForwardOrdered
└──────────────────────────────────────────────────────────────────────────────┘
 "a"  x[2,a] + x[3,a] + x[4,a] ≤ 1
 "b"  x[2,b] + x[3,b] + x[4,b] ≤ 1
```

## Documentation

See the [DimensionalData.jl documentation](https://rafaqz.github.io/DimensionalData.jl/stable/)
for more details on the syntax and features of `DimensionalData.DimArray`.

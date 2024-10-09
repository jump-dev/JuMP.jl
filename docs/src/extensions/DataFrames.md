# DataFrames.jl

[DataFrames.jl](https://github.com/JuliaData/DataFrames.jl) provides tools for
working with in-memory tabular data in Julia.

!!! compat
    Using the DataFrames extension with JuMP requires Julia v1.9 or later.

The DataFrames extension in JuMP lets you construct a `DataFrames.DataFrame` as
a container in the JuMP macros.

## License

DataFrames.jl is licensed under the [MIT license](https://github.com/JuliaData/DataFrames.jl/blob/main/LICENSE.md).

## Installation

Install DataFrames using `Pkg.add`:

```julia
import Pkg
Pkg.add("DataFrames")
```

## Use with JuMP

Activate the extension by loading both JuMP and DataFrames:

```jldoctest ext_data_frames
julia> using JuMP, DataFrames
```

Then, pass `container = DataFrames.DataFrame` in the [`@variable`](@ref),
[`@constraint`](@ref), or [`@expression`](@ref) macros:

```jldoctest ext_data_frames
julia> model = Model();

julia> @variable(
           model,
           x[i = 2:4, j = ["a", "b"]] >= i,
           container = DataFrames.DataFrame,
       )
6×3 DataFrame
 Row │ i      j       value
     │ Int64  String  GenericV…
─────┼──────────────────────────
   1 │     2  a       x[2,a]
   2 │     3  a       x[3,a]
   3 │     4  a       x[4,a]
   4 │     2  b       x[2,b]
   5 │     3  b       x[3,b]
   6 │     4  b       x[4,b]
```

Here `x` is a `DataFrames.DataFrame` array object, so operations use the
DataFrames syntax:

```jldoctest ext_data_frames
julia> x[x.j .== "a", [:i, :value]]
3×2 DataFrame
 Row │ i      value
     │ Int64  GenericV…
─────┼──────────────────
   1 │     2  x[2,a]
   2 │     3  x[3,a]
   3 │     4  x[4,a]

julia> DataFrames.unstack(x, :i, :j, :value)
3×3 DataFrame
 Row │ i      a           b
     │ Int64  GenericV…?  GenericV…?
─────┼───────────────────────────────
   1 │     2  x[2,a]      x[2,b]
   2 │     3  x[3,a]      x[3,b]
   3 │     4  x[4,a]      x[4,b]
```

You can use `container = DataFrames.DataFrame` in the [`@expression`](@ref)
macro:

```jldoctest ext_data_frames
julia> @expression(
           model,
           expr[j = ["a", "b"]],
           sum(x[x.j .== j, :value]),
           container = DataFrames.DataFrame,
       )
2×2 DataFrame
 Row │ j       value
     │ String  AffExpr
─────┼──────────────────────────────────
   1 │ a       x[2,a] + x[3,a] + x[4,a]
   2 │ b       x[2,b] + x[3,b] + x[4,b]
```

and in [`@constraint`](@ref):

```jldoctest ext_data_frames
julia> @constraint(
           model,
           [j = ["a", "b"]],
           sum(x[x.j .== j, :value]) <= 1,
           container = DataFrames.DataFrame,
       )
2×2 DataFrame
 Row │ j       value
     │ String  Constrai…
─────┼──────────────────────────────────────
   1 │ a       x[2,a] + x[3,a] + x[4,a] ≤ 1
   2 │ b       x[2,b] + x[3,b] + x[4,b] ≤ 1
```

### DataFrame-native syntax

While you can use indexing in JuMP's `@expression` and `@constraint` macros, it
may be more convenient to use DataFrames.jl split-apply-combine framework. For
example, `expr` can be equivalently written as:

```jldoctest ext_data_frames
julia> expr2 = model[:expr2] = DataFrames.combine(
           DataFrames.groupby(x, :j),
           :value => sum => :value,
       )
2×2 DataFrame
 Row │ j       value
     │ String  AffExpr
─────┼──────────────────────────────────
   1 │ a       x[2,a] + x[3,a] + x[4,a]
   2 │ b       x[2,b] + x[3,b] + x[4,b]
```

and the constraint could be written as

```jldoctest ext_data_frames
julia> df_constraint(v) = @constraint(model, sum(v) <= 1);

julia> DataFrames.combine(
           DataFrames.groupby(x, :j),
           :value => df_constraint => :value,
       )
2×2 DataFrame
 Row │ j       value
     │ String  Constrai…
─────┼──────────────────────────────────────
   1 │ a       x[2,a] + x[3,a] + x[4,a] ≤ 1
   2 │ b       x[2,b] + x[3,b] + x[4,b] ≤ 1
```

## Documentation

See the [DataFrames.jl documentation](https://dataframes.juliadata.org/stable/)
for more details on the syntax and features of `DataFrames.DataFrame`.

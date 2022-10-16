#  Copyright 2017, Iain Dunning, Joey Huchette, Miles Lubin, and contributors
#  This Source Code Form is subject to the terms of the Mozilla Public
#  License, v. 2.0. If a copy of the MPL was not distributed with this
#  file, You can obtain one at https://mozilla.org/MPL/2.0/.

_rows(x::Array) = zip(eachindex(x), Iterators.product(axes(x)...))

_rows(x::DenseAxisArray) = zip(vec(eachindex(x)), Iterators.product(axes(x)...))

_rows(x::SparseAxisArray) = zip(eachindex(x.data), keys(x.data))

"""
    table([f::Function=identity,] x, names::Symbol...)

Applies the function `f` to all elements of the variable container `x`,
returning the result as a `Vector` of `NamedTuple`s, where `names` are used for
the corresponding axis names. If `x` is an `N`-dimensional array, there must be
`N+1` names, so that the last name corresponds to the result of `f(x[i])`.

!!! info
    A `Vector` of `NamedTuple`s implements the [Tables.jl](https://github.com/JuliaData/Tables.jl)
    interface, and so the result can be used as input for any function
    that consumes a 'Tables.jl' compatible source.

## Example

```jldoctest; setup=:(using JuMP)
julia> model = Model();

julia> @variable(model, x[i=1:2, j=i:2] >= 0, start = i+j);

julia> Containers.table(start_value, x, :i, :j, :start)
3-element Vector{NamedTuple{(:i, :j, :start), Tuple{Int64, Int64, Float64}}}:
 (i = 1, j = 2, start = 3.0)
 (i = 1, j = 1, start = 2.0)
 (i = 2, j = 2, start = 4.0)
```
"""
function table(
    f::Function,
    x::Union{Array,DenseAxisArray,SparseAxisArray},
    names::Symbol...,
)
    got, want = length(names), ndims(x) + 1
    if got != want
        error("Invalid number column names provided: Got $got, expected $want.")
    end
    return [NamedTuple{names}((args..., f(x[i]))) for (i, args) in _rows(x)]
end

table(x, names::Symbol...) = table(identity, x, names...)

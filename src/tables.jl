#  Copyright 2017, Iain Dunning, Joey Huchette, Miles Lubin, and contributors
#  This Source Code Form is subject to the terms of the Mozilla Public
#  License, v. 2.0. If a copy of the MPL was not distributed with this
#  file, You can obtain one at https://mozilla.org/MPL/2.0/.

function _row_iterator(x::Array)
    return zip(eachindex(x), Iterators.product(axes(x)...))
end

function _row_iterator(x::Containers.DenseAxisArray)
    return zip(vec(eachindex(x)), Iterators.product(axes(x)...))
end

function _row_iterator(x::Containers.SparseAxisArray)
    return zip(eachindex(x.data), keys(x.data))
end

"""
    table([f::Function=identity,] x, value_name::Symbol, col_names::Symbol...)

Applies the function `f` to all elements of the variable container `x`,
returning the result as a `Vector` of `NamedTuple`s, where `col_names`
are used for the correspondig axis names, and `value_name` is used for the
result of `f(x[i])`.

!!! info
    A `Vector` of `NamedTuple`s implements the [Tables.jl](https://github.com/JuliaData/Tables.jl)
    interface, and so the result can be used as input for any function
    that consumes a 'Tables.jl' compatible source. 

## Example

```jldoctest; setup=:(using JuMP)
julia> model = Model();

julia> @variable(model, x[i=1:2, j=i:2] >= 0, start = i+j);

julia> table(start_value, x, :start, :I, :J)
3-element Vector{NamedTuple{(:I, :J, :start), Tuple{Int64, Int64, Float64}}}:
 (I = 1, J = 2, start = 3.0)
 (I = 1, J = 1, start = 2.0)
 (I = 2, J = 2, start = 4.0)
```
"""
function table(
    f::Function, 
    x::Union{Array,Containers.DenseAxisArray,Containers.SparseAxisArray}, 
    value_name::Symbol, 
    col_names::Symbol...,
)
    got, want = length(col_names), ndims(x)
    if got != want
        error("Invalid number column names provided: Got $got, expected $want.")
    end
    C = (col_names..., value_name)
    return [NamedTuple{C}((args..., f(x[i]))) for (i, args) in _row_iterator(x)]
end

function table(x, value_name::Symbol, col_names::Symbol...)
    return table(identity, x, value_name, col_names...)
end

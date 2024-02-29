#  Copyright 2017, Iain Dunning, Joey Huchette, Miles Lubin, and contributors
#  This Source Code Form is subject to the terms of the Mozilla Public
#  License, v. 2.0. If a copy of the MPL was not distributed with this
#  file, You can obtain one at https://mozilla.org/MPL/2.0/.

_rows(x::Array) = zip(eachindex(x), Iterators.product(axes(x)...))

_rows(x::DenseAxisArray) = zip(vec(eachindex(x)), Iterators.product(axes(x)...))

_rows(x::SparseAxisArray) = zip(eachindex(x.data), keys(x.data))

"""
    rowtable([f::Function=identity,] x; [header::Vector{Symbol} = Symbol[]])

Applies the function `f` to all elements of the variable container `x`,
returning the result as a `Vector` of `NamedTuple`s, where `header` is a vector
containing the corresponding axis names.

If `x` is an `N`-dimensional array, there must be `N+1` names, so that the last
name corresponds to the result of `f(x[i])`.

If `header` is left empty, then the default header is `[:x1, :x2, ..., :xN, :y]`.

!!! info
    A `Vector` of `NamedTuple`s implements the [Tables.jl](https://github.com/JuliaData/Tables.jl)
    interface, and so the result can be used as input for any function
    that consumes a 'Tables.jl' compatible source.

## Example

```jldoctest
julia> model = Model();

julia> @variable(model, x[i=1:2, j=i:2] >= 0, start = i+j);

julia> Containers.rowtable(start_value, x; header = [:i, :j, :start])
3-element Vector{@NamedTuple{i::Int64, j::Int64, start::Float64}}:
 (i = 1, j = 1, start = 2.0)
 (i = 1, j = 2, start = 3.0)
 (i = 2, j = 2, start = 4.0)

julia> Containers.rowtable(x)
3-element Vector{@NamedTuple{x1::Int64, x2::Int64, y::VariableRef}}:
 (x1 = 1, x2 = 1, y = x[1,1])
 (x1 = 1, x2 = 2, y = x[1,2])
 (x1 = 2, x2 = 2, y = x[2,2])
```
"""
function rowtable(
    f::Function,
    x::Union{Array,DenseAxisArray,SparseAxisArray};
    header::Vector{Symbol} = Symbol[],
)::Vector{<:NamedTuple}
    if isempty(header)
        header = Symbol[Symbol("x$i") for i in 1:ndims(x)]
        push!(header, :y)
    end
    got, want = length(header), ndims(x) + 1
    if got != want
        error(
            "Invalid number of column names provided: Got $got, expected $want.",
        )
    end
    elements = [(args..., f(x[i])) for (i, args) in _rows(x)]
    return NamedTuple{tuple(header...)}.(elements)
end

rowtable(x; kwargs...) = rowtable(identity, x; kwargs...)

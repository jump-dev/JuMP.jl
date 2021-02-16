#  Copyright 2017, Iain Dunning, Joey Huchette, Miles Lubin, and contributors
#  This Source Code Form is subject to the terms of the Mozilla Public
#  License, v. 2.0. If a copy of the MPL was not distributed with this
#  file, You can obtain one at https://mozilla.org/MPL/2.0/.

"""
    default_container(indices)

If `indices` is a [`NestedIterator`](@ref), return a
[`SparseAxisArray`](@ref). Otherwise, `indices` should be
a `VectorizedProductIterator` and the function returns
`Array` if all iterators of the product are `Base.OneTo` and retunrs
[`DenseAxisArray`](@ref) otherwise.
"""
function default_container end

"""
    container(f::Function, indices, ::Type{C})

Create a container of type `C` with indices `indices` and values at given
indices given by `f`.

    container(f::Function, indices)

Create a container with indices `indices` and values at given indices given by
`f`. The type of container used is determined by [`default_container`](@ref).

## Examples

```@jldoctest; setup=:(using JuMP)
julia> Containers.container((i, j) -> i + j, Containers.vectorized_product(Base.OneTo(3), Base.OneTo(3)))
3×3 Array{Int64,2}:
 2  3  4
 3  4  5
 4  5  6

julia> Containers.container((i, j) -> i + j, Containers.vectorized_product(1:3, 1:3))
2-dimensional DenseAxisArray{Int64,2,...} with index sets:
    Dimension 1, 1:3
    Dimension 2, 1:3
And data, a 3×3 Array{Int64,2}:
 2  3  4
 3  4  5
 4  5  6

julia> Containers.container((i, j) -> i + j, Containers.vectorized_product(2:3, Base.OneTo(3)))
2-dimensional DenseAxisArray{Int64,2,...} with index sets:
    Dimension 1, 2:3
    Dimension 2, Base.OneTo(3)
And data, a 2×3 Array{Int64,2}:
 3  4  5
 4  5  6

julia> Containers.container((i, j) -> i + j, Containers.nested(() -> 1:3, i -> i:3, condition = (i, j) -> isodd(i) || isodd(j)))
SparseAxisArray{Int64,2,Tuple{Int64,Int64}} with 5 entries:
  [1, 2]  =  3
  [2, 3]  =  5
  [3, 3]  =  6
  [1, 1]  =  2
  [1, 3]  =  4
```
"""
function container end

function container(f::Function, indices)
    return container(f, indices, default_container(indices))
end

const ArrayIndices{N} = VectorizedProductIterator{NTuple{N,Base.OneTo{Int}}}
default_container(::ArrayIndices) = Array
function container(f::Function, indices::ArrayIndices, ::Type{Array})
    return map(I -> f(I...), indices)
end

function _oneto(indices::AbstractVector{<:Integer})
    if indices == 1:length(indices)
        return Base.OneTo(length(indices))
    end
    return error("Index set for array is not one-based interval.")
end

function _oneto(::Any)
    return error("Index set for array is not one-based interval.")
end

_oneto(indices::Base.OneTo) = indices

function container(
    f::Function,
    indices::VectorizedProductIterator,
    ::Type{Array},
)
    return container(
        f,
        vectorized_product(_oneto.(indices.prod.iterators)...),
        Array,
    )
end
default_container(::VectorizedProductIterator) = DenseAxisArray
function container(
    f::Function,
    indices::VectorizedProductIterator,
    ::Type{DenseAxisArray},
)
    return DenseAxisArray(map(I -> f(I...), indices), indices.prod.iterators...)
end
default_container(::NestedIterator) = SparseAxisArray
# Returns the element type. If it is unknown but it is known to be `N`-tuples,
# returns `NTuple{N, Any}`.
_eltype_or_any(indices::Array) = eltype(indices)
function container(f::Function, indices, ::Type{SparseAxisArray})
    # Same as `map` but does not allocate the resulting vector.
    mappings = Base.Generator(I -> I => f(I...), indices)
    # Same as `Dict(mapping)` but it will error if two indices are the same.
    data = NoDuplicateDict(mappings)
    dict = data.dict
    if isempty(dict)
        # If `dict` is empty, it was not able to determine the type
        # of the key hence the type of `dict` is `Dict{Any, Any}`.
        # This is an issue since `SparseAxisArray` needs the key
        # type to be a tuple.
        dict = Dict{_eltype_or_any(indices),Any}()
    end
    return SparseAxisArray(dict)
end

#  Copyright 2017, Iain Dunning, Joey Huchette, Miles Lubin, and contributors
#  This Source Code Form is subject to the terms of the Mozilla Public
#  License, v. 2.0. If a copy of the MPL was not distributed with this
#  file, You can obtain one at https://mozilla.org/MPL/2.0/.

"""
    default_container(indices)

If `indices` is a [`NestedIterator`](@ref), return a
[`SparseAxisArray`](@ref). Otherwise, `indices` should be
a `VectorizedProductIterator` and the function returns
`Array` if all iterators of the product are `Base.OneTo` and returns
[`DenseAxisArray`](@ref) otherwise.
"""
function default_container end

"""
    AutoContainerType

Pass `AutoContainerType` to [`container`](@ref) to let the container type be
chosen based on the type of the indices using [`default_container`](@ref).
"""
struct AutoContainerType end

"""
    container(f::Function, indices[[, ::Type{C} = AutoContainerType], names])

Create a container of type `C` with index names `names`, indices `indices` and
values at given indices given by `f`.

If the method with `names` is not specialized on `Type{C}`, it falls back to
calling  `container(f, indices, c)` for backwards compatibility with containers
not supporting index names.

## Example

```jldoctest
julia> Containers.container((i, j) -> i + j, Containers.vectorized_product(Base.OneTo(3), Base.OneTo(3)))
3×3 Matrix{Int64}:
 2  3  4
 3  4  5
 4  5  6

julia> Containers.container((i, j) -> i + j, Containers.vectorized_product(1:3, 1:3))
2-dimensional DenseAxisArray{Int64,2,...} with index sets:
    Dimension 1, 1:3
    Dimension 2, 1:3
And data, a 3×3 Matrix{Int64}:
 2  3  4
 3  4  5
 4  5  6

julia> Containers.container((i, j) -> i + j, Containers.vectorized_product(2:3, Base.OneTo(3)))
2-dimensional DenseAxisArray{Int64,2,...} with index sets:
    Dimension 1, 2:3
    Dimension 2, Base.OneTo(3)
And data, a 2×3 Matrix{Int64}:
 3  4  5
 4  5  6

julia> Containers.container((i, j) -> i + j, Containers.nested(() -> 1:3, i -> i:3, condition = (i, j) -> isodd(i) || isodd(j)))
SparseAxisArray{Int64, 2, Tuple{Int64, Int64}} with 5 entries:
  [1, 1]  =  2
  [1, 2]  =  3
  [1, 3]  =  4
  [2, 3]  =  5
  [3, 3]  =  6
```
"""
function container(f::Function, indices, D, names)
    return container(f, indices, D)
end

function container(f::Function, indices, ::Type{AutoContainerType}, names)
    return container(f, indices, default_container(indices), names)
end

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
    names::Union{Nothing,AbstractVector} = nothing,
)
    if names !== nothing
        names = ntuple(i -> Symbol(names[i]), length(names))
    end
    return DenseAxisArray(
        map(I -> f(I...), indices),
        indices.prod.iterators...;
        names = names,
    )
end

default_container(::NestedIterator) = SparseAxisArray
# Returns the element type. If it is unknown but it is known to be `N`-tuples,
# returns `NTuple{N, Any}`.
_eltype_or_any(indices::Array) = eltype(indices)
function container(
    f::Function,
    indices,
    ::Type{SparseAxisArray},
    names::Union{Nothing,AbstractVector} = nothing,
)
    if names !== nothing
        names = ntuple(i -> Symbol(names[i]), length(names))
    end
    # Same as `map` but does not allocate the resulting vector.
    mappings = Base.Generator(I -> I => f(I...), indices)
    # Same as `OrderedCollections.OrderedDict(mapping)`, but it will error if
    # two indices are the same.
    data = NoDuplicateDict(mappings)
    return _sparseaxisarray(data.dict, f, indices, names)
end

# The NoDuplicateDict was able to infer the element type.
function _sparseaxisarray(
    dict::OrderedCollections.OrderedDict,
    ::Any,
    ::Any,
    names,
)
    return SparseAxisArray(dict, names)
end

# @default_eltype succeeded and inferred a tuple of the appropriate size!
# Use `return_types` to get the value type of the dictionary.
function _container_dict(
    K::Type{<:NTuple{N,Any}},
    f::Function,
    ::Type{<:NTuple{N,Any}},
) where {N}
    ret = Base.return_types(f, K)
    V = length(ret) == 1 ? first(ret) : Any
    return OrderedCollections.OrderedDict{K,V}()
end

# @default_eltype bailed and returned Any. Use an NTuple of Any of the
# appropriate size instead.
function _container_dict(::Any, ::Any, K::Type{<:NTuple{N,Any}}) where {N}
    return OrderedCollections.OrderedDict{K,Any}()
end

# @default_eltype bailed and returned Union{}. Use an NTuple of Any of the
# appropriate size instead. We need this method to avoid an ambiguity with
# `::Type{<:NTuple{N,Any}}` and `::Any`.
function _container_dict(
    ::Type{Union{}},
    ::Function,
    K::Type{<:NTuple{N,Any}},
) where {N}
    return OrderedCollections.OrderedDict{K,Any}()
end

# Calling `@default_eltye` on `x` isn't sufficient, because the iterator may
# skip every element based on the condition. Call it on an identical nested
# iterator, but this time without the condition.
_default_eltype(x::NestedIterator) = Base.@default_eltype nested(x.iterators...)
_default_eltype(x) = Base.@default_eltype x

# The NoDuplicateDict was not able to infer the element type. To make a
# best-guess attempt, collect all of the keys excluding the conditional
# statement (these must be defined, because the conditional applies to the
# lowest-level of the index loops), then get the eltype of the result.
function _sparseaxisarray(
    dict::OrderedCollections.OrderedDict{Any,Any},
    f,
    indices,
    names,
)
    @assert isempty(dict)
    d = _container_dict(_default_eltype(indices), f, _eltype_or_any(indices))
    return SparseAxisArray(d, names)
end

# Don't use length-1 tuples if there is only one index!
_container_key(i::Tuple) = i
_container_key(i::Tuple{T}) where {T} = i[1]

function container(f::Function, indices, D::Type{<:AbstractDict})
    return D(_container_key(i) => f(i...) for i in indices)
end

function container(::Function, ::Any, D::Type)
    return error(
        "Unable to build a container with the provided type $(D). Implement " *
        "`Containers.container(::Function, indices, ::Type{$(D)})`.",
    )
end

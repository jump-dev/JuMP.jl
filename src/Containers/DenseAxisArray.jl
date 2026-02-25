#  Copyright 2017, Iain Dunning, Joey Huchette, Miles Lubin, and contributors
#  This Source Code Form is subject to the terms of the Mozilla Public
#  License, v. 2.0. If a copy of the MPL was not distributed with this
#  file, You can obtain one at https://mozilla.org/MPL/2.0/.

# DenseAxisArray is inspired by the AxisArrays package.
# DenseAxisArray can be replaced with AxisArray once integer indices are no
# longer a special case. See discussions at:
# https://github.com/JuliaArrays/AxisArrays.jl/issues/117
# https://github.com/JuliaArrays/AxisArrays.jl/issues/84

struct _AxisLookup{D}
    data::D
end

# Default fallbacks.
Base.getindex(::_AxisLookup, key) = throw(KeyError(key))
Base.getindex(::_AxisLookup, key::Colon) = key

struct DenseAxisArray{T,N,Ax,L<:NTuple{N,_AxisLookup}} <: AbstractArray{T,N}
    data::Array{T,N}
    axes::Ax
    lookup::L
    names::NTuple{N,Symbol}
end

function Base.Array{T,N}(x::DenseAxisArray) where {T,N}
    return convert(Array{T,N}, copy(x.data))
end

function Base.hash(d::DenseAxisArray, h::UInt)
    return hash(d.data, hash(d.axes, hash(d.lookup, h)))
end

# Any -> _AxisLookup{<:Dict}: The most generic type of axis is a dictionary
# which maps keys to their index. This works for any AbstractVector-type axis.

function build_lookup(ax)
    d = Dict{eltype(ax),Int}()
    cnt = 1
    for el in ax
        if haskey(d, el)
            error("Repeated index $el. Index sets must have unique elements.")
        end
        d[el] = cnt
        cnt += 1
    end
    return _AxisLookup(d)
end

Base.getindex(x::_AxisLookup{Dict{K,Int}}, key) where {K} = x.data[key]
Base.getindex(::_AxisLookup{Dict{K,Int}}, key::Colon) where {K} = key

function Base.getindex(
    x::_AxisLookup{Dict{K,Int}},
    keys::AbstractVector,
) where {K}
    return [x[key] for key in keys]
end

function Base.getindex(
    x::_AxisLookup{Dict{K,Int}},
    key::T,
) where {T<:AbstractVector,K>:T}
    # In this method, we can't tell if `keys` is an actual key, or a vector of
    # keys. One common cause happens when `K` is `Any`.
    if haskey(x.data, key)
        if all(haskey(x.data, k) for k in key)
            error(
                "ambiguous use of getindex with key $key. We cannot tell if " *
                "you meant to return the single element corresponding to the " *
                "key, or a slice for each element in the key.",
            )
        end
        return x.data[key]
    end
    return [x[k] for k in key]
end

function Base.get(x::_AxisLookup{Dict{K,Int}}, key, default) where {K}
    return get(x.data, key, default)
end

# Base.OneTo -> _AxisLookup{<:Base.OneTo}: This one is an easy optimization, and
# avoids the unnecessary Dict lookup.

build_lookup(ax::Base.OneTo) = _AxisLookup(ax)
function Base.getindex(ax::_AxisLookup{<:Base.OneTo}, k::Integer)
    if !(k in ax.data)
        throw(KeyError(k))
    end
    return k
end

function Base.getindex(
    x::_AxisLookup{<:Base.OneTo},
    keys::AbstractVector{<:Integer},
)
    return [x[key] for key in keys]
end

function Base.get(ax::_AxisLookup{<:Base.OneTo}, k::Integer, default)
    return k in ax.data ? k : default
end

# AbstractUnitRange{<:Integer} -> _AxisLookup{Tuple{T,T}}: A related
# optimization to Base.OneTo.

function build_lookup(ax::AbstractUnitRange{T}) where {T<:Integer}
    return _AxisLookup{Tuple{T,T}}((first(ax), length(ax)))
end

function Base.getindex(
    x::_AxisLookup{Tuple{T,T}},
    key::Integer,
) where {T<:Integer}
    if !isequal(key, convert(T, key))
        throw(KeyError(key))
    end
    i = key - x.data[1] + 1
    if !(1 <= i <= x.data[2])
        throw(KeyError(key))
    end
    return i
end

function Base.getindex(
    x::_AxisLookup{Tuple{T,T}},
    keys::AbstractVector{<:Integer},
) where {T<:Integer}
    return [x[key] for key in keys]
end

function Base.get(
    x::_AxisLookup{Tuple{T,T}},
    key::Integer,
    default,
) where {T<:Integer}
    if !isequal(key, convert(T, key))
        return default
    end
    i = key - x.data[1] + 1
    if !(1 <= i <= x.data[2])
        return default
    end
    return i
end

# Implement a special case: If the axis is a vector of pairs, also allow tuples
# as indices. This is needed due to the behavior pairs and tuples when iterating
# through dictionaries, that is, `x[(k, v) in d]` gets added as `x[k => v]`, even
# though it looks to the user like they were tuples.

function Base.getindex(
    x::_AxisLookup{Dict{Pair{A,B},Int}},
    key::Tuple{A,B},
) where {A,B}
    return x.data[key[1]=>key[2]]
end

function Base.getindex(
    x::_AxisLookup{Dict{Pair{A,B},Int}},
    keys::AbstractVector{<:Tuple{A,B}},
) where {A,B}
    return [x[key] for key in keys]
end

function Base.get(
    x::_AxisLookup{Dict{Pair{A,B},Int}},
    key::Tuple{A,B},
    default,
) where {A,B}
    return get(x.data, key[1] => key[2], default)
end

_abstract_vector(x::AbstractVector) = x

function _abstract_vector(x::AbstractVector{<:CartesianIndex})
    return error(
        "Unsupported index type `CartesianIndex` in axis: $x. Cartesian " *
        "indices are restricted for indexing into and iterating over " *
        "multidimensional arrays.",
    )
end

_abstract_vector(x) = _abstract_vector([a for a in x])

_abstract_vector(x::AbstractArray) = vec(x)

function _abstract_vector(x::Number)
    @warn(
        "Axis contains one element: $x. If intended, you can safely " *
        "ignore this warning. To explicitly pass the axis with one " *
        "element, pass `[$x]` instead of `$x`.",
    )
    return _abstract_vector([x])
end

"""
    DenseAxisArray(data::Array{T, N}, axes...) where {T, N}

Construct a JuMP array with the underlying data specified by the `data` array
and the given axes. Exactly `N` axes must be provided, and their lengths must
match `size(data)` in the corresponding dimensions.

## Example

```jldoctest
julia> array = Containers.DenseAxisArray([1 2; 3 4], [:a, :b], 2:3)
2-dimensional DenseAxisArray{Int64,2,...} with index sets:
    Dimension 1, [:a, :b]
    Dimension 2, 2:3
And data, a 2×2 Matrix{Int64}:
 1  2
 3  4

julia> array[:b, 3]
4
```
"""
function DenseAxisArray(
    data::Array{T,N},
    axes...;
    names::Union{Nothing,NTuple{N,Symbol}} = nothing,
) where {T,N}
    @assert length(axes) == N
    new_axes = _abstract_vector.(axes)  # Force all axes to be AbstractVector!
    names = something(names, ntuple(n -> Symbol("#$n"), N))
    return DenseAxisArray(data, new_axes, build_lookup.(new_axes), names)
end

# A converter for different array types.
function DenseAxisArray(data::AbstractArray, axes...; kwargs...)
    return DenseAxisArray(collect(data), axes...; kwargs...)
end

"""
    DenseAxisArray{T}(undef, axes...) where T

Construct an uninitialized DenseAxisArray with element-type `T` indexed over the
given axes.

## Example

```jldoctest
julia> array = Containers.DenseAxisArray{Float64}(undef, [:a, :b], 1:2);

julia> fill!(array, 1.0)
2-dimensional DenseAxisArray{Float64,2,...} with index sets:
    Dimension 1, [:a, :b]
    Dimension 2, 1:2
And data, a 2×2 Matrix{Float64}:
 1.0  1.0
 1.0  1.0

julia> array[:a, 2] = 5.0
5.0

julia> array[:a, 2]
5.0

julia> array
2-dimensional DenseAxisArray{Float64,2,...} with index sets:
    Dimension 1, [:a, :b]
    Dimension 2, 1:2
And data, a 2×2 Matrix{Float64}:
 1.0  5.0
 1.0  1.0
```
"""
function DenseAxisArray{T}(::UndefInitializer, args...; kwargs...) where {T}
    return construct_undef_array(T, args; kwargs...)
end

function construct_undef_array(
    ::Type{T},
    args::Tuple{Vararg{Any,N}};
    kwargs...,
) where {T,N}
    data = Array{T,N}(undef, length.(args)...)
    return DenseAxisArray(data, args...; kwargs...)
end

Base.isempty(A::DenseAxisArray) = isempty(A.data)

# We specify `Ax` for the type of `axes` to avoid conflict where `axes` has type
# `Tuple{Vararg{Int,N}}`.
function Base.similar(
    A::DenseAxisArray{T,N,Ax},
    ::Type{S},
    axes::Ax,
) where {T,N,Ax<:Tuple{<:AbstractVector},S}
    return construct_undef_array(S, axes)
end

# Avoid conflict with method defined in Julia Base when the axes of the
# `DenseAxisArray` are all `Base.OneTo`:
function Base.similar(
    ::DenseAxisArray{T,N,Ax},
    ::Type{S},
    axes::Ax,
) where {T,N,Ax<:Tuple{Base.OneTo,Vararg{Base.OneTo}},S}
    return construct_undef_array(S, axes)
end

# AbstractArray interface

Base.size(A::DenseAxisArray) = size(A.data)
function Base.LinearIndices(::DenseAxisArray)
    return error("DenseAxisArray does not support this operation.")
end
Base.axes(A::DenseAxisArray) = A.axes
Base.CartesianIndices(a::DenseAxisArray) = CartesianIndices(a.data)

############
# Indexing #
############

function _is_assigned(A::DenseAxisArray{T,N}, idx...) where {T,N}
    if length(idx) == N
        keys = zeros(Int, N)
        for (i, v) in enumerate(idx)
            key = get(A.lookup[i], v, nothing)
            key === nothing && return false
            keys[i] = key
        end
        return isassigned(A.data, keys...)
    end
    return false
end

Base.isassigned(A::DenseAxisArray, idx...) = _is_assigned(A, idx...)

# For ambiguity with DenseAxisArray and Integer keys
Base.isassigned(A::DenseAxisArray, idx::Integer...) = _is_assigned(A, idx...)

# Disallow indexing with a mix of integers and Cartesian indices
Base.isassigned(A::DenseAxisArray, i::Union{Integer,CartesianIndex}...) = false

function Base.isassigned(A::DenseAxisArray, i::CartesianIndex)
    return isassigned(A.data, i)
end

Base.eachindex(A::DenseAxisArray) = eachindex(IndexStyle(A), A)

function Base.eachindex(::IndexCartesian, A::DenseAxisArray)
    return CartesianIndices(size(A.data))
end

function Base.eachindex(
    ::IndexCartesian,
    A::DenseAxisArray,
    B::DenseAxisArray...,
)
    ret = eachindex(A)
    for b in B
        if eachindex(b) != ret
            err = DimensionMismatch(
                "incompatible dimensions in eachindex. Got $(eachindex.((A, B...)))",
            )
            throw(err)
        end
    end
    return ret
end

# Use recursion over tuples to ensure the return-type of functions like
# `Base.to_index` are type-stable.
_getindex_recurse(::NTuple{0}, ::Tuple, ::Function) = ()
function _getindex_recurse(data::Tuple, keys::Tuple, condition::Function)
    d, d_rest = first(data), Base.tail(data)
    k, k_rest = first(keys), Base.tail(keys)
    remainder = _getindex_recurse(d_rest, k_rest, condition)
    return condition(k) ? tuple(d[k], remainder...) : remainder
end

function Base.to_index(A::DenseAxisArray{T,N}, idx) where {T,N}
    if length(idx) < N
        throw(BoundsError(A, idx))
    elseif any(i -> !isone(idx[i]), (N+1):length(idx))
        throw(KeyError(idx))
    end
    return _getindex_recurse(A.lookup, idx, x -> true)
end

_is_range(::Any) = false
_is_range(::Colon) = true
_is_range(::AbstractVector{<:Integer}) = true

function _kwargs_to_args(A::DenseAxisArray{T,N}; kwargs...) where {T,N}
    return ntuple(N) do i
        kw = keys(kwargs)[i]
        if A.names[i] != kw
            error(
                "Invalid index $kw in position $i. When using keyword " *
                "indexing, the indices must match the exact name and order " *
                "used when creating the container.",
            )
        end
        return kwargs[i]
    end
end

function Base.getindex(A::DenseAxisArray, args...; kwargs...)
    if !isempty(kwargs)
        if !isempty(args)
            error("Cannot index with mix of positional and keyword arguments")
        end
        return _getindex_inner(A, _kwargs_to_args(A; kwargs...)...)
    end
    return _getindex_inner(A, args...)
end

function _getindex_inner(A::DenseAxisArray{T}, args::Vararg{Any,N}) where {T,N}
    new_indices = Base.to_index(A, args)
    if !any(_is_range, new_indices)
        return A.data[new_indices...]::T
    end
    new_axes = _getindex_recurse(A.axes, new_indices, _is_range)
    names = A.names[findall(_is_range, new_indices)]
    return DenseAxisArray(A.data[new_indices...], new_axes...; names = names)
end

Base.getindex(A::DenseAxisArray, idx::CartesianIndex) = A.data[idx]

function Base.setindex!(
    A::DenseAxisArray{T,N},
    v,
    args...;
    kwargs...,
) where {T,N}
    if !isempty(kwargs)
        if !isempty(args)
            error("Cannot index with mix of positional and keyword arguments")
        end
        return setindex!(A, v, _kwargs_to_args(A; kwargs...)...)
    end
    return A.data[Base.to_index(A, args)...] = v
end

function Base.setindex!(
    A::DenseAxisArray{T},
    v::T,
    idx::CartesianIndex,
) where {T}
    A.data[idx] = v
    return
end

function Base.setindex!(
    A::DenseAxisArray{T,N},
    value::DenseAxisArray{T,N},
    args...,
) where {T,N}
    for key in Base.product(args...)
        A[key...] = value[key...]
    end
    return A
end

########
# Keys #
########

"""
    DenseAxisArrayKey

Structure to hold a DenseAxisArray key when it is viewed as key-value collection.
"""
struct DenseAxisArrayKey{T<:Tuple}
    I::T
end
Base.getindex(k::DenseAxisArrayKey, args...) = getindex(k.I, args...)
Base.getindex(a::DenseAxisArray, k::DenseAxisArrayKey) = a[k.I...]

function Base.setindex!(
    A::DenseAxisArray{T},
    value::T,
    key::DenseAxisArrayKey,
) where {T}
    return setindex!(A, value, key.I...)
end

struct DenseAxisArrayKeys{T<:Tuple,S<:DenseAxisArrayKey,N} <: AbstractArray{S,N}
    product_iter::Base.Iterators.ProductIterator{T}
    function DenseAxisArrayKeys(a::DenseAxisArray{TT,N,Ax}) where {TT,N,Ax}
        product_iter = Base.Iterators.product(a.axes...)
        return new{Ax,DenseAxisArrayKey{eltype(product_iter)},N}(product_iter)
    end
end

Base.size(iter::DenseAxisArrayKeys) = size(iter.product_iter)

function Base.eltype(iter::DenseAxisArrayKeys)
    return DenseAxisArrayKey{eltype(iter.product_iter)}
end

function Base.iterate(iter::DenseAxisArrayKeys)
    next = iterate(iter.product_iter)
    return next === nothing ? nothing : (DenseAxisArrayKey(next[1]), next[2])
end

function Base.iterate(iter::DenseAxisArrayKeys, state)
    next = iterate(iter.product_iter, state)
    return next === nothing ? nothing : (DenseAxisArrayKey(next[1]), next[2])
end

function Base.keys(a::DenseAxisArray)
    return DenseAxisArrayKeys(a)
end
Base.getindex(a::DenseAxisArrayKeys, idx::CartesianIndex) = a[idx.I...]

function Base.getindex(
    a::DenseAxisArrayKeys{T,S,N},
    args::Vararg{Int,N},
) where {T,S,N}
    key = _getindex_recurse(a.product_iter.iterators, args, x -> true)
    return DenseAxisArrayKey(key)
end

function Base.IndexStyle(::Type{DenseAxisArrayKeys{T,N,Ax}}) where {T,N,Ax}
    return IndexCartesian()
end

################
# Broadcasting #
################

# This implementation follows the instructions at
# https://docs.julialang.org/en/v1/manual/interfaces/#man-interfaces-broadcasting
# for implementing broadcast.

function Base.BroadcastStyle(::Type{<:DenseAxisArray})
    return Broadcast.ArrayStyle{DenseAxisArray}()
end

function _broadcast_axes_check(x::NTuple{N}) where {N}
    axes = first(x)
    for i in 2:N
        if x[i][1] != axes[1]
            error(
                "Unable to broadcast over DenseAxisArrays with different axes.",
            )
        end
    end
    return axes
end

_broadcast_args(f, x::Tuple) = _broadcast_args(f, first(x), Base.tail(x))

_broadcast_args(f, ::Tuple{}) = ()

_broadcast_args(f::Val{:axes}, x::Any, tail) = _broadcast_args(f, tail)

function _broadcast_args(f::Val{:axes}, x::DenseAxisArray, tail)
    return ((x.axes, x.lookup), _broadcast_args(f, tail)...)
end

_broadcast_args(f::Val{:data}, x::Any, tail) = (x, _broadcast_args(f, tail)...)

function _broadcast_args(f::Val{:data}, x::DenseAxisArray, tail)
    return (x.data, _broadcast_args(f, tail)...)
end

_broadcast_args(f::Val{:names}, x::Any, tail) = _broadcast_args(f, tail)

function _broadcast_args(f::Val{:names}, x::DenseAxisArray, tail)
    return (x.names, _broadcast_args(f, tail)...)
end

function Base.Broadcast.broadcasted(
    ::Broadcast.ArrayStyle{DenseAxisArray},
    f,
    args...,
)
    axes, lookup = _broadcast_axes_check(_broadcast_args(Val(:axes), args))
    new_args = _broadcast_args(Val(:data), args)
    names = _broadcast_args(Val(:names), args)
    return DenseAxisArray(broadcast(f, new_args...), axes, lookup, first(names))
end

########
# Show #
########

# Adapted printing from Julia's show.jl

# Copyright (c) 2009-2016: Jeff Bezanson, Stefan Karpinski, Viral B. Shah,
# and other contributors:
#
# https://github.com/JuliaLang/julia/contributors
#
# Permission is hereby granted, free of charge, to any person obtaining
# a copy of this software and associated documentation files (the
# "Software"), to deal in the Software without restriction, including
# without limitation the rights to use, copy, modify, merge, publish,
# distribute, sublicense, and/or sell copies of the Software, and to
# permit persons to whom the Software is furnished to do so, subject to
# the following conditions:
#
# The above copyright notice and this permission notice shall be
# included in all copies or substantial portions of the Software.

function Base.summary(io::IO, A::DenseAxisArray)
    _summary(io, A)
    for (k, ax) in enumerate(A.axes)
        print(io, "    Dimension $k, ")
        show(IOContext(io, :limit => true), ax)
        println(io)
    end
    return print(io, "And data, a ", summary(A.data))
end

function _summary(io, ::DenseAxisArray{T,N}) where {T,N}
    return println(
        io,
        "$N-dimensional DenseAxisArray{$T,$N,...} with index sets:",
    )
end

function Base.print_array(io::IO, X::DenseAxisArray{T,1}) where {T}
    return Base.print_matrix(io, X.data)
end

function Base.print_array(io::IO, X::DenseAxisArray{T,2}) where {T}
    return Base.print_matrix(io, X.data)
end

function Base.show_nd(
    io::IO,
    a::DenseAxisArray,
    print_matrix::Function,
    label_slices::Bool,
)
    limit::Bool = get(io, :limit, false)
    if isempty(a)
        return
    end
    tailinds = Base.tail(Base.tail(axes(a.data)))
    nd = ndims(a) - 2
    for I in CartesianIndices(tailinds)
        idxs = I.I
        if limit
            for i in 1:nd
                ii = idxs[i]
                ind = tailinds[i]
                if length(ind) > 10
                    if ii == ind[4] &&
                       all(d -> idxs[d] == first(tailinds[d]), 1:(i-1))
                        for j in (i+1):nd
                            szj = size(a.data, j + 2)
                            indj = tailinds[j]
                            if szj > 10 &&
                               first(indj) + 2 < idxs[j] <= last(indj) - 3
                                @goto skip
                            end
                        end
                        print(io, "...\n\n")
                        @goto skip
                    end
                    if ind[3] < ii <= ind[end-3]
                        @goto skip
                    end
                end
            end
        end
        if label_slices
            print(io, "[:, :, ")
            for i in 1:(nd-1)
                show(io, a.axes[i+2][idxs[i]])
                print(io, ", ")
            end
            show(io, a.axes[end][idxs[end]])
            println(io, "] =")
        end
        slice = view(a.data, axes(a.data, 1), axes(a.data, 2), idxs...)
        print_matrix(io, slice)
        print(io, idxs == map(last, tailinds) ? "" : "\n\n")
        @label skip
    end
end

function Base.show(io::IO, array::DenseAxisArray)
    summary(io, array)
    isempty(array) && return
    println(io, ":")
    return Base.print_array(io, array)
end

# TODO(odow): deprecate this at some point? We have to implement it here because
# it used to work in Julia 1.5. In Julia 1.6, the Base implementation changed to
# assume `x` was 1-indexed. It doesn't make sense to repeat a DenseAxisArray,
# but some users may depend on it's functionality so we have a work-around
# instead of just breaking code.
Base.repeat(x::DenseAxisArray; kwargs...) = repeat(x.data; kwargs...)

###
### view
###

_get_subaxis(::Colon, b::AbstractVector) = b

function _get_subaxis(a::AbstractVector, b::AbstractVector)
    for ai in a
        if !(ai in b)
            throw(KeyError(ai))
        end
    end
    return a
end

function _get_subaxis(a::T, b::AbstractVector{T}) where {T}
    if !(a in b)
        throw(KeyError(a))
    end
    return a
end

struct DenseAxisArrayView{T,N,D,A} <: AbstractArray{T,N}
    data::D
    axes::A
    function DenseAxisArrayView(x::DenseAxisArray{T}, args...) where {T}
        axis = _get_subaxis.(args, axes(x))
        N = length(_type_stable_axes(axis))
        return new{T,N,typeof(x),typeof(axis)}(x, axis)
    end
end

Base.view(A::DenseAxisArray, args...) = DenseAxisArrayView(A, args...)

Base.size(x::DenseAxisArrayView) = length.(axes(x))

_type_stable_axes(x::Tuple) = _type_stable_axes(first(x), Base.tail(x))
_type_stable_axes(::Tuple{}) = ()
_type_stable_axes(::Any, tail) = _type_stable_axes(tail)
function _type_stable_axes(x::AbstractVector, tail)
    return (x, _type_stable_axes(tail)...)
end

Base.axes(x::DenseAxisArrayView) = _type_stable_axes(x.axes)

_is_subaxis(key::K, axis::AbstractVector{K}) where {K} = key in axis

function _is_subaxis(key::AbstractVector{K}, axis::AbstractVector{K}) where {K}
    return all(k -> k in axis, key)
end

function _type_stable_args(axis::AbstractVector, ::Colon, axes, args)
    return (axis, _type_stable_args(axes, args)...)
end

function _type_stable_args(axis::AbstractVector, arg, axes, args)
    if !_is_subaxis(arg, axis)
        throw(KeyError(arg))
    end
    return (arg, _type_stable_args(axes, args)...)
end

function _type_stable_args(axis::Any, arg, axes, args)
    return (axis, _type_stable_args(axes, tuple(arg, args...))...)
end

function _type_stable_args(axes::Tuple, args::Tuple)
    return _type_stable_args(
        first(axes),
        first(args),
        Base.tail(axes),
        Base.tail(args),
    )
end

_type_stable_args(axes::Tuple, ::Tuple{}) = axes

function _fixed_indices(view_axes::Tuple, axes::Tuple)
    return filter(ntuple(i -> i, length(view_axes))) do i
        return !(typeof(view_axes[i]) <: eltype(axes[i]))
    end
end

function _kwargs_to_args(A::DenseAxisArrayView{T,N}; kwargs...) where {T,N}
    non_default_indices = _fixed_indices(A.axes, A.data.axes)
    return ntuple(N) do i
        kw = keys(kwargs)[i]
        if A.data.names[non_default_indices[i]] != kw
            error(
                "Invalid index $kw in position $i. When using keyword " *
                "indexing, the indices must match the exact name and order " *
                "used when creating the container.",
            )
        end
        return kwargs[i]
    end
end

function Base.getindex(x::DenseAxisArrayView, args...; kwargs...)
    if !isempty(kwargs)
        if !isempty(args)
            error("Cannot index with mix of positional and keyword arguments")
        end
        return getindex(x, _kwargs_to_args(x; kwargs...)...)
    end
    indices = _type_stable_args(x.axes, args)
    return getindex(x.data, indices...)
end

Base.getindex(a::DenseAxisArrayView, k::DenseAxisArrayKey) = a[k.I...]

function Base.setindex!(
    a::DenseAxisArrayView{T},
    value::T,
    k::DenseAxisArrayKey,
) where {T}
    return setindex!(a, value, k.I...)
end

function Base.setindex!(
    x::DenseAxisArrayView{T},
    value::T,
    args...;
    kwargs...,
) where {T}
    if !isempty(kwargs)
        if !isempty(args)
            error("Cannot index with mix of positional and keyword arguments")
        end
        return setindex!(x, value, _kwargs_to_args(x; kwargs...)...)
    end
    indices = _type_stable_args(x.axes, args)
    return setindex!(x.data, value, indices...)
end

function Base.eachindex(A::DenseAxisArrayView)
    # Return a generator so that we lazily evaluate the product instead of
    # collecting into a vector.
    #
    # In future, we might want to return the appropriate matrix of
    # `CartesianIndex` to avoid having to do the lookups with
    # `DenseAxisArrayKey`.
    return (DenseAxisArrayKey(k) for k in Base.product(axes(A)...))
end

Base.show(io::IO, x::DenseAxisArrayView) = print(io, x.data)

Base.print_array(io::IO, x::DenseAxisArrayView) = show(io, x)

function Base.summary(io::IO, x::DenseAxisArrayView)
    return print(io, "view(::DenseAxisArray, ", join(x.axes, ", "), "), over")
end

struct _InitNotProvided end

function Base.sum(
    f::F,
    x::Union{DenseAxisArray{T},DenseAxisArrayView{T}};
    dims = Colon(),
    init = _InitNotProvided(),
) where {F<:Function,T}
    if dims != Colon()
        return error(
            "`sum(x::DenseAxisArray; dims)` is not supported. Convert the array " *
            "to an `Array` using `sum(Array(x); dims=$dims)`, or use an explicit " *
            "for-loop summation instead.",
        )
    end
    if init == _InitNotProvided()
        return sum(f(xi) for xi in x)
    else
        return sum(f(xi) for xi in x; init)
    end
end

function Base.sum(x::Union{DenseAxisArray,DenseAxisArrayView}; kwargs...)
    return sum(identity, x; kwargs...)
end

function Base.promote_shape(a::DenseAxisArray, b::DenseAxisArray)
    if axes(a) != axes(b)
        msg = "Dimension and axes of a DenseAxisArray must match"
        throw(DimensionMismatch(msg))
    end
    return axes(a)
end

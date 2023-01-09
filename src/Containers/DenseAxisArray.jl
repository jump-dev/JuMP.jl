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
Base.:(==)(x::_AxisLookup{D}, y::_AxisLookup{D}) where {D} = x.data == y.data

# Default fallbacks.
Base.getindex(::_AxisLookup, key) = throw(KeyError(key))
Base.getindex(::_AxisLookup, key::Colon) = key

struct DenseAxisArray{T,N,Ax,L<:NTuple{N,_AxisLookup}} <: AbstractArray{T,N}
    data::Array{T,N}
    axes::Ax
    lookup::L
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
    key::K,
) where {K<:AbstractVector}
    return x.data[key]
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

function build_lookup(ax::AbstractUnitRange{<:Integer})
    return _AxisLookup((first(ax), length(ax)))
end

function Base.getindex(
    x::_AxisLookup{Tuple{T,T}},
    key::Integer,
) where {T<:Integer}
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
    i = key - x.data[1] + 1
    if !(1 <= i <= x.data[2])
        return default
    end
    return i
end

# Implement a special case: If the axis is a vector of pairs, also allow tuples
# as indices. This is needed due to the behavior pairs and tuples when iterating
# through dictionaries, i.e., `x[(k, v) in d]` gets added as `x[k => v]`, even
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

# Example
```jldoctest; setup=:(using JuMP)
julia> array = JuMP.Containers.DenseAxisArray([1 2; 3 4], [:a, :b], 2:3)
2-dimensional DenseAxisArray{Int64,2,...} with index sets:
    Dimension 1, Symbol[:a, :b]
    Dimension 2, 2:3
And data, a 2×2 Array{Int64,2}:
 1  2
 3  4

julia> array[:b, 3]
4
```
"""
function DenseAxisArray(data::Array{T,N}, axs...) where {T,N}
    @assert length(axs) == N
    new_axes = _abstract_vector.(axs)  # Force all axes to be AbstractVector!
    return DenseAxisArray(data, new_axes, build_lookup.(new_axes))
end

# A converter for different array types.
function DenseAxisArray(data::AbstractArray, axes...)
    return DenseAxisArray(collect(data), axes...)
end

"""
    DenseAxisArray{T}(undef, axes...) where T

Construct an uninitialized DenseAxisArray with element-type `T` indexed over the
given axes.

# Example
```jldoctest; setup=:(using JuMP)
julia> array = JuMP.Containers.DenseAxisArray{Float64}(undef, [:a, :b], 1:2);

julia> fill!(array, 1.0)
2-dimensional DenseAxisArray{Float64,2,...} with index sets:
    Dimension 1, Symbol[:a, :b]
    Dimension 2, 1:2
And data, a 2×2 Array{Float64,2}:
 1.0  1.0
 1.0  1.0

julia> array[:a, 2] = 5.0
5.0

julia> array[:a, 2]
5.0

julia> array
2-dimensional DenseAxisArray{Float64,2,...} with index sets:
    Dimension 1, Symbol[:a, :b]
    Dimension 2, 1:2
And data, a 2×2 Array{Float64,2}:
 1.0  5.0
 1.0  1.0
```
"""
function DenseAxisArray{T}(::UndefInitializer, axs...) where {T}
    return construct_undef_array(T, axs)
end

function construct_undef_array(::Type{T}, axs::Tuple{Vararg{Any,N}}) where {T,N}
    return DenseAxisArray(Array{T,N}(undef, length.(axs)...), axs...)
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
function Base.LinearIndices(A::DenseAxisArray)
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

if VERSION >= v"1.9.0-DEV"
    # Disallow indexing with a mix of integers and Cartesian indices
    function Base.isassigned(
        A::DenseAxisArray,
        i::Union{Integer,CartesianIndex}...,
    )
        return false
    end

    function Base.isassigned(A::DenseAxisArray, i::CartesianIndex)
        return isassigned(A.data, i)
    end
end

Base.eachindex(A::DenseAxisArray) = CartesianIndices(size(A.data))

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
_is_range(::Union{Vector{Int},Colon}) = true

function Base.getindex(A::DenseAxisArray{T,N}, idx...) where {T,N}
    new_indices = Base.to_index(A, idx)
    if !any(_is_range, new_indices)
        return A.data[new_indices...]::T
    end
    new_axes = _getindex_recurse(A.axes, new_indices, _is_range)
    return DenseAxisArray(A.data[new_indices...], new_axes...)
end

Base.getindex(A::DenseAxisArray, idx::CartesianIndex) = A.data[idx]

function Base.setindex!(A::DenseAxisArray{T,N}, v, idx...) where {T,N}
    return A.data[Base.to_index(A, idx)...] = v
end

function Base.setindex!(
    A::DenseAxisArray{T},
    v::T,
    idx::CartesianIndex,
) where {T}
    A.data[idx] = v
    return
end

function Base.IndexStyle(::Type{DenseAxisArray{T,N,Ax}}) where {T,N,Ax}
    return IndexAnyCartesian()
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

_broadcast_axes(x::Tuple) = _broadcast_axes(first(x), Base.tail(x))
_broadcast_axes(::Tuple{}) = ()
_broadcast_axes(::Any, tail) = _broadcast_axes(tail)
function _broadcast_axes(x::DenseAxisArray, tail)
    return ((x.axes, x.lookup), _broadcast_axes(tail)...)
end

_broadcast_args(x::Tuple) = _broadcast_args(first(x), Base.tail(x))
_broadcast_args(::Tuple{}) = ()
_broadcast_args(x::Any, tail) = (x, _broadcast_args(tail)...)
_broadcast_args(x::DenseAxisArray, tail) = (x.data, _broadcast_args(tail)...)

function Base.Broadcast.broadcasted(
    ::Broadcast.ArrayStyle{DenseAxisArray},
    f,
    args...,
)
    axes_lookup = _broadcast_axes_check(_broadcast_axes(args))
    new_args = _broadcast_args(args)
    return DenseAxisArray(broadcast(f, new_args...), axes_lookup...)
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
function _summary(io, A::DenseAxisArray{T,N}) where {T,N}
    return println(
        io,
        "$N-dimensional DenseAxisArray{$T,$N,...} with index sets:",
    )
end

if isdefined(Base, :print_array) # 0.7 and later
    function Base.print_array(io::IO, X::DenseAxisArray{T,1}) where {T}
        return Base.print_matrix(io, X.data)
    end
    function Base.print_array(io::IO, X::DenseAxisArray{T,2}) where {T}
        return Base.print_matrix(io, X.data)
    end
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
                       all(d -> idxs[d] == first(tailinds[d]), 1:i-1)
                        for j in i+1:nd
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

function Base.view(A::Containers.DenseAxisArray, args...)
    return DenseAxisArrayView(A, args...)
end

Base.size(x::DenseAxisArrayView) = length.(axes(x))

_type_stable_axes(x::Tuple) = _type_stable_axes(first(x), Base.tail(x))
_type_stable_axes(::Tuple{}) = ()
_type_stable_axes(::Any, tail) = _type_stable_axes(tail)
function _type_stable_axes(x::AbstractVector, tail)
    return (x, _type_stable_axes(tail)...)
end

Base.axes(x::DenseAxisArrayView) = _type_stable_axes(x.axes)

function _type_stable_args(axis::AbstractVector, ::Colon, axes, args)
    return (axis, _type_stable_args(axes, args)...)
end

function _type_stable_args(axis::AbstractVector, arg, axes, args)
    if !(arg in axis)
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

function Base.getindex(x::DenseAxisArrayView, args...)
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

function Base.setindex!(x::DenseAxisArrayView{T}, value::T, args...) where {T}
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

#  Copyright 2017, Iain Dunning, Joey Huchette, Miles Lubin, and contributors
#  This Source Code Form is subject to the terms of the Mozilla Public
#  License, v. 2.0. If a copy of the MPL was not distributed with this
#  file, You can obtain one at https://mozilla.org/MPL/2.0/.

"""
    struct SparseAxisArray{T,N,K<:NTuple{N, Any}} <: AbstractArray{T,N}
        data::Dict{K,T}
    end

`N`-dimensional array with elements of type `T` where only a subset of the
entries are defined. The entries with indices `idx = (i1, i2, ..., iN)` in
`keys(data)` has value `data[idx]`. Note that as opposed to
`SparseArrays.AbstractSparseArray`, the missing entries are not assumed to be
`zero(T)`, they are simply not part of the array. This means that the result of
`map(f, sa::SparseAxisArray)` or `f.(sa::SparseAxisArray)` has the same sparsity
structure than `sa` even if `f(zero(T))` is not zero.

## Examples

```jldoctest; setup=:(using JuMP)
julia> dict = Dict((:a, 2) => 1.0, (:a, 3) => 2.0, (:b, 3) => 3.0)
Dict{Tuple{Symbol,Int64},Float64} with 3 entries:
  (:b, 3) => 3.0
  (:a, 2) => 1.0
  (:a, 3) => 2.0

julia> array = JuMP.Containers.SparseAxisArray(dict)
JuMP.Containers.SparseAxisArray{Float64,2,Tuple{Symbol,Int64}} with 3 entries:
  [b, 3]  =  3.0
  [a, 2]  =  1.0
  [a, 3]  =  2.0

julia> array[:b, 3]
3.0
```
"""
struct SparseAxisArray{T,N,K<:NTuple{N,Any}} <: AbstractArray{T,N}
    data::Dict{K,T}
end

Base.length(sa::SparseAxisArray) = length(sa.data)

Base.IteratorSize(::Type{<:SparseAxisArray}) = Base.HasLength()

Base.iterate(sa::SparseAxisArray, args...) = iterate(values(sa.data), args...)

Base.hash(s::SparseAxisArray, h::UInt) = hash(s.data, h)

function Base.size(::SparseAxisArray)
    return error(
        "`Base.size` is not implemented for `SparseAxisArray` because " *
        "although it is a subtype of `AbstractArray`, it is conceptually " *
        "closer to a dictionary with `N`-dimensional keys. If you encounter " *
        "this error and you didn't call `size` explicitly, it is because " *
        "you called a method that is unsupported for `SparseAxisArray`s. " *
        "Consult the JuMP documentation for a list of supported operations.",
    )
end

# A `length` argument can be given because `IteratorSize` is `HasLength`
function Base.similar(
    ::SparseAxisArray{S,N,K},
    ::Type{T},
    length::Integer = 0,
) where {S,T,N,K}
    d = Dict{K,T}()
    if !iszero(length)
        sizehint!(d, length)
    end
    return SparseAxisArray(d)
end

Base.mapreduce(f, op, sa::SparseAxisArray) = mapreduce(f, op, values(sa.data))

Base.:(==)(sa1::SparseAxisArray, sa2::SparseAxisArray) = sa1.data == sa2.data

############
# Indexing #
############

Base.haskey(sa::SparseAxisArray, idx) = haskey(sa.data, idx)

function Base.haskey(sa::SparseAxisArray{T,1,Tuple{I}}, idx::I) where {T,I}
    return haskey(sa.data, (idx,))
end

function Base.setindex!(
    d::SparseAxisArray{T,N,K},
    value,
    idx::K,
) where {T,N,K<:NTuple{N,Any}}
    return setindex!(d, value, idx...)
end

function Base.setindex!(d::SparseAxisArray{T,N,K}, value, idx...) where {T,N,K}
    if length(idx) < N
        throw(BoundsError(d, idx))
    elseif _sliced_key_type(K, idx...) !== nothing
        throw(
            ArgumentError(
                "Slicing is not supported when calling `setindex!` on a" *
                " SparseAxisArray",
            ),
        )
    end
    return setindex!(d.data, value, idx)
end

function Base.getindex(
    d::SparseAxisArray{T,N,K},
    idx::K,
) where {T,N,K<:NTuple{N,Any}}
    return getindex(d, idx...)
end

function Base.getindex(d::SparseAxisArray{T,N,K}, idx...) where {T,N,K}
    if length(idx) < N
        throw(BoundsError(d, idx))
    end
    K2 = _sliced_key_type(K, idx...)
    if K2 !== nothing
        new_data = Dict{K2,T}(
            _sliced_key(k, idx) => v for (k, v) in d.data if _filter(k, idx)
        )
        return SparseAxisArray(new_data)
    end
    return getindex(d.data, idx)
end

# Method to check whether an index is an attempt at a slice.
@generated function _sliced_key_type(::Type{K}, idx...) where {K<:Tuple}
    expr = Expr(:curly, :Tuple)
    for i in 1:length(idx)
        Ki = K.parameters[i]
        if idx[i] <: Colon || idx[i] <: AbstractVector{<:Ki}
            push!(expr.args, Ki)
        end
    end
    return length(expr.args) == 1 ? :(nothing) : expr
end

# Methods to check whether a key `k` is a valid subset of `idx`.
_filter(::Any, ::Colon) = true
_filter(ki::Any, i::Any) = ki == i
_filter(ki::K, i::AbstractVector{<:K}) where {K} = ki in i
_filter(::Tuple{}, ::Tuple{}) = true
function _filter(k::Tuple, idx::Tuple)
    return _filter(k[1], idx[1]) && _filter(Base.tail(k), Base.tail(idx))
end

# Methods to subset the key into a new key, dropping all singleton axes.
_sliced_key(k, ::Any) = (k,)
_sliced_key(::K, ::K) where {K} = ()
_sliced_key(::Tuple{}, ::Tuple{}) = ()
function _sliced_key(k::Tuple, idx::Tuple)
    return tuple(
        _sliced_key(k[1], idx[1])...,
        _sliced_key(Base.tail(k), Base.tail(idx))...,
    )
end

Base.eachindex(d::SparseAxisArray) = keys(d.data)

################
# Broadcasting #
################

struct BroadcastStyle{N,K} <: Broadcast.BroadcastStyle end

function Base.BroadcastStyle(::Type{<:SparseAxisArray{T,N,K}}) where {T,N,K}
    return BroadcastStyle{N,K}()
end

# Disallow mixing broadcasts.
function Base.BroadcastStyle(::BroadcastStyle, ::Base.BroadcastStyle)
    return throw(
        ArgumentError(
            "Cannot broadcast Containers.SparseAxisArray with" *
            " another array of different type",
        ),
    )
end

# Allow broadcasting over scalars.
function Base.BroadcastStyle(
    style::BroadcastStyle,
    ::Base.Broadcast.DefaultArrayStyle{0},
)
    return style
end

# The fallback uses `axes` but recommend in the docstring to create a custom
# method for custom style if needed.
function Base.Broadcast.instantiate(
    bc::Base.Broadcast.Broadcasted{<:BroadcastStyle},
)
    return bc
end

# We use a couple of lisp-y tricks in `copy`.
#
# * `_indices` searches for the first SparseAxisArray in the broadcast and
#   returns the `keys(x.data)` as an iterator. It also checks that any other
#   SparseAxisArray's have the same set of keys.
# * Given an `index` from `_indices`, `_get_arg` walks the arguments of
#   broadcast and returns the equivalent of
#   `tuple([arg[index...] for arg in bc.args]...)` in a type-stable way.
function Base.copy(
    bc::Base.Broadcast.Broadcasted{BroadcastStyle{N,K}},
) where {N,K}
    dict = Dict(index => _getindex(bc, index) for index in _indices(bc.args...))
    if isempty(dict) && dict isa Dict{Any,Any}
        # If `dict` is empty (e.g., because there are no indices), then
        # inference will produce a `Dict{Any,Any}`, and we won't have enough
        # type information to call SparseAxisArray(dict). As a work-around, we
        # explicitly construct the type of the resulting SparseAxisArray.
        # For more, see JuMP issue #2867.
        return SparseAxisArray{Any,N,K}(dict)
    end
    return SparseAxisArray(dict)
end

# Recursion end for `_check_same_indices`. No more arguments to search.
_check_same_indices(::Any) = nothing

# This argument is not a `SparseAxisArray`. Keep searching.
function _check_same_indices(indices, ::Any, args...)
    return _check_same_indices(indices, args...)
end

# This argument is a SparseAxisArray. Check it has the same keys as the first
# argument.
function _check_same_indices(indices, x::SparseAxisArray, args...)
    if length(x.data) != length(indices) ||
       any(i -> !haskey(x.data, i), indices)
        throw(
            ArgumentError(
                "Cannot broadcast Containers.SparseAxisArray with" *
                " different indices",
            ),
        )
    end
    return _check_same_indices(indices, args...)
end

# `_indices` returns the set of keys for the SparseAxisArray being broadcasted
# over. It might be in any position of the arguments, so we keep stripping away
# arguments until we find one.
_indices(::Any, args...) = _indices(args...)

function _indices(x::SparseAxisArray, args...)
    indices = keys(x.data)
    # Now that we have some keys, check the remaining arugments for
    # SparseAxisArrays, and if we find one, check it has the same set of keys.
    _check_same_indices(indices, args...)
    return indices
end

function _indices(
    bc::Base.Broadcast.Broadcasted{BroadcastStyle{N,K}},
    args...,
) where {N,K}
    return _indices(bc.args...)
end

"""
    _get_arg(args::Tuple, index::Tuple)

Return a tuple corresponding to `tuple([arg[index] for arg in args]...)` in a
type-stable way.
"""
function _get_arg(args::Tuple, index::Tuple)
    return (_getindex(first(args), index), _get_arg(Base.tail(args), index)...)
end
_get_arg(::Tuple{}, ::Tuple) = ()

# We need this `getindex` lookalike because some of the `x` may be scalars that
# we are broadcasting over!
_getindex(x::SparseAxisArray, index) = getindex(x, index...)
_getindex(x::Any, ::Any) = x
_getindex(x::Ref, ::Any) = x[]

function _getindex(
    bc::Base.Broadcast.Broadcasted{BroadcastStyle{N,K}},
    index,
) where {N,K}
    return bc.f(_get_arg(bc.args, index)...)
end

@static if VERSION >= v"1.3"
    # `broadcast_preserving_zero_d` calls `axes(A)` which calls `size(A)` which
    # is not defined. When at least one argument is a `SparseAxisArray`, we can
    # simply redirect `broadcast_preserving_zero_d` to `broadcast` since we know
    # the result won't be zero dimensional.

    # Called by `A * 2`
    function Base.Broadcast.broadcast_preserving_zero_d(
        f,
        A::SparseAxisArray,
        As...,
    )
        return broadcast(f, A, As...)
    end
    # Called by `2 * A`
    function Base.Broadcast.broadcast_preserving_zero_d(
        f,
        x,
        A::SparseAxisArray,
        As...,
    )
        return broadcast(f, x, A, As...)
    end
end

########
# Show #
########

# Inspired from Julia SparseArrays stdlib package
# `Base.summary` is also called from `showerror` on `BoundsError`.
function Base.summary(io::IO, sa::SparseAxisArray)
    num_entries = length(sa.data)
    return print(
        io,
        typeof(sa),
        " with ",
        num_entries,
        isone(num_entries) ? " entry" : " entries",
    )
end

function Base.show(io::IO, ::MIME"text/plain", sa::SparseAxisArray)
    summary(io, sa)
    if !isempty(sa.data)
        println(io, ":")
        show(io, sa)
    end
    return
end

Base.show(io::IO, x::SparseAxisArray) = show(convert(IOContext, io), x)

function Base.show(io::IOContext, x::SparseAxisArray)
    if isempty(x)
        return show(io, MIME("text/plain"), x)
    end
    limit = get(io, :limit, false)::Bool
    half_screen_rows = limit ? div(displaysize(io)[1] - 8, 2) : typemax(Int)
    if !haskey(io, :compact)
        io = IOContext(io, :compact => true)
    end
    key_strings = [
        (join(key, ", "), value) for
        (i, (key, value)) in enumerate(x.data) if
        i < half_screen_rows || i > length(x) - half_screen_rows
    ]
    sort!(key_strings; by = x -> x[1])
    pad = maximum(length(x[1]) for x in key_strings)
    for (i, (key, value)) in enumerate(key_strings)
        print(io, "  [", rpad(key, pad), "]  =  ", value)
        if i != length(key_strings)
            println(io)
            if i == half_screen_rows
                println(io, "   ", " "^pad, "   \u22ee")
            end
        end
    end
    return
end

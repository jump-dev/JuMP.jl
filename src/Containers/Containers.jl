#  Copyright 2017, Iain Dunning, Joey Huchette, Miles Lubin, and contributors
#  This Source Code Form is subject to the terms of the Mozilla Public
#  License, v. 2.0. If a copy of the MPL was not distributed with this
#  file, You can obtain one at http://mozilla.org/MPL/2.0/.
#############################################################################
# JuMP
# An algebraic modeling language for Julia
# See http://github.com/JuliaOpt/JuMP.jl
#############################################################################

"""
    Containers

Module defining containers `DenseArray` and `SparseArray` that behaves as
regular `AbstractArray` but with custom indexes that are not necessarily
integers.
"""
module Containers

"""
    struct SparseArray{T,N} <: AbstractArray{T,N}
        data::Dict{Tuple,T}
    end

`N`-dimensional array with elements of type `T` where only a subset of the
entries are defined. The entries with indices `idx = (i1, i2, ..., iN)` in
`keys(data)` has value `data[idx]`. Note that as opposed to
`SparseArrays.AbstractSparseArray`, the missing entries are not assumed to be
`zero(T)`, they are simply not part of the array. This means that the result of
`map(f, sa::SparseArray)` or `f.(sa::SparseArray)` has the same sparsity
structure than `sa` even if `f(zero(T))` is not zero.
"""
struct SparseArray{T,N} <: AbstractArray{T,N}
    data::Dict{Tuple,T}
end

Base.length(sa::SparseArray) = length(sa.data)
Base.IteratorSize(::Type{<:SparseArray}) = Base.HasLength()
# By default `IteratorSize` for `Generator{<:AbstractArray{T,N}}` is
# `HasShape{N}`
Base.IteratorSize(::Type{Base.Generator{<:SparseArray}}) = Base.HasLength()
# Needed in `collect_to_with_first!`
Base.eachindex(g::Base.Generator{<:SparseArray}) = eachindex(g.iter)
@static if VERSION < v"0.7-"
    Base.start(sa::SparseArray) = start(values(sa.data))
    Base.next(sa::SparseArray, state) = start(values(sa.data), state)
    Base.done(sa::SparseArray, state) = start(values(sa.data), state)
else
    Base.iterate(sa::SparseArray, args...) = iterate(values(sa.data), args...)
end

# A `length` argument can be given because `IteratorSize` is `HasLength`
function Base.similar(sa::SparseArray{S,N}, ::Type{T},
                      length::Integer=0) where {S, T, N}
    d = Dict{Tuple,T}()
    if !iszero(length)
        sizehint!(d, length)
    end
    return SparseArray{T,N}(d)
end
# The generic implementation uses `LinearIndices`
function Base.collect_to_with_first!(dest::SparseArray, first_value, iterator,
                                     state)
    indices = eachindex(iterator)
    dest[first(indices)] = first_value
    for index in Iterators.drop(indices, 1)
        @static if VERSION < v"0.7-"
            element, state = next(iterator, state)
        else
            element, state = iterate(iterator, state)
        end
        dest[index] = element
    end
    return dest
end

function Base.mapreduce(f, op, sa::SparseArray)
    mapreduce(f, op, values(sa.data))
end
Base.:(==)(sa1::SparseArray, sa2::SparseArray) = sa1.data == sa2.data

########
# Show #
########

# Inspired from Julia SparseArrays stdlib package
function Base.show(io::IO, ::MIME"text/plain", sa::SparseArray)
    num_entries = length(sa.data)
    print(io, typeof(sa), " with ", num_entries, " stored ",
              isone(num_entries) ? "entry" : "entries")
    if !iszero(num_entries)
        println(io, ":")
        show(io, sa)
    end
end
Base.show(io::IO, x::SparseArray) = show(convert(IOContext, io), x)
function Base.show(io::IOContext, x::SparseArray)
    # TODO: make this a one-line form
    if isempty(x)
        return show(io, MIME("text/plain"), x)
    end
    limit::Bool = get(io, :limit, false)
    half_screen_rows = limit ? div(displaysize(io)[1] - 8, 2) : typemax(Int)
    key_string(key::Tuple) = join(key, ", ")
    print_entry(i) = i < half_screen_rows || i > length(x) - half_screen_rows
    pad = maximum(Int[print_entry(i) ? length(key_string(key)) : 0 for (i, key) in enumerate(keys(x.data))])
    if !haskey(io, :compact)
        io = IOContext(io, :compact => true)
    end
    for (i, (key, value)) = enumerate(x.data)
        if print_entry(i)
            print(io, "  ", '[', rpad(key_string(key), pad), "]  =  ", value)
            if i != length(x)
                println(io)
            end
        elseif i == half_screen_rows
            println(io, "   ", " "^pad, "   \u22ee")
        end
    end
end

############
# Indexing #
############

# Error for sa[..., :, ...]
function _colon_error() end
function _colon_error(::Colon, args...)
    throw(ArgumentError("Indexing with `:` is not supported by" *
                        " Containers.SparseArray"))
end
_colon_error(arg, args...) = _colon_error(args...)
function Base.setindex!(d::SparseArray, value, idx...)
    _colon_error(idx...)
    setindex!(d.data, value, idx)
end
function Base.getindex(d::SparseArray, idx...)
    _colon_error(idx...)
    getindex(d.data, idx)
end
Base.eachindex(d::SparseArray) = keys(d.data)

# Need to define it as indices may be non-integers
Base.to_index(d::SparseArray, idx) = idx

# Arbitrary typed indices. Linear indexing not supported.
struct IndexAnyCartesian <: Base.IndexStyle end
Base.IndexStyle(::IndexAnyCartesian, ::IndexAnyCartesian) = IndexAnyCartesian()
Base.IndexStyle(::Type{<:SparseArray}) = IndexAnyCartesian()
# eachindex redirect to keys
Base.keys(::IndexAnyCartesian, d::SparseArray) = keys(d)

################
# Broadcasting #
################

# Need to define it as indices may be non-integers
Base.Broadcast.newindex(d::SparseArray, idx) = idx

struct BroadcastStyle{N} <: Broadcast.BroadcastStyle end
function Base.BroadcastStyle(::BroadcastStyle, ::Base.BroadcastStyle)
    throw(ArgumentError("Cannot broadcast Containers.SparseArray with another" *
                        " array of different type"))
end
# Scalars can be used with SparseArray in broadcast
Base.BroadcastStyle(::BroadcastStyle{N}, ::Base.Broadcast.DefaultArrayStyle{0}) where {N} = BroadcastStyle{N}()
Base.BroadcastStyle(::Type{<:SparseArray{T, N}}) where {T, N} = BroadcastStyle{N}()
function Base.similar(b::Base.Broadcast.Broadcasted{BroadcastStyle{N}}, ::Type{T}) where {T, N, Ax}
    SparseArray{T, N}(Dict{Tuple, T}())
end

# Check that all SparseArrays involved have the same indices. The other
# arguments are scalars
function check_same_eachindex(each_index) end
check_same_eachindex(each_index, not_sa, args...) = check_same_eachindex(eachindex, args...)
function check_same_eachindex(each_index, sa::SparseArray, args...)
    if Set(each_index) != Set(eachindex(sa))
        throw(ArgumentError("Cannot broadcast Containers.SparseArray with" *
                            " different indices"))
    end
    check_same_eachindex(eachindex, args...)
end
_eachindex(not_sa, args...) = _eachindex(args...)
function _eachindex(sa::SparseArray, args...)
    each_index = eachindex(sa)
    check_same_eachindex(each_index, args...)
    return each_index
end
# Need to define it as it falls back to `axes` by default
function Base.eachindex(bc::Base.Broadcast.Broadcasted{<:BroadcastStyle})
    return _eachindex(bc.args...)
end

# The fallback uses `axes` but recommend in the docstring to create a custom
# method for custom style if needed.
Base.Broadcast.instantiate(bc::Base.Broadcast.Broadcasted{<:BroadcastStyle}) = bc

# The generic method in `Base` is `getindex(::Broadcasted, ::Union{Integer, CartesianIndex})`
# which is not applicable here since the index is not integer
# TODO make a change in `Base` so that we don't have to call a function starting
# with an `_`.
function Base.getindex(bc::Base.Broadcast.Broadcasted{<:BroadcastStyle}, I)
    return Base.Broadcast._broadcast_getindex(bc, I)
end

# The generic implementation fall back to converting `bc` to
# `Broadcasted{Nothing}`. It is advised in `Base` to define a custom method for
# custom styles. The fallback for `Broadcasted{Nothing}` is not appropriate as
# indices are not integers for `SparseArray`.
function Base.copyto!(dest::SparseArray{T, N}, bc::Base.Broadcast.Broadcasted{BroadcastStyle{N}}) where {T, N}
    for key in eachindex(bc)
        dest[key] = bc[key]
    end
    return dest
end

end

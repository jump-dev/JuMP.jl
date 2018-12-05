#  Copyright 2017, Iain Dunning, Joey Huchette, Miles Lubin, and contributors
#  This Source Code Form is subject to the terms of the Mozilla Public
#  License, v. 2.0. If a copy of the MPL was not distributed with this
#  file, You can obtain one at http://mozilla.org/MPL/2.0/.
#############################################################################
# JuMP
# An algebraic modeling language for Julia
# See http://github.com/JuliaOpt/JuMP.jl
#############################################################################

module Containers

"""
    struct SparseArray{T,N} <: AbstractArray{T,N}
        data::Dict{Tuple,T}
    end

`N`-dimensional array with elements of type `T` where only a subset of the entries are defined.
The entries with indices `idx = (i1, i2, ..., iN)` in `keys(data)` has value `data[idx]`.
Note that as opposed to `SparseArrays.AbstractSparseArray`, the missing entries are not assumed
to be `zero(T)`. The indexing structure of the array is not rectangular and the missing entries
are simply not part of the array. This means that the result of `map(f, sa::SparseArray)` or
`f.(sa::SparseArray)` has the same sparsity structure than `sa` even if `f(zero(T))` is not zero.
"""
struct SparseArray{T,N} <: AbstractArray{T,N}
    data::Dict{Tuple,T}
end

Base.length(sa::SparseArray) = length(sa.data)
Base.iterate(sa::SparseArray, args...) = iterate(sa, args...)
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

Base.setindex!(d::SparseArray, value, idx) = setindex!(d.data, value, idx)
Base.getindex(d::SparseArray, idx) = getindex(d.data, idx)
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
# TODO throw an helpful error when trying to mix SparseArray with a non-scalar
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
        throw(ArgumentError("Cannot broadcast Containers.SparseArray with different indices"))
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
# Need to define it as `axes` is not defined
Base.Broadcast.instantiate(bc::Base.Broadcast.Broadcasted{<:BroadcastStyle}) = bc
# Need to define it as indices are not integers
function Base.copyto!(dest::SparseArray{T, N}, bc::Base.Broadcast.Broadcasted{BroadcastStyle{N}}) where {T, N}
    for key in eachindex(bc)
        dest[key] = Base.Broadcast._broadcast_getindex(bc, key)
    end
    return dest
end

end

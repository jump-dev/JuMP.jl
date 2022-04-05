#  Copyright 2017, Iain Dunning, Joey Huchette, Miles Lubin, and contributors
#  This Source Code Form is subject to the terms of the Mozilla Public
#  License, v. 2.0. If a copy of the MPL was not distributed with this
#  file, You can obtain one at https://mozilla.org/MPL/2.0/.

# Lightweight unsafe view for vectors. It seems that the only way to avoid
# triggering allocations is to have only bitstype fields, so we store a pointer.
struct _VectorView{T} <: DenseVector{T}
    offset::Int
    len::Int
    ptr::Ptr{T}
end

Base.getindex(v::_VectorView, idx::Integer) = unsafe_load(v.ptr, idx + v.offset)

function Base.setindex!(v::_VectorView, value, idx::Integer)
    unsafe_store!(v.ptr, value, idx + v.offset)
    return v
end

function Base.setindex!(v::_VectorView{T}, value::T, idx::Vector{Int}) where {T}
    for i in idx
        v[i] = value
    end
    return v
end

Base.length(v::_VectorView) = v.len

function Base.fill!(v::_VectorView{T}, value) where {T}
    val = convert(T, value)
    for i in 1:length(v)
        v[i] = val
    end
    return v
end

function _rmul!(v::_VectorView{T}, value::T) where {T<:Number}
    for i in 1:length(v)
        v[i] *= value
    end
    return
end

function _reinterpret_unsafe(::Type{T}, x::Vector{R}) where {T,R}
    # how many T's fit into x?
    @assert isbitstype(T) && isbitstype(R)
    len = length(x) * sizeof(R)
    p = reinterpret(Ptr{T}, pointer(x))
    return _VectorView(0, div(len, sizeof(T)), p)
end

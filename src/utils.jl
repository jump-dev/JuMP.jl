#  Copyright 2017, Iain Dunning, Joey Huchette, Miles Lubin, and contributors
#  This Source Code Form is subject to the terms of the Mozilla Public
#  License, v. 2.0. If a copy of the MPL was not distributed with this
#  file, You can obtain one at http://mozilla.org/MPL/2.0/.

mutable struct IndexedVector{T}
    elts::Vector{T}
    nzidx::Vector{Int}
    nnz::Int
    empty::BitArray{1}
end

IndexedVector(::Type{T},n::Integer) where {T} = IndexedVector(zeros(T,n),zeros(Int,n),0,trues(n))

function addelt!(v::IndexedVector{T},i::Integer,val::T) where T
    if val != zero(T)
        if v.empty[i]  # new index
            v.elts[i] = val
            v.nzidx[v.nnz += 1] = i
            v.empty[i] = false
        else
            v.elts[i] += val
        end
    end
    return nothing
end

function rmz!(v::IndexedVector{T}) where T
    i = 1
    while i <= v.nnz
        j = v.nzidx[i]
        if v.elts[j] == zero(T)
            v.empty[j] = true
            # If i == v.nnz then this has no effect but it would be inefficient to branch
            v.nzidx[i] = v.nzidx[v.nnz]
            v.nnz -= 1
        else
            i += 1
        end
    end
end

function Base.empty!(v::IndexedVector{T}) where T
    elts = v.elts
    nzidx = v.nzidx
    empty = v.empty
    for i in 1:v.nnz
        elts[nzidx[i]] = zero(T)
        empty[nzidx[i]] = true
    end
    v.nnz = 0
end

Base.length(v::IndexedVector) = length(v.elts)
function Base.resize!(v::IndexedVector, n::Integer)
    if n > length(v)
        @assert v.nnz == 0 # only resize empty vector
        resize!(v.elts, n)
        fill!(v.elts,0)
        resize!(v.nzidx, n)
        resize!(v.empty, n)
        fill!(v.empty,true)
    end
end

# lightweight unsafe view for vectors
# it seems that the only way to avoid triggering
# allocations is to have only bitstype fields, so
# we store a pointer.
struct VectorView{T} <: DenseVector{T}
    offset::Int
    len::Int
    ptr::Ptr{T}
end

Base.getindex(v::VectorView,idx::Integer) = unsafe_load(v.ptr, idx+v.offset)
Base.setindex!(v::VectorView,value,idx::Integer) = unsafe_store!(v.ptr,value,idx+v.offset)
function Base.setindex!(v::VectorView{T},value::T,idx::Vector{Int}) where T
    for i in idx
        v[i] = value
    end
end
Base.length(v::VectorView) = v.len
function Base.fill!(v::VectorView{T},value) where T
    val = convert(T,value)
    for i in 1:length(v)
        v[i] = val
    end
    nothing
end
function Compat.rmul!(v::VectorView{T},value::T) where T<:Number
    for i in 1:length(v)
        v[i] *= value
    end
    nothing
end
function reinterpret_unsafe(::Type{T},x::Vector{R}) where {T,R}
    # how many T's fit into x?
    @assert isbitstype(T) && isbitstype(R)
    len = length(x)*sizeof(R)
    p = reinterpret(Ptr{T},pointer(x))
    return VectorView(0,div(len,sizeof(T)),p)
end

# For non-zero index arrays on 0.5; see
# https://github.com/alsam/OffsetArrays.jl#special-note-for-julia-05.
_size(A::AbstractArray) = map(length, axes(A))
_size(A) = size(A)
_size(A::AbstractArray, d) = d <= ndims(A) ? _size(A)[d] : 1

one_indexed(A) = all(x -> isa(x, Base.OneTo), axes(A))

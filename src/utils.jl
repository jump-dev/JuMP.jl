#  Copyright 2015, Iain Dunning, Joey Huchette, Miles Lubin, and contributors
#  This Source Code Form is subject to the terms of the Mozilla Public
#  License, v. 2.0. If a copy of the MPL was not distributed with this
#  file, You can obtain one at http://mozilla.org/MPL/2.0/.

type IndexedVector{T}
    elts::Vector{T}
    nzidx::Vector{Int}
    nnz::Int
    empty::BitArray{1}
end

IndexedVector{T}(::Type{T},n::Integer) = IndexedVector(zeros(T,n),zeros(Int,n),0,trues(n))

function addelt!{T}(v::IndexedVector{T},i::Integer,val::T)
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

function Base.empty!{T}(v::IndexedVector{T})
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

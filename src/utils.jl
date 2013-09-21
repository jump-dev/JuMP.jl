type IndexedVector{T}
    elts::Vector{T}
    nzidx::Vector{Int}
    nnz::Int
end

IndexedVector(T::Type,n::Integer) = IndexedVector(zeros(T,n),zeros(Int,n),0)

function addelt{T}(v::IndexedVector{T},i::Integer,val::T)
    if v.elts[i] == zero(T) # new index
        v.elts[i] = val
        v.nzidx[v.nnz += 1] = i
    else
        v.elts[i] += val
        if v.elts[i] == zero(T)
            # set to tiny value.

            # if two instances of T can sum to zero(T),
            # type must also implement one().
            # if not, code will never be reached
            v.elts[i] = 1e-50*one(T)
        end
    end
end

import Base.empty!

function empty!(v::IndexedVector)
    elts = v.elts
    nzidx = v.nzidx
    for i in 1:v.nnz
        elts[nzidx[i]] = zero(T)
    end
    v.nnz = 0
end


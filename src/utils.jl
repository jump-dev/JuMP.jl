type IndexedVector{T}
    elts::Vector{T}
    nzidx::Vector{Int}
    nnz::Int
end

IndexedVector(T::Type,n::Integer) = IndexedVector(zeros(T,n),zeros(Int,n),0)

function addelt{T}(v::IndexedVector{T},i::Integer,val::T)
    if v.elts[i] == 0 # new index
        v.elts[i] = val
        v.nzidx[v.nnz += 1] = i
    else
        v.elts[i] += val
        if v.elts[i] == 0
            # set to tiny value
            v.elts[i] = 1e-50
        end
    end
end

import Base.empty!

function empty!(v::IndexedVector)
    for i in 1:v.nnz
        v.elts[v.nzidx[i]] = 0
    end
    v.nnz = 0
end


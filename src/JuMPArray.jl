type JuMPArray{T,N,R<:OrdinalRange} <: JuMPContainer
    innerArray::Array{T,N}
    name::Symbol
    indexsets::NTuple{N,R}
end

function Base.getindex{T,N}(d::JuMPArray{T,N,UnitRange{Int}},indices::Int...)
    length(indices) == N || error("Wrong number of indices for ",d.name,", expected ",length(d.indexsets))
    idx = Array(Int, N)
    for i in 1:N
        idx[i] = indices[i] - start(d.indexsets[i]) + 1
    end
    return d.innerArray[idx...]
end

function Base.getindex{T,N}(d::JuMPArray{T,N,StepRange{Int,Int}},indices::Int...)
    length(indices) == N || error("Wrong number of indices for ",d.name,", expected ",length(d.indexsets))
    idx = Array(Int, N)
    steps  = 0
    starts = 0
    for i in 1:N
        steps  = step(d.indexsets[i])
        starts = start(d.indexsets[i])
        idx[i] = convert(Int,(indices[i]-starts)/steps) + 1
    end
    return d.innerArray[idx...]
end

function Base.setindex!{T,N}(d::JuMPArray{T,N,UnitRange{Int}},val::T,indices::Int...)
    length(indices) == N || error("Wrong number of indices for ",d.name,", expected ",length(d.indexsets))
    idx = Array(Int, N)
    for i in 1:N
        idx[i] = indices[i] - start(d.indexsets[i]) + 1
    end
    d.innerArray[idx...] = val
end

function Base.setindex!{T,N}(d::JuMPArray{T,N,StepRange{Int,Int}},val::T,indices::Int...)
    length(indices) == N || error("Wrong number of indices for ",d.name,", expected ",length(d.indexsets))
    idx = Array(Int, N)
    steps  = 0
    starts = 0
    for i in 1:N
        steps  = step(d.indexsets[i])
        starts = start(d.indexsets[i])
        idx[i] = convert(Int,(indices[i]-starts)/steps) + 1
    end
    d.innerArray[idx...] = val
end

Base.map{T,N,R}(f::Function,d::JuMPArray{T,N,R}) = 
    JuMPArray{T,N,R}(map(f,d.innerArray), d.name, d.indexsets)

Base.eltype{T,N,R}(x::JuMPArray{T,N,R}) = T


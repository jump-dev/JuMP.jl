#  Copyright 2017, Iain Dunning, Joey Huchette, Miles Lubin, and contributors
#  This Source Code Form is subject to the terms of the Mozilla Public
#  License, v. 2.0. If a copy of the MPL was not distributed with this
#  file, You can obtain one at http://mozilla.org/MPL/2.0/.

using Base.Meta
# multivarate "dictionary" used for collections of variables/constraints

abstract type JuMPContainer{T,N} end

include("JuMPArray.jl")

type IndexPair
    idxvar
    idxset
end

type JuMPDict{T,N} <: JuMPContainer{T,N}
    tupledict::Dict{NTuple{N,Any},T}
    meta::Dict{Symbol,Any}

    (::Type{JuMPDict{T,N}}){T,N}() = new{T,N}(Dict{NTuple{N,Any},T}(), Dict{Symbol,Any}())
end

function JuMPDict{T,N}(d::Dict{NTuple{N,Any},T})
    tmp = JuMPDict{T,N}()
    tmp.tupledict = d
    tmp
end

function JuMPDict{T,N}(d::Dict{NTuple{N,Any},T}, meta::Dict{Symbol,Any})
    tmp = JuMPDict{T,N}()
    tmp.tupledict = d
    tmp.meta = meta
    tmp
end

type JuMPContainerData
    name
    indexsets
    indexexprs::Vector{IndexPair}
    condition
end

# Needed by getvaluewarn when called by _mapInner
getname(data::JuMPContainerData) = data.name

#JuMPDict{T,N}(name::AbstractString) =
#    JuMPDict{T,N}(Dict{NTuple{N},T}(), name)

Base.getindex(d::JuMPDict, t...) = d.tupledict[t]
Base.setindex!(d::JuMPDict, value, t...) = (d.tupledict[t] = value)

Base.map(f::Function, d::JuMPArray) =
    JuMPArray(map(f, d.innerArray), d.indexsets, d.lookup, copy(d.meta))

function Base.map{T,N}(f::Function, d::JuMPDict{T,N})
    ret = Base.return_types(f, Tuple{T})
    R = (length(ret) == 1 ? ret[1] : Any)
    x = JuMPDict{R,N}()
    for (k,v) in d.tupledict
        x.tupledict[k] = f(v)
    end
    x.meta = copy(d.meta)
    return x
end

Base.isempty(d::JuMPContainer) = isempty(_innercontainer(d))

coloncheck(args::Number...) = nothing

function coloncheck(args...)
    if any(t -> t == Colon(), args)
        error("Colons not allowed as keys in JuMP containers.")
    end
end

# generate and instantiate a type which is indexed by the given index sets
# the following types of index sets are allowed:
# 0:K -- range with compile-time starting index
# S -- general iterable set
function gendict(instancename,T,idxsets...)
    N = length(idxsets)
    truearray = true
    for idxset in idxsets
        s = isexpr(idxset,:escape) ? idxset.args[1] : idxset
        if !(isexpr(s,:(:)) && length(s.args) == 2 && s.args[1] == 1)
            truearray = false
            break
        end
    end
    sizes = Expr(:tuple, [:(length($rng)) for rng in idxsets]...)
    if truearray
        :($instancename = Array{$T}($sizes...))
    else
        indexsets = Expr(:tuple, idxsets...)
        :($instancename = JuMPArray(Array{$T}($sizes...), $indexsets))
    end
end

metadata(x::Union{JuMPArray,JuMPDict}) = x.meta
metadata{T<:JuMPContainer}(::T) = error("Type $T has no field meta. This field is used to store metadata such as the JuMP.Model at the key :model.")
pushmeta!(x::JuMPContainer, sym::Symbol, val) = (metadata(x)[sym] = val)
getmeta(x::JuMPContainer, sym::Symbol) = metadata(x)[sym]
hasmeta(x::JuMPContainer, sym::Symbol) = haskey(metadata(x), sym)

# TODO do we want this or should we use broadcast syntax?
# for accessor in (:getdual, :getlowerbound, :getupperbound, :getstart)
#     @eval $accessor(x::AbstractArray) = map($accessor,x)
#     @eval $accessor(x::JuMPContainer) = map($accessor,x)
# end


_similar(x::Array) = Array{Float64}(size(x))
_similar{T}(x::Dict{T}) = Dict{T,Float64}()

_innercontainer(x::JuMPArray) = x.innerArray
_innercontainer(x::JuMPDict)  = x.tupledict

JuMPContainer_from(x::JuMPDict,inner) = JuMPDict(inner)
JuMPContainer_from(x::JuMPArray,inner) = JuMPArray(inner, x.indexsets)

# delegate zero-argument functions
for f in (:(Base.ndims), :(Base.length), :(Base.abs))
    @eval $f(x::JuMPArray) = $f(x.innerArray)
end

Base.length(x::JuMPDict) = length(x.tupledict)

Base.ndims{T,N}(x::JuMPDict{T,N}) = N
Base.abs(x::JuMPDict) = map(abs, x)
# avoid dangerous behavior with "end" (#730)
Base.endof(x::JuMPArray) = error("endof() (and \"end\" syntax) not implemented for JuMPArray objects.")
Base.size(x::JuMPArray) = error(string("size (and \"end\" syntax) not implemented for JuMPArray objects.",
"Use JuMP.size if you want to access the dimensions."))
Base.size(x::JuMPArray,k) = error(string("size (and \"end\" syntax) not implemented for JuMPArray objects.",
" Use JuMP.size if you want to access the dimensions."))
size(x::JuMPArray) = size(x.innerArray)
size(x::JuMPArray,k) = size(x.innerArray,k)
# for uses of size() within JuMP
size(x) = Base.size(x)
size(x,k) = Base.size(x,k)
# delegate one-argument functions
Base.issymmetric(x::JuMPArray) = issymmetric(x.innerArray)

Base.eltype{T}(x::JuMPContainer{T}) = T

Base.full(x::JuMPContainer) = x

# keys/vals iterations for JuMPContainers
Base.keys(d::JuMPDict)    = keys(d.tupledict)
Base.values(d::JuMPDict)  = values(d.tupledict)

Base.keys(d::JuMPArray)   = KeyIterator(d)
Base.values(d::JuMPArray) = ValueIterator(d.innerArray)

# Wrapper type so that you can't access the values directly
type ValueIterator{T,N}
    x::Array{T,N}
end
Base.start(it::ValueIterator)   =  start(it.x)
Base.next(it::ValueIterator, k) =   next(it.x, k)
Base.done(it::ValueIterator, k) =   done(it.x, k)
Base.length(it::ValueIterator)  = length(it.x)

type KeyIterator{JA<:JuMPArray}
    x::JA
    dim::Int
    next_k_cache::Array{Any,1}
    function (::Type{KeyIterator{JA}}){JA}(d)
        n = ndims(d.innerArray)
        new{JA}(d, n, Array{Any}(n+1))
    end
end

KeyIterator{JA}(d::JA) = KeyIterator{JA}(d)

function indexability(x::JuMPArray)
    for i in  1:length(x.indexsets)
        if !method_exists(getindex, (typeof(x.indexsets[i]),))
            return false
        end
    end

    return true
end

function Base.start(it::KeyIterator)
    if indexability(it.x)
        return start(it.x.innerArray)
    else
        return notindexable_start(it.x)
    end
end

@generated function notindexable_start{T,N,NT}(x::JuMPArray{T,N,NT})
    quote
        $(Expr(:tuple, 0, [:(start(x.indexsets[$i])) for i in 1:N]...))
    end
end

@generated function _next{T,N,NT}(x::JuMPArray{T,N,NT}, k::Tuple)
    quote
        $(Expr(:tuple, [:(next(x.indexsets[$i], k[$i+1])[1]) for i in 1:N]...))
    end
end

function Base.next(it::KeyIterator, k::Tuple)
    cartesian_key = _next(it.x, k)
    pos = -1
    for i in 1:it.dim
        if !done(it.x.indexsets[i], next(it.x.indexsets[i], k[i+1])[2] )
            pos = i
            break
        end
    end
    if pos == - 1
        it.next_k_cache[1] = 1
        return cartesian_key, tuple(it.next_k_cache...)
    end
    it.next_k_cache[1] = 0
    for i in 1:it.dim
        if i < pos
            it.next_k_cache[i+1] = start(it.x.indexsets[i])
        elseif i == pos
            it.next_k_cache[i+1] = next(it.x.indexsets[i], k[i+1])[2]
        else
            it.next_k_cache[i+1] = k[i+1]
        end
    end
    cartesian_key, tuple(it.next_k_cache...)
end

Base.done(it::KeyIterator, k::Tuple) = (k[1] == 1)

@generated __next{T,N,NT}(x::JuMPArray{T,N,NT}, k::Integer) =
    quote
        subidx = ind2sub(size(x),k)
        $(Expr(:tuple, [:(x.indexsets[$i][subidx[$i]]) for i in 1:N]...)), next(x.innerArray,k)[2]
    end
Base.next(it::KeyIterator, k) = __next(it.x,k::Integer)
Base.done(it::KeyIterator, k) = done(it.x.innerArray, k::Integer)

Base.length(it::KeyIterator)  = length(it.x.innerArray)

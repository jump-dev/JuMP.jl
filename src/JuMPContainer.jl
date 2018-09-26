#  Copyright 2017, Iain Dunning, Joey Huchette, Miles Lubin, and contributors
#  This Source Code Form is subject to the terms of the Mozilla Public
#  License, v. 2.0. If a copy of the MPL was not distributed with this
#  file, You can obtain one at http://mozilla.org/MPL/2.0/.

using Base.Meta
# multivarate "dictionary" used for collections of variables/constraints

abstract type JuMPContainer{T,N} end

include("JuMPArray.jl")

mutable struct IndexPair
    idxvar
    idxset
end

mutable struct JuMPDict{T,N} <: JuMPContainer{T,N}
    tupledict::Dict{NTuple{N,Any},T}
    meta::Dict{Symbol,Any}

    JuMPDict{T,N}() where {T,N} = new{T,N}(Dict{NTuple{N,Any},T}(), Dict{Symbol,Any}())
end

function JuMPDict(d::Dict{NTuple{N,Any},T}) where {T,N}
    tmp = JuMPDict{T,N}()
    tmp.tupledict = d
    tmp
end

function JuMPDict(d::Dict{NTuple{N,Any},T}, meta::Dict{Symbol,Any}) where {T,N}
    tmp = JuMPDict{T,N}()
    tmp.tupledict = d
    tmp.meta = meta
    tmp
end

mutable struct JuMPContainerData
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

function Base.map(f::Function, d::JuMPDict{T,N}) where {T,N}
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
if VERSION < v"0.7-"
    function _gendict_checkgen(s::Expr)
        isexpr(s,:(:)) && length(s.args) == 2 && s.args[1] == 1
    end
else
    function _gendict_checkgen(s::Expr)
        isexpr(s,:call) && s.args[1] == :(:) && length(s.args) == 3 && s.args[2] == 1
    end
end
_gendict_checkgen(::Any) = false
function gendict(instancename,T,idxsets...)
    N = length(idxsets)
    truearray = true
    for idxset in idxsets
        s = isexpr(idxset,:escape) ? idxset.args[1] : idxset
        if !_gendict_checkgen(s)
            truearray = false
            break
        end
    end
    sizes = Expr(:tuple, [:(length($rng)) for rng in idxsets]...)
    if truearray
        :($instancename = Array{$T}(undef, $sizes...))
    else
        indexsets = Expr(:tuple, idxsets...)
        :($instancename = JuMPArray(Array{$T}(undef, $sizes...), $indexsets))
    end
end


metadata(x::Union{JuMPArray,JuMPDict}) = x.meta
metadata(::T) where {T<:JuMPContainer} = error("Type $T has no field meta. This field is used to store metadata such as the JuMP.Model at the key :model.")
pushmeta!(x::JuMPContainer, sym::Symbol, val) = (metadata(x)[sym] = val)
getmeta(x::JuMPContainer, sym::Symbol) = metadata(x)[sym]
hasmeta(x::JuMPContainer, sym::Symbol) = haskey(metadata(x), sym)

# duck typing approach -- if eltype(innerArray) doesn't support accessor, will fail
for accessor in (:getdual, :getlowerbound, :getupperbound, :getvalue)
    @eval $accessor(x::AbstractArray) = map($accessor,x)
end
# With JuMPContainer, we take care in _mapInner of the warning if NaN values are returned
# by the accessor so we use the inner accessor that does not generate warnings
for (accessor, inner) in ((:getdual, :_getDual), (:getlowerbound, :getlowerbound), (:getupperbound, :getupperbound), (:getvalue, :_getValue))
    @eval $accessor(x::JuMPContainer) = _map($inner,x)
end


_similar(x::Array) = Array{Float64}(undef, size(x))
_similar(x::Dict{T}) where {T} = Dict{T,Float64}()

_innercontainer(x::JuMPArray) = x.innerArray
_innercontainer(x::JuMPDict)  = x.tupledict

# Warning for getter returning NaN
function _warnnan(f, data)
    if f === _getValue
        getvaluewarn(data)
    elseif f === _getDual
        getdualwarn(data)
    end
end

function _mapInner(f, x::JuMPContainer{T}) where T
    vars = _innercontainer(x)
    vals = _similar(vars)
    warnedyet = false
    for I in eachindex(vars)
        tmp = f(vars[I])
        if isnan(tmp) && !warnedyet
            _warnnan(f, getname(x))
            warnedyet = true
        end
        vals[I] = tmp
    end
    vals
end

JuMPContainer_from(x::JuMPDict,inner) = JuMPDict(inner)
JuMPContainer_from(x::JuMPArray,inner) = JuMPArray(inner, x.indexsets)

# The name _map is used instead of map so that this function is called instead of map(::Function, ::JuMPArray)
function _map(f, x::JuMPContainer{T}) where T
    mapcontainer_warn(f, x)
    ret = JuMPContainer_from(x, _mapInner(f, x))
    # I guess copy!(::Dict, ::Dict) isn't defined, so...
    for (key,val) in metadata(x)
        pushmeta!(ret, key, val)
    end
    if T == Variable
        m = _getmodel(x)
        # cache indexing info for new container for printing purposes
        m.varData[ret] = printdata(x)
    end
    ret
end

# delegate zero-argument functions
for f in (:(Base.ndims), :(Base.length), :(Base.abs))
    @eval $f(x::JuMPArray) = $f(x.innerArray)
end

Base.length(x::JuMPDict) = length(x.tupledict)

Base.ndims(x::JuMPDict{T,N}) where {T,N} = N
Base.abs(x::JuMPDict) = map(abs, x)
# avoid dangerous behavior with "end" (#730)
Compat.lastindex(x::JuMPArray,d=1) = error("lastindex() (and \"end\" syntax) not implemented for JuMPArray objects.")
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
LinearAlgebra.issymmetric(x::JuMPArray) = issymmetric(x.innerArray)

Base.eltype(x::JuMPContainer{T}) where {T} = T

# keys/vals iterations for JuMPContainers
Base.keys(d::JuMPDict)    = keys(d.tupledict)
Base.values(d::JuMPDict)  = values(d.tupledict)

Base.keys(d::JuMPArray)   = KeyIterator(d)
Base.values(d::JuMPArray) = ValueIterator(d.innerArray)

# Wrapper type so that you can't access the values directly
mutable struct ValueIterator{T,N}
    x::Array{T,N}
end
if VERSION < v"0.7-"
    Base.start(it::ValueIterator)   =  start(it.x)
    Base.next(it::ValueIterator, k) =   next(it.x, k)
    Base.done(it::ValueIterator, k) =   done(it.x, k)
else
    Base.iterate(it::ValueIterator) = iterate(it.x)
    Base.iterate(it::ValueIterator, k) = iterate(it.x, k)
end
Base.length(it::ValueIterator) = length(it.x)

mutable struct KeyIterator{JA<:JuMPArray}
    x::JA
    dim::Int
    next_k_cache::Array{Any,1}
    function KeyIterator{JA}(d) where JA
        n = ndims(d.innerArray)
        new{JA}(d, n, Array{Any}(undef, VERSION < v"0.7-" ? n+1 : n))
    end
end

KeyIterator(d::JA) where {JA} = KeyIterator{JA}(d)

function indexability(x::JuMPArray)
    for i in  1:length(x.indexsets)
        if !hasmethod(getindex, (typeof(x.indexsets[i]),))
            return false
        end
    end

    return true
end

if VERSION < v"0.7-"
    function Base.start(it::KeyIterator)
        if indexability(it.x)
            return start(it.x.innerArray)
        else
            return notindexable_start(it.x)
        end
    end
else
    function Base.iterate(it::KeyIterator)
        if indexability(it.x)
            return iterate(it.x.innerArray)
        else
            return notindexable_start(it.x)
        end
    end
end

if VERSION < v"0.7-"
    @generated function notindexable_start(x::JuMPArray{T,N,NT}) where {T,N,NT}
        quote
            $(Expr(:tuple, 0, [:(start(x.indexsets[$i])) for i in 1:N]...))
        end
    end
else
    _add_zero(item_state::Tuple) = (item_state[1], (0, item_state[2]...))
    function notindexable_start(x::JuMPArray{T,N,NT}) where {T,N,NT}
        item_states = ntuple(i -> iterate(x.indexsets[i]), Val(N))
        map(item_state -> item_state[1], item_states), item_states
    end
end

if VERSION < v"0.7-"
    @generated function _next(x::JuMPArray{T,N,NT}, k::Tuple) where {T,N,NT}
        quote
            $(Expr(:tuple, [:(next(x.indexsets[$i], k[$i+1])[1]) for i in 1:N]...))
        end
    end
else
    function _next(x::JuMPArray{T,N,NT}, k::Tuple) where {T,N,NT}
        map(item_state -> item_state[1], k)
    end
end

if VERSION < v"0.7-"
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
else
    function Base.iterate(it::KeyIterator, k::Tuple)
        pos = -1
        for i in 1:it.dim
            if iterate(it.x.indexsets[i], k[i][2]) ≠ nothing
                pos = i
                break
            end
        end
        if pos == -1
            return nothing
        end
        for i in 1:it.dim
            if i < pos
                it.next_k_cache[i] = iterate(it.x.indexsets[i])
            elseif i == pos
                it.next_k_cache[i] = iterate(it.x.indexsets[i], k[i][2])
            else
                it.next_k_cache[i] = k[i]
            end
        end
        next_k = tuple(it.next_k_cache...)
        _next(it.x, next_k), next_k
    end
end

@generated __next(x::JuMPArray{T,N,NT}, k::Integer) where {T,N,NT} =
    quote
        subidx = _ind2sub(size(x), k)
        $(Expr(:tuple, [:(x.indexsets[$i][subidx[$i]]) for i in 1:N]...)), next(x.innerArray,k)[2]
    end

if VERSION < v"0.7-"
    Base.next(it::KeyIterator, k) = __next(it.x,k::Integer)
    Base.done(it::KeyIterator, k) = done(it.x.innerArray, k::Integer)
else
    function Base.iterate(it::KeyIterator, k)
        iterate(it.x.innerArray, k::Integer) ≡ nothing && (return nothing)
        __next(it.x, k::Integer)
    end
end

Base.length(it::KeyIterator)  = length(it.x.innerArray)

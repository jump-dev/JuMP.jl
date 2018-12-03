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
struct ValueIterator{T,N}
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

struct KeyIterator{T}
    inner::T
end
KeyIterator(d::JuMPArray) = KeyIterator(Iterators.product(d.indexsets...))

if VERSION < v"0.7-"
    Base.start(it::KeyIterator) = start(it.inner)
    Base.next(it::KeyIterator, state) = next(it.inner, state)
    Base.done(it::KeyIterator, state) = done(it.inner, state)
else
    Base.iterate(it::KeyIterator) = iterate(it.inner)
    Base.iterate(it::KeyIterator, state) = iterate(it.inner, state)
end
Base.length(it::KeyIterator)  = length(it.inner)

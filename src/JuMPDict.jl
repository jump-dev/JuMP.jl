#  Copyright 2015, Iain Dunning, Joey Huchette, Miles Lubin, and contributors
#  This Source Code Form is subject to the terms of the Mozilla Public
#  License, v. 2.0. If a copy of the MPL was not distributed with this
#  file, You can obtain one at http://mozilla.org/MPL/2.0/.

using Base.Meta
# multivarate "dictionary" used for collections of variables/constraints

abstract JuMPContainer{T,N}
abstract JuMPArray{T,N,Q} <: JuMPContainer{T,N} # Q is true if all index sets are of the form 1:n
typealias OneIndexedArray{T,N} JuMPArray{T,N,true}

type IndexPair
    idxvar
    idxset
end
#= Generated on the fly
type JuMPArray{T}
    innerArray::Array{T,N}
    name::String
    indexsets
end
=#
type JuMPDict{T,N} <: JuMPContainer{T,N}
    tupledict::Dict{NTuple{N,Any},T}
    name::Symbol
    indexsets
    indexexprs::Vector{IndexPair}
    condition
end


#JuMPDict{T,N}(name::String) =
#    JuMPDict{T,N}(Dict{NTuple{N},T}(), name)

Base.getindex(d::JuMPDict, t...) = d.tupledict[t]
Base.setindex!(d::JuMPDict, value, t...) = (d.tupledict[t] = value)

function Base.map{T,N}(f::Function, d::JuMPDict{T,N})
    ret = Base.return_types(f, @compat(Tuple{T}))
    R = (length(ret) == 1 ? ret[1] : Any)
    x = JuMPDict(Dict{NTuple{N,Any},R}(), d.name, copy(d.indexsets), copy(d.indexexprs), copy(d.condition))
    for (k,v) in d.tupledict
        x.tupledict[k] = f(v)
    end
    return x
end

Base.isempty(d::JuMPArray) = (isempty(d.innerArray))
Base.isempty(d::JuMPDict)  = (isempty(d.tupledict))

# generate and instantiate a type which is indexed by the given index sets
# the following types of index sets are allowed:
# 0:K -- range with compile-time starting index
# S -- general iterable set
macro gendict(instancename,T,idxpairs,idxsets...)
    N = length(idxsets)
    allranges = all(s -> (isexpr(s,:(:)) && length(s.args) == 2), idxsets)
    if allranges
        # JuMPArray
        typename = symbol(string("JuMPArray",gensym()))
        offset = Array(Int,N)
        dictnames = Array(Symbol,N)
        for i in 1:N
            if isa(idxsets[i].args[1],Int)
                offset[i] = 1 - idxsets[i].args[1]
            else
                error("Currently only ranges with integer compile-time starting values are allowed as index sets. $(idxsets[i].args[1]) is not an integer in range $(idxsets[i]).")
            end
        end
        truearray = all(idxsets) do rng
            isexpr(rng,:(:)) && # redundant, but might as well
            rng.args[1] == 1 && # starts from 1
            length(rng.args) == 2 || rng.args[2] == 1 # steps of one
        end
        typecode = :(type $(typename){T} <: JuMPArray{T,$N,$truearray}; innerArray::Array{T,$N}; name::String;
                            indexsets; indexexprs::Vector{IndexPair} end)
        getidxlhs = :(Base.getindex(d::$(typename)))
        setidxlhs = :(setindex!(d::$(typename),val))
        getidxrhs = :(Base.getindex(d.innerArray))
        setidxrhs = :(setindex!(d.innerArray,val))
        maplhs = :(Base.map(f::Function,d::$(typename)))
        if truearray # return a julia Array here
            maprhs = :(map(f,d.innerArray))
        else
            maprhs = :($(typename)(map(f,d.innerArray),d.name,d.indexsets,d.indexexprs))
        end
        wraplhs = :(JuMPContainer_from(d::$(typename),inner)) # helper function that wraps array into JuMPArray of similar type
        wraprhs = :($(typename)(inner, d.name, d.indexsets, d.indexexprs))
        for i in 1:N
            varname = symbol(string("x",i))

            push!(getidxlhs.args,:($varname))
            push!(setidxlhs.args,:($varname))

            push!(getidxrhs.args,:(isa($varname, Int) ? $varname+$(offset[i]) : $varname ))
            push!(setidxrhs.args,:($varname+$(offset[i])))

        end

        badgetidxlhs = :(Base.getindex(d::$(typename),wrong...))
        badgetidxrhs = :(error("Wrong number of indices for ",d.name,
                               ", expected ",length(d.indexsets)))

        funcs = :($getidxlhs = $getidxrhs; $setidxlhs = $setidxrhs;
                  $maplhs = $maprhs; $badgetidxlhs = $badgetidxrhs)
        if !truearray
            funcs = :($funcs; $wraplhs = $wraprhs)
        end
        geninstance = :($(esc(instancename)) = $(typename)(Array($T),$(string(instancename)),$(esc(Expr(:tuple,idxsets...))),$(idxpairs)))
        for i in 1:N
            push!(geninstance.args[2].args[2].args, :(length($(esc(idxsets[i])))))
        end
        eval(Expr(:toplevel, typecode))
        eval(Expr(:toplevel, funcs))

        return geninstance

        #= TODO: Use this with code in JuMPArray.jl once Julia can make it efficient
        if all([length(s.args) == 3 for s in idxsets])
            # has step
            for i in 1:N
                if length(idxsets[i].args) == 2
                    push!(idxsets[i].args, idxsets[i].args[2])
                    idxsets[i].args[2] = 1
                end
            end
        end
        geninstance = :($(esc(instancename)) = JuMPArray(Array($T),$(quot(instancename)),$(esc(Expr(:tuple,idxsets...)))))
        for i in 1:N
            push!(geninstance.args[2].args[2].args, :(length($(esc(idxsets[i])))))
        end
        return geninstance
        =#
    else
        # JuMPDict
        return :(
            $(esc(instancename)) = JuMPDict{$T,$N}(Dict{NTuple{$N},$T}(),$(quot(instancename)), $(esc(Expr(:tuple,idxsets...))), $idxpairs, :())
        )
    end
end

# duck typing approach -- if eltype(innerArray) doesn't support accessor, will fail
for accessor in (:getDual, :getLower, :getUpper)
    @eval $accessor(x::JuMPContainer) = map($accessor,x)
end

_similar(x::Array) = Array(Float64,size(x))
_similar{T}(x::Dict{T}) = Dict{T,Float64}()

_innercontainer(x::JuMPArray) = x.innerArray
_innercontainer(x::JuMPDict)  = x.tupledict

function _getValueInner(x)
    vars = _innercontainer(x)
    vals = _similar(vars)
    warnedyet = false
    for I in eachindex(vars)
        tmp = _getValue(vars[I])
        if isnan(tmp) && !warnedyet
            warn("Variable value not defined for entry of $(x.name). Check that the model was properly solved.")
            warnedyet = true
        end
        vals[I] = tmp
    end
    vals
end


JuMPContainer_from(x::JuMPDict,inner) =
    JuMPDict(inner, x.name, x.indexsets, x.indexexprs, x.condition)
JuMPContainer_from(x::OneIndexedArray, inner) = inner

function getValue(x::JuMPContainer)
    getvalue_warn(x)
    JuMPContainer_from(x,_getValueInner(x))
end

# delegate zero-argument functions
for f in (:(Base.endof), :(Base.ndims), :(Base.length), :(Base.abs), :(Base.start))
    @eval $f(x::JuMPArray) = $f(x.innerArray)
end

for f in (:(Base.first), :(Base.length), :(Base.start))
    @eval $f(x::JuMPDict)  = $f(x.tupledict)
end
Base.ndims{T,N}(x::JuMPDict{T,N}) = N
Base.abs(x::JuMPDict) = map(abs, x)
# delegate one-argument functions
Base.size(x::JuMPArray)   = size(x.innerArray)
Base.size(x::JuMPArray,k) = size(x.innerArray,k)
Base.issym(x::JuMPArray) = issym(x.innerArray)
Base.trace(x::OneIndexedArray) = trace(x.innerArray)
Base.diag(x::OneIndexedArray) = diag(x.innerArray)
Base.diagm{T}(x::JuMPArray{T,1,true}) = diagm(x.innerArray)

function _local_index(indexsets, dim, k)
    n = length(indexsets)
    cprod = 1
    for i in 1:(dim-1)
        cprod *= length(indexsets[i])
    end
    idx = Compat.ceil(Integer, mod1(k,cprod*length(indexsets[dim])) / cprod)
    return indexsets[dim][idx]
end

function Base.next(x::JuMPArray,k)
    (var, gidx) = next(x.innerArray, k)
    key = map(i->_local_index(x.indexsets, i, k), 1:length(x.indexsets))
    return (tuple(key..., var), gidx)
end
function Base.next(x::JuMPDict,k)
    ((idx,var),gidx) = next(x.tupledict,k)
    return (tuple(idx..., var), gidx)
end
Base.done(x::JuMPArray,k) = done(x.innerArray,k)
Base.done(x::JuMPDict,k)  = done(x.tupledict,k)

Base.eltype{T}(x::JuMPContainer{T}) = T

Base.full(x::JuMPContainer) = x
Base.full(x::OneIndexedArray) = x.innerArray

export @gendict

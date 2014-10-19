using Base.Meta
# multivarate "dictionary" used for collections of variables/constraints

abstract JuMPContainer{T}
abstract JuMPArray{T} <: JuMPContainer{T}

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
type JuMPDict{T,N} <: JuMPContainer{T}
    tupledict::Dict{NTuple{N},T}
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
    ret = Base.return_types(f, (T,))
    R = (length(ret) == 1 ? ret[1] : Any)
    x = JuMPDict(Dict{NTuple{N},R}(), d.name, copy(d.indexsets), copy(d.indexexprs), copy(d.condition))
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
        typecode = :(type $(typename){T} <: JuMPArray{T}; innerArray::Array{T,$N}; name::String;
                            indexsets; indexexprs::Vector{IndexPair} end)
        getidxlhs = :(Base.getindex(d::$(typename)))
        setidxlhs = :(setindex!(d::$(typename),val))
        getidxrhs = :(Base.getindex(d.innerArray))
        setidxrhs = :(setindex!(d.innerArray,val))
        maplhs = :(Base.map(f::Function,d::$(typename)))
        maprhs = :($(typename)(map(f,d.innerArray),d.name,d.indexsets,d.indexexprs))
        for i in 1:N
            varname = symbol(string("x",i))
            
            push!(getidxlhs.args,:($varname))
            push!(setidxlhs.args,:($varname))

            push!(getidxrhs.args,:((isa($varname, Int) || isa($varname, Range)) ? $varname+$(offset[i]) : $varname ))
            push!(setidxrhs.args,:($varname+$(offset[i])))

        end

        badgetidxlhs = :(Base.getindex(d::$(typename),wrong...))
        badgetidxrhs = :(error("Wrong number of indices for ",d.name,
                               ", expected ",length(d.indexsets)))

        funcs = :($getidxlhs = $getidxrhs; $setidxlhs = $setidxrhs;
                  $maplhs = $maprhs; $badgetidxlhs = $badgetidxrhs)
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
for accessor in (:getValue, :getDual, :getLower, :getUpper)
    @eval $accessor(x::JuMPContainer) = map($accessor,x)
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
Base.size(x::JuMPArray,k) = size(x.innerArray,k)
function Base.next(x::JuMPArray,k)
    @assert (ndim = ndims(x)) > 0
    idxsets = x.indexsets
    lengths = map(length, idxsets)
    cprod = cumprod([lengths...])
    key = Array(Int, ndim)
    key[1] = idxsets[1][mod1(k,lengths[1])]
    for i in 2:ndim
        key[i] = idxsets[i][int(ceil(mod1(k,cprod[i])/cprod[i-1]))]
    end
    (var, gidx) = next(x.innerArray, k)
    return (tuple(key..., var), gidx)
end
function Base.next(x::JuMPDict,k)
    ((idx,var),gidx) = next(x.tupledict,k)
    return (tuple(idx..., var), gidx)
end
Base.done(x::JuMPArray,k) = done(x.innerArray,k)
Base.done(x::JuMPDict,k)  = done(x.tupledict,k)

(-)(x::JuMPArray,y::Array) = x.innerArray-y

Base.eltype{T}(x::JuMPContainer{T}) = T

Base.full(x::JuMPContainer) = x

export @gendict

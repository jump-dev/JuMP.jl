using Base.Meta
# multivarate "dictionary" used for collections of variables/constraints

abstract JuMPContainer

type JuMPDict{T,N} <: JuMPContainer
    tupledict::Dict{NTuple{N},T}
    name::Symbol
end

#JuMPDict{T,N}(name::String) =
#    JuMPDict{T,N}(Dict{NTuple{N},T}(), name)

Base.getindex(d::JuMPDict, t...) = d.tupledict[t]
Base.setindex!(d::JuMPDict, value, t...) = (d.tupledict[t] = value)

function Base.map(f::Function, d::JuMPDict)
    x = JuMPDict(d.name)
    for (k,v) in d.tupledict
        x.tupledict[k] = f(v)
    end
    return x
end

# generate and instantiate a type which is indexed by the given index sets
# the following types of index sets are allowed:
# 0:K -- range with compile-time starting index
# S -- general iterable set
macro gendict(instancename,T,idxsets...)
    N = length(idxsets)
    allranges = all([isexpr(s,:(:)) for s in idxsets])
    if allranges
        # JuMPArray
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
    else
        # JuMPDict
        return :(
            $(esc(instancename)) = JuMPDict{$T,$N}(Dict{NTuple{$N},$T}(),$(quot(instancename)))
        )
    end
end

# duck typing approach -- if eltype(innerArray) doesn't support accessor, will fail
for accessor in (:getValue, :getDual, :getLower, :getUpper)
    @eval $accessor(x::JuMPContainer) = map($accessor,x)
end

# delegate zero-argument functions
for f in (:(Base.endof), :(Base.ndims), :(Base.length), :(Base.abs))
    @eval $f(x::JuMPContainer)  = $f(x.innerArray)
end
# delegate one-argument functions
for f in (:(Base.size),)
    @eval $f(x::JuMPContainer,k) = $f(x.innerArray,k)
end

(-)(x::JuMPContainer,y::Array) = x.innerArray-y

Base.eltype{T}    (x::JuMPDict{T}) = T

Base.full(x::JuMPContainer) = x

export @gendict

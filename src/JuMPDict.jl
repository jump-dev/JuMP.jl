# multivarate "dictionary" used for collections of variables/constraints

abstract JuMPDict

# generate and instantiate a type which is indexed by the given index sets
# the following types of index sets are allowed:
# 0:K -- range with compile-time starting index
# S -- general iterable set
macro gendict(instancename,T,idxsets...)
    N = length(idxsets)
    typename = symbol(string("JuMPDict",gensym()))
    isrange = Array(Bool,N)
    offset = Array(Int,N)
    dictnames = Array(Symbol,N)
    for i in 1:N
        if isa(idxsets[i],Expr) && idxsets[i].head == :(:)
            isrange[i] = true
            if isa(idxsets[i].args[1],Int)
                offset[i] = 1 - idxsets[i].args[1]
            else
                error("Currently only ranges with integer compile-time starting values are allowed as index sets. $(idxsets[i].args[1]) is not an integer in range $(idxsets[i]).")
            end
        else
            isrange[i] = false
            dictnames[i] = gensym()
        end
    end
    #typecode = :(type $(esc(typename)); innerArray::Array{$T,$N}; name::String; end)
    typecode = :(type $(typename){T} <: JuMPDict; innerArray::Array{T,$N}; name::String; end)
    builddicts = quote end
    for i in 1:N
        if !isrange[i]
            push!(typecode.args[3].args,:($(dictnames[i])::Dict))
            push!(builddicts.args, quote 
                $(esc(dictnames[i])) = Dict(); 
                for (j,k) in enumerate($(esc(idxsets[i])))
                    $(esc(dictnames[i]))[k] = j
                end 
                end)
        end
    end
    getidxlhs = :(getindex(d::$(typename)))
    setidxlhs = :(setindex!(d::$(typename),val))
    getidxrhs = :(getindex(d.innerArray))
    setidxrhs = :(setindex!(d.innerArray,val))
    maplhs = :(mapvals(f,d::$(typename)))
    maprhs = :($(typename)(map(f,d.innerArray),d.name))
    for i in 1:N
        varname = symbol(string("x",i))
        
        if isrange[i]
            push!(getidxlhs.args,:($varname))
            push!(setidxlhs.args,:($varname))

            push!(getidxrhs.args,:($varname+$(offset[i])))
            push!(setidxrhs.args,:($varname+$(offset[i])))
        else
            push!(getidxlhs.args,varname)
            push!(setidxlhs.args,varname)

            push!(getidxrhs.args,:(d.($(Expr(:quote,dictnames[i])))[$varname]))
            push!(setidxrhs.args,:(d.($(Expr(:quote,dictnames[i])))[$varname]))
            push!(maprhs.args,:(d.($(Expr(:quote,dictnames[i])))))
        end
    end

    #funcs = :($(esc(getidxlhs)) = $(esc(getidxrhs)); $(esc(setidxlhs)) = $(esc(setidxrhs)))
    funcs = :($getidxlhs = $getidxrhs; $setidxlhs = $setidxrhs; $maplhs = $maprhs)
    geninstance = :($(esc(instancename)) = $(typename)(Array($T),$(string(instancename))))
    for i in 1:N
        push!(geninstance.args[2].args[2].args, :(length($(esc(idxsets[i])))))
        if !isrange[i]
            push!(geninstance.args[2].args, esc(dictnames[i]))
        end
    end
    eval(Expr(:toplevel, typecode))
    eval(Expr(:toplevel, funcs))

    quote
        $builddicts
        $geninstance
    end

end

# duck typing approach -- if eltype(innerArray) doesn't support accessor, will fail
for accessor in (:getValue, :getDual, :getLower, :getUpper)
    @eval $accessor(x::JuMPDict) = mapvals($accessor,x)
end

# delegate zero-argument functions
for f in (:endof, :ndims, :length, :abs)
    @eval $f(x::JuMPDict) = $f(x.innerArray)
end
# delegate one-argument functions
for f in (:size,)
    @eval $f(x::JuMPDict,k) = $f(x.innerArray,k)
end

(-)(x::JuMPDict,y::Array) = x.innerArray-y

export @gendict

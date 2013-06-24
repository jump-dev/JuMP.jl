# multivarate "dictionary" used for collections of variables/constraints

abstract MathProgDict

# generate and instantiate a type which is indexed by the given index sets
# the following types of index sets are allowed:
# 0:K -- range with compile-time starting index
# S -- general iterable set
macro gendict(instancename,T,idxsets...)
    N = length(idxsets)
    typename = symbol(string("MathProgDict",gensym()))
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
    typecode = :(type $(typename){T} <: MathProgDict; innerArray::Array{T,$N}; name::String; end)
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
    maplhs = :(map(f,d::$(typename)))
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

getValue(x::MathProgDict) = map(getValue,x)
endof(x::MathProgDict) = endof(x.innerArray)
ndims(x::MathProgDict) = ndims(x.innerArray)
size(x::MathProgDict,n) = size(x.innerArray,n)
length(x::MathProgDict) = length(x.innerArray)
abs(x::MathProgDict) = abs(x.innerArray)
(-)(x::MathProgDict,y::Array) = x.innerArray-y

export @gendict

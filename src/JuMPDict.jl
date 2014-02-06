using Base.Meta
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
        if isexpr(idxsets[i],:(:)) && length(idxsets[i].args) == 2 # don't yet optimize ranges with steps
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
    typecode = :(type $(typename){T} <: JuMPDict; innerArray::Array{T,$N}; name::String;
                        indexsets end)
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
    maprhs = :($(typename)(map(f,d.innerArray),d.name,d.indexsets))
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

    funcs = :($getidxlhs = $getidxrhs; $setidxlhs = $setidxrhs; $maplhs = $maprhs)
    geninstance = :($(esc(instancename)) = $(typename)(Array($T),$(string(instancename)),$(esc(Expr(:tuple,idxsets...)))))
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


function pretty(dict::JuMPDict, mode=:Con)
    m = dict.innerArray[1].m

    dimensions = length(dict.indexsets)
    if dimensions >= 5
        return ""  # Not enough indices!
    end

    # Check that bounds are same throughout
    colLow = m.colLower[dict.innerArray[1].col]
    colUp  = m.colUpper[dict.innerArray[1].col]
    all_same = true
    for v in dict.innerArray[2:end]
        all_same &= m.colLower[v.col] == colLow
        all_same &= m.colUpper[v.col] == colUp
        !all_same && break
    end
    if !all_same
        return ""  # The variables have different bounds, so can't handle
    end

    dim_names = ("i","j","k","l")

    # First, the central bit of the expression
    name_and_indices = "$(dict.name)"
    name_and_indices *= (mode == :Con) ? "[i" : "_{i"
    for dim = 2:dimensions
        name_and_indices *= ",$(dim_names[dim])"
    end
    name_and_indices *= (mode == :Con) ? "]" : "}"

    # Then the tail list of sets
    tail_str = (mode == :Con) ? ", for all " : "\\quad  \\forall "
    for dim = 1:dimensions
        tail_str *= "$(dim_names[dim])"
        tail_str *= (mode == :Con) ? " in {" : " \\in \\{ "
        tail_str *= "$(dict.indexsets[dim][1])..$(dict.indexsets[dim][end])"
        tail_str *= (mode == :Con) ? "}" : " \\}"
        if dim != dimensions
            tail_str *= ", "
        end
    end
    
    colCat = m.colCat[dict.innerArray[1].col]
    greater = (mode == :Con) ? "\u2265" : "\\geq"
    less = (mode == :Con) ? "\u2264" : "\\leq"
    if colCat == INTEGER && colLow == 0 && colUp == 1
        return name_and_indices * tail_str * ", binary"
    elseif colLow == -Inf && colUp == Inf
        return name_and_indices * tail_str * 
                (colCat == INTEGER ? ", integer" : ", free")
    elseif colLow == -Inf
        return name_and_indices * " $less $colUp" * tail_str * 
            (colCat == INTEGER ? ", integer" : "")
    elseif colUp == Inf
        return name_and_indices * " $greater $colLow" * tail_str * 
            (colCat == INTEGER ? ", integer" : "")
    else
        return "$colLow $less " * name_and_indices * " $less $colUp" * tail_str * 
            (colCat == INTEGER ? ", integer" : "")
    end
end
#############################################################################
# JuMP
# An algebraic modelling langauge for Julia
# See http://github.com/JuliaOpt/JuMP.jl
#############################################################################
# print.jl
# All "pretty printers" for JuMP types, including IJulia LaTeX support
# where possible.
# Note that, in general, these methods are not "fast" - in fact, they are
# quite slow. This is OK, as they are normally only used on smaller models
# and don't appear in code that needs to be fast, but they should not affect
# the speed of code materially that doesn't use it.
#############################################################################

#############################################################################
#### type MODEL

# checkNameStatus
# Not exported. Initially variables defined as indexed with @defVar do not
# have names because generating them is relatively slow compared to just
# creating the variable. We lazily generate them only when someone requires
# them. This function checks if we have done that already, and if not,
# generates them
function checkNameStatus(m::Model)
    for i in 1:m.numCols
        if m.colNames[i] == ""
            fillVarNames(m)
            return
        end
    end
end

function fillVarNames(m::Model)
    for dict in m.dictList
        fillVarNames(dict)
    end
end

function fillVarNames(v::JuMPArray{Variable})
    idxsets = v.indexsets
    lengths = map(length, idxsets)
    N = length(idxsets)
    name = v.name
    cprod = cumprod([lengths...])
    for (ind,var) in enumerate(v.innerArray)
        setName(var,string("$name[$(idxsets[1][mod1(ind,lengths[1])])", [ ",$(idxsets[i][int(ceil(mod1(ind,cprod[i]) / cprod[i-1]))])" for i=2:N ]..., "]"))
    end
end

function fillVarNames(v::JuMPDict{Variable})
    name = v.name
    for tmp in v
        ind, var = tmp[1:end-1], tmp[end]
        setName(var,string("$name[", join([string(i) for i in ind],","), "]"))
    end
end

# Default REPL
function Base.show(io::IO, m::Model)
    print(io, m.objSense == :Max ? "Maximization" : ((m.objSense == :Min && !isempty(m.obj)) ? "Minimization" : "Feasibility"))
    println(io, " problem with:")
    println(io, " * $(length(m.linconstr)) linear constraints")
    nquad = length(m.quadconstr)
    if nquad > 0
        println(io, " * $(nquad) quadratic constraints")
    end
    nlp = m.nlpdata
    if nlp != nothing && length(nlp.nlconstr) > 0
        println(io, " * $(length(nlp.nlconstr)) nonlinear constraints")
    end
    print(io, " * $(m.numCols) variables")
    nbin = sum(m.colCat .== :Bin)
    nint = sum(m.colCat .== :Int)
    nsc = sum(m.colCat .== :SemiCont)
    nsi = sum(m.colCat .== :SemiInt)
    varstr = {}
    nbin == 0 || push!(varstr, "$nbin binary")
    nint == 0 || push!(varstr, "$nint integer")
    nsc  == 0 || push!(varstr, "$nsc semicontinuous")
    nsi  == 0 || push!(varstr, "$nsi semi-integer")
    if isempty(varstr)
        println(io,)
    else
        println(io, ": $(join(varstr, ","))")
    end
    print(io, "Solver set to ")
    if isa(m.solver, UnsetSolver)
        solver = "Default"
    else
        solver = string(m.solver)
    end
    print(io, split(solver, "Solver")[1])
end

#############################################################################
# JuMPDict for values
Base.show(io::IO, dict::JuMPContainer{Float64}) = print(io, dict)
function Base.print(io::IO, dict::JuMPContainer{Float64})
    println(io, dict.name)
    print_values(io, dict, 1, {}, "")
end

function print_values(io::IO, dict::JuMPContainer{Float64}, depth::Int, 
                        parent_index::Vector{Any}, parent_str::String)
    dims = length(dict.indexsets)
    indexset = dict.indexsets[depth]

    # Turn index set into strings
    index_strs = map(string, indexset)

    # Determine longest index so we can align columns
    max_index_len = 0
    for index_str in index_strs
        max_index_len = max(max_index_len, strwidth(index_str))
    end

    # If have recursed, we need to prepend the parent's index strings
    # accumulated, as well as white space so the alignment works.
    for i = 1:length(index_strs)
        index_strs[i] = parent_str * lpad(index_strs[i],max_index_len," ")
    end

    # Create a string for the number of spaces we need to indent
    indent = " "^(2*(depth-1))

    # Determine the need to recurse
    if depth == dims
        # Deepest level
        for i = 1:length(indexset)
            value = length(parent_index) == 0 ? 
                        dict[indexset[i]] :
                        dict[parent_index...,indexset[i]]
            println(io, indent * "[" * index_strs[i] * "] = ", value)
        end
    else
        # At least one more layer to go
        for i = 1:length(indexset)
            index = indexset[i]
            # Print the ":" version of indices we will recurse over
            println(io, indent * "[" * index_strs[i] * ",:"^(dims-depth) * "]")
            print_values(io, dict, depth+1,
                 length(parent_index) == 0 ? {index} : {parent_index...,index},
                index_strs[i] * ",")
        end
    end
end

# support types that don't have built-in comparison
function _isless(t1::Tuple, t2::Tuple)
    n1, n2 = length(t1), length(t2)
    for i = 1:min(n1, n2)
        a, b = t1[i], t2[i]
        if !isequal(a, b)
            return applicable(isless,a,b) ? isless(a, b) : isless(hash(a),hash(b))
        end
    end
    return n1 < n2
end

# adapted from showdict in base/dict.jl
function Base.print(io::IO, dict::JuMPDict{Float64})
    rows, cols = Base.tty_size()[1] - 3, Base.tty_size()[2]
    nelem = length(dict.tupledict)
    print(io, "$(length(dict.indexsets))-dimensional JuMPDict with $nelem ")
    print(io, nelem == 1 ? "entry" : "entries")
    isempty(dict) && return 
    print(io, ":")

    rows < 2   && (print(io, " …"); return)
    cols < 12  && (cols = 12) # Minimum widths of 2 for key, 4 for value
    cols -= 6 # Subtract the widths of prefix "  " separator " => "
    rows -= 2 # Subtract the summary and final ⋮ continuation lines

    sortedkeys = sort(collect(keys(dict.tupledict)), lt = _isless)

    ks = Array(String, min(rows, length(dict)))
    keylen = 0
    for (i, key) in enumerate(sortedkeys)
        i > rows && break
        ks[i] = join(map(string,key),",")
        keylen = clamp(length(ks[i]), keylen, div(cols, 3))
    end

    for (i,key) in enumerate(sortedkeys)
        print(io, "\n ")
        i > rows && (print(io, rpad("⋮", keylen), " = ⋮"); break)
        v = dict[key...]
        print(io, dict.name, "[")
        print(io, rpad("$(ks[i])]", keylen+1))
        print(io, " = ")
        print(io, v) 
    end
end

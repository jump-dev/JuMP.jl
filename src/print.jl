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
        idxsets = dict.indexsets
        lengths = map(length, idxsets)
        N = length(idxsets)
        name = dict.name
        cprod = cumprod([lengths...])
        for (ind,var) in enumerate(dict.innerArray)
            setName(var,string("$name[$(idxsets[1][mod1(ind,lengths[1])])", [ ",$(idxsets[i][int(ceil(mod1(ind,cprod[i]) / cprod[i-1]))])" for i=2:N ]..., "]"))
        end
    end
end

# REPL printing
function print(io::IO, m::Model)
    checkNameStatus(m)

    println(io, string(m.objSense," ",quadToStr(m.obj)))
    println(io, "Subject to ")
    for c in m.linconstr
        println(io, conToStr(c))
    end
    for c in m.quadconstr
        println(io, conToStr(c))
    end

    # Handle special case of indexed variables
    in_dictlist = zeros(Bool, m.numCols)
    for dict in m.dictList
        out_str = dictstring(dict, :REPL)
        if out_str != ""
            println(io, out_str)
            # Don't repeat this variable
            for v in dict.innerArray
                in_dictlist[v.col] = true
            end
        end
    end

    # Handle all other variables
    for i in 1:m.numCols
        in_dictlist[i] && continue
        println(io, boundstring(m.colNames[i], m.colLower[i], m.colUpper[i],
                                m.colCat[i], "", :REPL))
    end
end

# IJulia printing
function writemime(io::IO, ::MIME"text/latex", m::Model)
    checkNameStatus(m)
  
    println(io, "\$\$")  # Begin MathJax mode
    println(io, "\\begin{alignat*}{1}")
    
    # Objective
    print(io, m.objSense == :Max ? "\\max \\quad &" : "\\min \\quad &")
    print(io, quadToStr(m.obj))
    println(io, "\\\\")
    # Constraints
    print(io, "\\text{Subject to} \\quad")
    for c in m.linconstr
        println(io, "& $(conToStr(c,true)) \\\\")
    end
    for c in m.quadconstr
        println(io, conToStr(c))
    end
    
    # Handle special cases
    in_dictlist = zeros(Bool, m.numCols)
    for dict in m.dictList
        out_str = dictstring(dict, :IJulia)
        if out_str != ""
            println(io, "& $out_str \\\\")
            # Don't repeat this variable
            for v in dict.innerArray
                in_dictlist[v.col] = true
            end
        end
    end
    for i in 1:m.numCols
        in_dictlist[i] && continue
        print(io, "& ")
        println(io, boundstring(m.colNames[i], m.colLower[i], m.colUpper[i],
                                m.colCat[i], "", :IJulia))
    end
    print(io, "\\end{alignat*}\n\$\$")
end

# Default REPL
function show(io::IO, m::Model)
    print(io, m.objSense == :Max ? "Maximization" : ((m.objSense == :Min && !isempty(m.obj)) ? "Minimization" : "Feasibility"))
    println(io, " problem with:")
    println(io, " * $(length(m.linconstr)) linear constraints")
    nquad = length(m.quadconstr)
    if nquad > 0
        println(io, " * $(nquad) quadratic constraints")
    end
    print(io, " * $(m.numCols) variables")  
    nint = sum(m.colCat .== INTEGER)
    println(io, nint == 0 ? "" : " ($nint integer)")
    print(io, "Solver set to ")
    if typeof(m.solver) == MissingSolver
        solver = nquad > 0 ? string(MathProgBase.defaultQPsolver) : (nint > 0 ? string(MathProgBase.defaultMIPsolver) : string(MathProgBase.defaultLPsolver))
    else
        solver = string(m.solver)
    end
    println(io, split(solver, "Solver")[1])
end

#############################################################################
#### type VARIABLE

# boundstring
# Not exported. Provides an output-specifc string of a variable and its
# bounds. iterate_over is only used for variables that are indexed (via the
# dictstring)
function boundstring(var_name, colLow, colUp, colCat, iterate_over="", mode=:REPL)
    greater = (mode == :REPL) ? "\u2265" : "\\geq"
    less    = (mode == :REPL) ? "\u2264" : "\\leq"
    int_str = (colCat == INTEGER) ? ", integer" : ""

    if colCat == INTEGER && colLow == 0 && colUp == 1
        return "$(var_name)$(iterate_over), binary"
    elseif colLow == -Inf && colUp == Inf
        return "$(var_name)$(iterate_over) free" * 
                (colCat == INTEGER ? ", integer" : "")
    elseif colLow == -Inf
        return "$var_name $less $(colUp)$(iterate_over)$(int_str)"
    elseif colUp == Inf
        return "$var_name $greater $(colLow)$(iterate_over)$(int_str)"
    else  # both bounds
        return "$colLow $less $var_name $less $(colUp)$(iterate_over)$(int_str)"
    end
end

# dictstring
# Not exported. Takes a JuMPDict and, if it can, builds a string that
# summarizes the variables into one line. If not, it will return an empty
# string and these variables should be printed one-by-one. Mode should be
# :REPL or :IJulia
function dictstring(dict::JuMPDict, mode=:REPL)
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

    name_and_indices, tail_str = dictnameindices(dict, mode)
    
    colCat = m.colCat[dict.innerArray[1].col]
    return boundstring(name_and_indices, colLow, colUp, colCat, tail_str, mode)
end

# dictnameindices
# Not exported. Builds the x[i,j,k,l] part and the "for all" parts. This is also
# used for printing JuMPDict so thats why its seperated out from dictstring
function dictnameindices(dict::JuMPDict, mode=:REPL)
    dimensions = length(dict.indexsets)
    dim_names = ("i","j","k","l")

    # The central bit of the expression
    name_and_indices = "$(dict.name)"
    name_and_indices *= (mode == :REPL) ? "[i" : "_{i"
    for dim = 2:dimensions
        name_and_indices *= ",$(dim_names[dim])"
    end
    name_and_indices *= (mode == :REPL) ? "]" : "}"

    # Then the tail list of sets
    tail_str = (mode == :REPL) ? ", for all " : " \\quad \\forall "
    for dim = 1:dimensions
        tail_str *= "$(dim_names[dim])"
        tail_str *= (mode == :REPL) ? " in {" : " \\in \\{ "
        if typeof(dict.indexsets[dim]) <: Range1
            # Range with increments of 1, easy
            tail_str *= "$(dict.indexsets[dim][1])..$(dict.indexsets[dim][end])"
        elseif typeof(dict.indexsets[dim]) <: Range
            # Range with increment that might not be 1
            r_first = first(dict.indexsets[dim])
            r_end   = last(dict.indexsets[dim])
            if length(dict.indexsets[dim]) == 1
                tail_str *= "$r_first"
            elseif length(dict.indexsets[dim]) == 2
                tail_str *= "$r_first,$r_end"
            elseif length(dict.indexsets[dim]) == 3
                r_second = dict.indexsets[dim][2]
                tail_str *= "$r_first,$r_second,$r_end"
            else
                r_second = dict.indexsets[dim][2]
                r_penult = dict.indexsets[dim][end-1]
                tail_str *= "$r_first,$r_second..$r_penult,$r_end"
            end
        else
            # Arbitrary set
            MAXCHAR = 15
            cur_str = ""
            for i in dict.indexsets[dim]
                str_i = string(i)
                if length(str_i) + length(cur_str) >= MAXCHAR
                    # Stop here
                    cur_str *= ".."
                    break
                else
                    # It will fit
                    if length(cur_str) > 0
                        cur_str *= ","
                    end
                    cur_str *= str_i 
                end
            end
            tail_str *= cur_str
        end
        tail_str *= (mode == :REPL) ? "}" : " \\}"
        if dim != dimensions
            tail_str *= ", "
        end
    end

    return name_and_indices, tail_str
end

show(io::IO, v::Variable) = print(io, getName(v))
print(io::IO, v::Variable) = print(io, getName(v))
writemime(io::IO, ::MIME"text/latex", v::Variable) = print(io, getName(v))

show(io::IO, dict::JuMPDict) = print(io, dict)
function print(io::IO, dict::JuMPDict)
    # Best case: bounds and all dims
    str = dictstring(dict::JuMPDict, :REPL)
    if str != ""
        print(io, str)    
        return
    end
    # Maybe too many dims?
    dimensions = length(dict.indexsets)
    if dimensions >= 5
        # Can't handle that many
        print(io, "$(dict.name)[...]")
        return
    end
    # Must have been inconsistent bounds - just don't print them
    name_and_indices, tail_str = dictnameindices(dict, :REPL)
    print(io, ".. \u2264 $(name_and_indices) \u2264 ..$(tail_str)")
end
function writemime(io::IO, ::MIME"text/latex", dict::JuMPDict)
    # Best case: bounds and all dims
    str = dictstring(dict::JuMPDict, :IJulia)
    if str != ""
        print(io, str)    
        return
    end
    # Maybe too many dims?
    dimensions = length(dict.indexsets)
    if dimensions >= 5
        # Can't handle that many
        print(io, "$(dict.name)[...]")
        return
    end
    # Must have been inconsistent bounds - just don't print them
    name_and_indices, tail_str = dictnameindices(dict, :REPL)
    print(io, "\$\$ \\dots \\leq $(name_and_indices) \\leq \\dots$(tail_str) \$\$")
end

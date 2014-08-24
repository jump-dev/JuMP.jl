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

# Helper function that rounds carefully for the purposes of printing
# e.g.   5.3  =>  5.3
#        1.0  =>  1
function string_intclamp(f::Float64)
    str = string(f)
    length(str) >= 2 && str[end-1:end] == ".0" ? str[1:end-2] : str
end

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
    for (ind,var) in v.tupledict
        setName(var,string("$name[", join([string(i) for i in ind],","), "]"))
    end
end

# REPL printing
function Base.print(io::IO, m::Model)
    checkNameStatus(m)

    nlp = m.nlpdata

    qobj_str = quadToStr(m.obj)
    if m.nlpdata != nothing && nlp.nlobj != nothing
        qobj_str = (qobj_str == "0" ? "" : qobj_str*" + ")
        println(io, string(m.objSense," ",qobj_str,"(nonlinear expression)"))
    else
        println(io, string(m.objSense," ",qobj_str))
    end
    println(io, "Subject to ")
    for c in m.linconstr
        println(io, conToStr(c))
    end
    for c in m.quadconstr
        println(io, conToStr(c))
    end
    for c in m.sosconstr
        println(io, conToStr(c))
    end
    if nlp != nothing && length(nlp.nlconstr) > 0
        if length(nlp.nlconstr) == 1
            println(io, "1 nonlinear constraint")
        else
            println(io, "$(length(m.nlpdata.nlconstr)) nonlinear constraints")
        end
    end

    # Handle special case of indexed variables
    in_dictlist = zeros(Bool, m.numCols)
    for dict in m.dictList
        if length(dict) == 0
            continue
        end
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
function Base.writemime(io::IO, ::MIME"text/latex", m::Model)
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
        println(io, "& $(conToStr(c)) \\\\")
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
        print(io, boundstring(m.colNames[i], m.colLower[i], m.colUpper[i],
                              m.colCat[i], "", :IJulia))
        println(io, "\\\\")
    end
    print(io, "\\end{alignat*}\n\$\$")
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
    nint = sum(m.colCat .== INTEGER)
    println(io, nint == 0 ? "" : " ($nint integer)")
    print(io, "Solver set to ")
    if isa(m.solver, UnsetSolver)
        solver = nquad > 0 ? string(MathProgBase.defaultQPsolver) : (nint > 0 ? string(MathProgBase.defaultMIPsolver) : string(MathProgBase.defaultLPsolver))
    else
        solver = string(m.solver)
    end
    print(io, split(solver, "Solver")[1])
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
        return "$var_name $less $(string_intclamp(colUp))$(iterate_over)$(int_str)"
    elseif colUp == Inf
        return "$var_name $greater $(string_intclamp(colLow))$(iterate_over)$(int_str)"
    else  # both bounds
        return "$(string_intclamp(colLow)) $less $var_name $less $(string_intclamp(colUp))$(iterate_over)$(int_str)"
    end
end

# dictstring
# Not exported. Takes a JuMPDict and, if it can, builds a string that
# summarizes the variables into one line. If not, it will return an empty
# string and these variables should be printed one-by-one. Mode should be
# :REPL or :IJulia
function dictstring(dict::JuMPContainer{Variable}, mode=:REPL)

    isempty(dict) && return ""

    v = first(dict)[2]

    m = v.m
    colCat = m.colCat[v.col]

    dimensions = length(dict.indexsets)
    if dimensions >= 5
        return ""  # Not enough indices!
    end

    # Check that bounds are same throughout
    colLow = m.colLower[v.col]
    colUp  = m.colUpper[v.col]
    all_same = true
    for el in dict
        v = el[2]
        all_same &= m.colLower[v.col] == colLow
        all_same &= m.colUpper[v.col] == colUp
        !all_same && break
    end
    if !all_same
        return ""  # The variables have different bounds, so can't handle
    end

    name_and_indices, tail_str = dictnameindices(dict, mode)
    
    return boundstring(name_and_indices, colLow, colUp, colCat, tail_str, mode)
end

# fillranges
# Not exported. Constructs the description of an arbitrary range. Output looks
# something like {2,4..18,20}
function fillranges(idx)
    r_first = first(idx)
    r_end   = last(idx)
    str = ""
    if length(idx) == 1
        str *= "$r_first"
    elseif length(idx) == 2
        str *= "$r_first,$r_end"
    elseif length(idx) == 3
        r_second = idx[2]
        str *= "$r_first,$r_second,$r_end"
    else
        r_second = idx[2]
        r_penult = idx[end-1]
        str *= "$r_first,$r_second..$r_penult,$r_end"
    end
    return str
end

# dictnameindices
# Not exported. Builds the x[i,j,k,l] part and the "for all" parts. This is also
# used for printing JuMPDict so thats why its separated out from dictstring
function dictnameindices(dict::JuMPContainer{Variable}, mode=:REPL)
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
            tail_str *= fillranges(dict.indexsets[dim])
        else
            try # try to detect ranges in disguise
                elem = dict.indexsets[dim][1]
                off = dict.indexsets[dim][2] - elem
                for k in dict.indexsets[dim][2:end]
                    if (k - elem) != off
                        error("Internal error")
                    end
                    elem = k
                end
                if off == 1
                    tail_str *= "$(dict.indexsets[dim][1])..$(dict.indexsets[dim][end])"
                else
                    tail_str *= fillranges(dict.indexsets[dim])
                end
            catch # Arbitrary set
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
        end
        tail_str *= (mode == :REPL) ? "}" : " \\}"
        if dim != dimensions
            tail_str *= ", "
        end
    end
    if isa(dict, JuMPDict) && !isempty(dict.conditions)
        tail_str *= " s.t. $(join(dict.conditions, " and "))"
    end

    return name_and_indices, tail_str
end

##########
# Variable
Base.show(io::IO, v::Variable) = print(io, getName(v))
Base.print(io::IO, v::Variable) = print(io, getName(v))
Base.writemime(io::IO, ::MIME"text/latex", v::Variable) = print(io, getName(v))

#########
# AffExpr
Base.print(io::IO, a::GenericAffExpr) = print(io, affToStr(a))
Base.show(io::IO, a::GenericAffExpr) = print(io, affToStr(a))

function affToStr(a::AffExpr, showConstant=true)
    if length(a.vars) == 0
        if showConstant
            return string_intclamp(a.constant)
        else
            return "0"
        end
    end

    # Get reference to models
    moddict = Dict{Model,IndexedVector}()
    for var in a.vars
        mod = var.m
        if !haskey(moddict, mod)
            checkNameStatus(mod)
            moddict[var.m] = IndexedVector(Float64,mod.numCols)
        end
    end

    # Collect like terms
    for ind in 1:length(a.vars)
        addelt!(moddict[a.vars[ind].m], a.vars[ind].col, a.coeffs[ind])
    end

    elm = 0
    termStrings = Array(UTF8String, 2*length(a.vars))
    for m in keys(moddict)
        indvec = moddict[m]
        for i in 1:indvec.nnz
            idx = indvec.nzidx[i]
            if abs(abs(indvec.elts[idx])-1) < 1e-20
                if elm == 0
                    elm += 1
                    if indvec.elts[idx] < 0
                        termStrings[1] = "-$(getName(m,idx))"
                    else
                        termStrings[1] = "$(getName(m,idx))"
                    end
                else 
                    if indvec.elts[idx] < 0
                        termStrings[2*elm] = " - "
                    else
                        termStrings[2*elm] = " + "
                    end
                    termStrings[2*elm+1] = "$(getName(m,idx))"
                    elm += 1
                end
            elseif abs(indvec.elts[idx]) > 1e-20
                if elm == 0
                    elm += 1
                    termStrings[1] = "$(string_intclamp(indvec.elts[idx])) $(getName(m,idx))"
                else 
                    if indvec.elts[idx] < 0
                        termStrings[2*elm] = " - "
                    else
                        termStrings[2*elm] = " + "
                    end
                    termStrings[2*elm+1] = "$(string_intclamp(abs(indvec.elts[idx]))) $(getName(m,idx))"
                    elm += 1
                end
            end
        end
    end

    if elm == 0
        ret = "0"
    else
        # And then connect them up with +s
        ret = join(termStrings[1:(2*elm-1)])
    end
    
    if abs(a.constant) >= 0.000001 && showConstant
        if a.constant < 0
            ret = string(ret, " - ", string_intclamp(abs(a.constant)))
        else
            ret = string(ret, " + ", string_intclamp(a.constant))
        end
    end
    return ret
end

##########
# QuadExpr
##########
Base.print(io::IO, q::GenericQuadExpr) = print(io, quadToStr(q))
Base.show(io::IO, q::GenericQuadExpr)  = print(io, quadToStr(q))

function quadToStr(q::QuadExpr)
    if length(q.qvars1) == 0
        return affToStr(q.aff)
    end

    # canonicalize and merge duplicates
    for ind in 1:length(q.qvars1)
            if q.qvars1[ind].m != q.qvars2[ind].m
                error("You cannot have a quadratic term with variables from different models")
            end
            if q.qvars2[ind].col < q.qvars1[ind].col
                    q.qvars1[ind],q.qvars2[ind] = q.qvars2[ind],q.qvars1[ind]
            end
    end
    Q = sparse([v.col for v in q.qvars1], [v.col for v in q.qvars2], q.qcoeffs)
    I,J,V = findnz(Q)
    Qnnz = length(V)

    termStrings = Array(UTF8String, 2*Qnnz)
    if Qnnz > 0
        if V[1] < 0
            termStrings[1] = "-"
        else
            termStrings[1] = ""
        end
        for ind in 1:Qnnz
            if ind >= 2
                if V[ind] < 0
                    termStrings[2*ind-1] = " - "
                else 
                    termStrings[2*ind-1] = " + "
                end
            end
            x = Variable(q.qvars1[ind].m,I[ind])
            if I[ind] == J[ind]
                # Squared term
                if abs(V[ind]) == 1.0
                    termStrings[2*ind] = string(getName(x),"\u00B2")
                else
                    termStrings[2*ind] = string(string_intclamp(abs(V[ind]))," ",
                                                getName(x),"\u00B2")
                end
            else
                # Normal term
                y = Variable(q.qvars1[ind].m,J[ind])
                if abs(V[ind]) == 1.0
                    termStrings[2*ind] = string(getName(x),"*",getName(y))
                else
                    termStrings[2*ind] = string(string_intclamp(abs(V[ind]))," ",
                                                getName(x),"*",getName(y))
                end
            end
        end
    end
    ret = join(termStrings)

    if q.aff.constant == 0 && length(q.aff.vars) == 0
        return ret
    else
        aff = affToStr(q.aff)
        if aff[1] == '-'
            return string(ret, " - ", aff[2:end])
        else
            return string(ret, " + ", aff)
        end
    end
end

#############################################################################
# JuMPDict for variables
Base.show(io::IO, dict::JuMPContainer{Variable}) = print(io, dict)
function Base.print(io::IO, dict::JuMPContainer{Variable})
    # Best case: bounds and all dims
    str = dictstring(dict, :REPL)
    if str != ""
        print(io, str)    
        return
    end
    # Easy case: empty JuMPDict
    isempty(dict) && return nothing
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

function Base.writemime(io::IO, ::MIME"text/latex", dict::JuMPContainer{Variable})
    # Best case: bounds and all dims
    str = dictstring(dict::JuMPDict, :IJulia)
    if str != ""
        print(io, "\\( $str \\)")
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
    name_and_indices, tail_str = dictnameindices(dict, :IJulia)
    print(io, "\\( \\dots \\leq $(name_and_indices) \\leq \\dots$(tail_str) \\)")
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



###################
# Linear Constraint (or rather, the general version of it)
Base.print(io::IO, c::GenericRangeConstraint) = print(io, conToStr(c))
Base.show(io::IO,  c::GenericRangeConstraint) = print(io, conToStr(c))

function conToStr(c::GenericRangeConstraint)
    s = sense(c)
    if s == :range
        return string(string_intclamp(c.lb)," <= ",affToStr(c.terms,false)," <= ",string_intclamp(c.ub))
    else
        return string(affToStr(c.terms,false)," ",s," ",string_intclamp(rhs(c)))
    end
end

#################
# Quad Constraint
Base.print(io::IO, c::QuadConstraint) = print(io, conToStr(c))
Base.show(io::IO, c::QuadConstraint)  = print(io, conToStr(c))

conToStr(c::QuadConstraint) = string(quadToStr(c.terms), " ", c.sense, " 0")

################
# SOS Constraint
function conToStr(c::SOSConstraint) 
    nvar = length(c.terms)
    termStrings = Array(UTF8String, nvar)
    # termStrings = Array(UTF8String, nvar+2)
    # termStrings[1] = "$(c.sostype): {"
    for i in 1:nvar
        termStrings[i] = "$(c.weights[i]) $(c.terms[i])"
    end
    # if nvar > 0
    #     termStrings[2] = "$(c.weights[1]) $(c.terms[1])"
    #     for i in 2:nvar
    #         termStrings[i+1] = ", $(c.weights[i]) $(c.terms[i])"
    #     end
    # end
    # termStrings[end] = "}"
    return string("$(c.sostype): {", join(termStrings,", "), "}")
end

Base.print(io::IO, c::SOSConstraint) = print(io, conToStr(c))
Base.show(io::IO, c::SOSConstraint)  = print(io, conToStr(c))


################
# Constraint Ref

Base.print(io::IO, c::ConstraintRef{LinearConstraint}) = print(io, conToStr(c.m.linconstr[c.idx]))
Base.print(io::IO, c::ConstraintRef{QuadConstraint}) = print(io, conToStr(c.m.quadconstr[c.idx]))
Base.print(io::IO, c::ConstraintRef{SOSConstraint}) = print(io, conToStr(c.m.sosconstr[c.idx]))
Base.show{T}(io::IO, c::ConstraintRef{T}) = print(io, c)

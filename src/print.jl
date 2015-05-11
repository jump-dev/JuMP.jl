#  Copyright 2015, Iain Dunning, Joey Huchette, Miles Lubin, and contributors
#  This Source Code Form is subject to the terms of the Mozilla Public
#  License, v. 2.0. If a copy of the MPL was not distributed with this
#  file, You can obtain one at http://mozilla.org/MPL/2.0/.
#############################################################################
# JuMP
# An algebraic modelling langauge for Julia
# See http://github.com/JuliaOpt/JuMP.jl
#############################################################################
# print.jl
# All "pretty printers" for JuMP types.
# - Delegates to appropriate handler methods for REPL or IJulia.
# - These handler methods then pass the correct symbols to use into a
#   generic string builder. The IJulia handlers will also wrap in MathJax
#   start/close tags.
# - To find printing code for a type in this file, search for `## TypeName`
# - Code here does not need to be fast, in fact simplicity trumps speed
#   within reason as this code is thorny enough as it is.
# - Corresponding tests are in test/print.jl, although test/operator.jl
#   is also testing the constraint/expression code extensively as well.
#############################################################################

# Used for dispatching
abstract PrintMode
abstract REPLMode <: PrintMode
abstract IJuliaMode <: PrintMode

# Whether something is zero or not for the purposes of printing it
const PRINT_ZERO_TOL = 1e-10

# List of indices available for variable printing
const DIMS = ["i","j","k","l","m","n"]

# Helper function that rounds carefully for the purposes of printing
# e.g.   5.3  =>  5.3
#        1.0  =>  1
function str_round(f::Float64)
    abs(f) == 0.0 && return "0" # strip sign off zero
    str = string(f)
    length(str) >= 2 && str[end-1:end] == ".0" ? str[1:end-2] : str
end

# REPL-specific symbols
const repl_leq = @windows? "<=" : "≤"
const repl_geq = @windows? ">=" : "≥"
const repl_eq  = @windows? "==" : "="
const repl_times = "*"
const repl_sq    = "\u00B2"  # Superscript 2
const repl_ind_open  = "["
const repl_ind_close = "]"
const repl_for_all   = "for all"
const repl_in        = "in"
const repl_open_set  = "{"
const repl_mid_set   = ".."
const repl_close_set = "}"
const repl_union     = "or"
const repl_infty     = "Inf"
const repl_open_rng  = "["
const repl_close_rng = "]"
const repl_integer   = "integer"

# IJulia-specific symbols
const ijulia_leq        = "\\leq"
const ijulia_geq        = "\\geq"
const ijulia_eq         = "="
const ijulia_times      = "\\times"
const ijulia_sq         = "^2"
const ijulia_ind_open   = "_{"
const ijulia_ind_close  = "}"
const ijulia_for_all    = "\\quad\\forall"
const ijulia_in         = "\\in"
const ijulia_open_set   = "\\{"
const ijulia_mid_set    = ",\\dots,"
const ijulia_close_set  = "\\}"
const ijulia_union      = "\\cup"
const ijulia_infty      = "\\intfy"
const ijulia_open_rng  = "\\["
const ijulia_close_rng = "\\]"
const ijulia_integer   = "\\in \\mathbb{Z}"

# If not already mathmode, then wrap in MathJax start/close tags
math(s,mathmode) = mathmode ? s : "\$\$ $s \$\$"

#------------------------------------------------------------------------
## Model
#------------------------------------------------------------------------
function Base.print(io::IO, m::Model; ignore_print_hook=(m.printhook==nothing))
    ignore_print_hook || return m.printhook(io, m)
    print(io, model_str(REPLMode,m))
end

function Base.show(io::IO, m::Model)
    plural(n) = (n==1 ? "" : "s")
    print(io, m.objSense == :Max ? "Maximization" : ((m.objSense == :Min && (!isempty(m.obj) || (m.nlpdata != nothing && isa(m.nlpdata.nlobj, ReverseDiffSparse.SymbolicOutput)))) ? "Minimization" : "Feasibility"))
    println(io, " problem with:")
    nlin = length(m.linconstr)
    println(io, " * $(nlin) linear constraint$(plural(nlin))")
    nquad = length(m.quadconstr)
    if nquad > 0
        println(io, " * $(nquad) quadratic constraint$(plural(nquad))")
    end
    nlp = m.nlpdata
    if nlp != nothing && length(nlp.nlconstr) > 0
        println(io, " * $(length(nlp.nlconstr)) nonlinear constraint$(plural(length(nlp.nlconstr)))")
    end
    print(io, " * $(m.numCols) variable$(plural(m.numCols))")
    nbin = sum(m.colCat .== :Bin)
    nint = sum(m.colCat .== :Int)
    nsc = sum(m.colCat .== :SemiCont)
    nsi = sum(m.colCat .== :SemiInt)
    varstr = Any[]
    nbin == 0 || push!(varstr, "$nbin binary")
    nint == 0 || push!(varstr, "$nint integer")
    nsc  == 0 || push!(varstr, "$nsc semicontinuous")
    nsi  == 0 || push!(varstr, "$nsi semi-integer")
    if isempty(varstr)
        println(io,)
    else
        println(io, ": $(join(varstr, ", "))")
    end
    print(io, "Solver set to ")
    if isa(m.solver, UnsetSolver)
        solver = "Default"
    else
        solver = string(m.solver)
    end
    print(io, split(solver, "Solver")[1])
end
Base.writemime(io::IO, ::MIME"text/latex", m::Model) =
    print(io, model_str(IJuliaMode,m))
function model_str(mode, m::Model, leq, geq, in_set,
                            open_set, mid_set, close_set, union, infty,
                            open_rng, close_rng, integer)
    ijl = mode == IJuliaMode
    sep = ijl ? " & " : " "
    eol = ijl ? "\\\\\n" : "\n"
    nlp = m.nlpdata


    # Objective
    qobj_str = quad_str(mode, m.obj)
    obj_sense = ijl ? (m.objSense == :Max ? "\\max" : "\\min")*"\\quad" :
                      (m.objSense == :Max ? "Max" : "Min")
    str = obj_sense * sep
    if nlp != nothing && nlp.nlobj != nothing
        str *= (qobj_str=="0"?"":"$qobj_str + ") * "(nonlinear expression)"
    else
        str *= qobj_str
    end
    str *= eol

    # Constraints
    str *= ijl ? "\\text{Subject to} \\quad" : "Subject to" * eol
    for c in m.linconstr
        str *= sep * con_str(mode,c,mathmode=true) * eol
    end
    for c in m.quadconstr
        str *= sep * con_str(mode,c,mathmode=true) * eol
    end
    for c in m.sosconstr
        str *= sep * con_str(mode,c,mathmode=true) * eol
    end
    if nlp != nothing && length(nlp.nlconstr) > 0
        num = length(nlp.nlconstr)
        str *= sep * string("$num nonlinear constraint", num>1?"s":"") * eol
    end

    # Display indexed variables
    in_dictlist = falses(m.numCols)
    for d in m.dictList
        str *= sep * cont_str(mode,d,mathmode=true)  * eol

        # make sure that you haven't changed a variable type in the collection
        cat = getCategory(first(d)[end])
        allsame = true
        for v in d
            if getCategory(v[end]) != cat
                allsame = false
                break
            end
        end
        if allsame
            for it in d  # Mark variables in JuMPContainer as printed
                in_dictlist[it[end].col] = true
            end
        end
    end

    # Display non-indexed variables
    for i in 1:m.numCols
        in_dictlist[i] && continue
        var_name = var_str(mode,m,i)
        var_lb, var_ub = m.colLower[i], m.colUpper[i]
        str_lb, str_ub = str_round(var_lb), str_round(var_ub)
        var_cat = m.colCat[i]
        if var_cat == :Bin  # x binary
            str *= sep * "$var_name $in_set $(open_set)0,1$close_set"
        elseif var_cat == :SemiInt  # x in union of 0 and {lb,...,ub}
            str *= sep * "$var_name $in_set $open_set$str_lb$mid_set$str_ub$close_set $union $(open_set)0$close_set"
        elseif var_cat == :SemiCont  # x in union of 0 and [lb,ub]
            str *= sep * "$var_name $in_set $open_rng$str_lb,$str_ub$close_rng $union $(open_set)0$close_set"
        elseif var_cat == :Fixed
            str *= sep * "$var_name = $str_lb"
        elseif var_lb == -Inf && var_ub == +Inf # Free variable
            str *= sep * "$var_name free"
        elseif var_lb == -Inf  # No lower bound
            str *= sep * "$var_name $leq $str_ub"
        elseif var_ub == +Inf  # No upper bound
            str *= sep * "$var_name $geq $str_lb"
        else
            str *= sep * "$str_lb $leq $var_name $leq $str_ub"
        end
        if var_cat == :Int
            str *= ", $integer"
        end
        str *= eol
    end

    ijl ? "\$\$ \\begin{alignat*}{1}"*str*"\\end{alignat*}\n \$\$" :
          str
end

# Handlers to use correct symbols
model_str(::Type{REPLMode}, m::Model) =
    model_str(REPLMode, m, repl_leq, repl_geq, repl_in,
                        repl_open_set, repl_mid_set, repl_close_set,
                        repl_union, repl_infty, repl_open_rng, repl_close_rng,
                        repl_integer)
model_str(::Type{IJuliaMode}, m::Model; mathmode=true) =
    math(model_str(IJuliaMode, m, ijulia_leq, ijulia_geq, ijulia_in,
                        ijulia_open_set, ijulia_mid_set, ijulia_close_set,
                        ijulia_union, ijulia_infty, ijulia_open_rng, ijulia_close_rng,
                        ijulia_integer), mathmode)


#------------------------------------------------------------------------
## Variable
#------------------------------------------------------------------------
Base.print(io::IO, v::Variable) = print(io, var_str(REPLMode,v))
Base.show( io::IO, v::Variable) = print(io, var_str(REPLMode,v))
Base.writemime(io::IO, ::MIME"text/latex", v::Variable) =
    print(io, var_str(IJuliaMode,v,mathmode=false))
function var_str(mode, m::Model, col::Int, ind_open, ind_close)
    colNames = mode == REPLMode ? m.colNames : m.colNamesIJulia
    if colNames[col] == ""
        for cont in m.dictList
            fill_var_names(mode, colNames, cont)
        end
    end
    return colNames[col] == "" ? "col_$col" : colNames[col]
end
function fill_var_names(mode, colNames, v::JuMPArray{Variable})
    idxsets = v.indexsets
    lengths = map(length, idxsets)
    N = length(idxsets)
    name = v.name
    cprod = cumprod([lengths...])
    for (ind,var) in enumerate(v.innerArray)
        idx_strs = [string( idxsets[1][mod1(ind,lengths[1])] )]
        for i = 2:N
            push!(idx_strs, string(idxsets[i][int(ceil(mod1(ind,cprod[i]) / cprod[i-1]))]))
        end
        if mode == IJuliaMode
            colNames[var.col] = string(name, "_{", join(idx_strs,",") , "}")
        else
            colNames[var.col] = string(name,  "[", join(idx_strs,",") , "]")
        end
    end
end
function fill_var_names(mode, colNames, v::JuMPDict{Variable})
    name = v.name
    for tmp in v
        ind, var = tmp[1:end-1], tmp[end]
        if mode == IJuliaMode
            colNames[var.col] = string(name, "_{", join([string(i) for i in ind],","), "}")
        else
            colNames[var.col] = string(name,  "[", join([string(i) for i in ind],","), "]")
        end
    end
end

# Handlers to use correct symbols
var_str(::Type{REPLMode}, v::Variable) =
    var_str(REPLMode, v.m, v.col)
var_str(::Type{IJuliaMode}, v::Variable; mathmode=true) =
    var_str(IJuliaMode, v.m, v.col, mathmode=mathmode)

var_str(::Type{REPLMode}, m::Model, col::Int) =
    var_str(REPLMode, m, col, repl_ind_open, repl_ind_close)
var_str(::Type{IJuliaMode}, m::Model, col::Int; mathmode=true) =
    math(var_str(IJuliaMode, m, col, ijulia_ind_open, ijulia_ind_close), mathmode)

#------------------------------------------------------------------------
## JuMPContainer{Variable}
#------------------------------------------------------------------------
Base.print(io::IO, j::JuMPContainer{Variable}) = print(io, cont_str(REPLMode,j))
Base.show( io::IO, j::JuMPContainer{Variable}) = print(io, cont_str(REPLMode,j))
Base.writemime(io::IO, ::MIME"text/latex", j::JuMPContainer{Variable}) =
    print(io, cont_str(IJuliaMode,j,mathmode=false))
# Generic string converter, called by mode-specific handlers
function cont_str(mode, j::JuMPContainer{Variable}, leq, eq, geq,
                            ind_open, ind_close, for_all, in_set,
                            open_set, mid_set, close_set, union, infty,
                            open_rng, close_rng, integer)
    # Check if anything in the container
    isempty(j) && return string(j.name, " (no indices)")

    # 1. construct the part with variable name and indexing
    locvars = map(j.indexexprs) do tmp
        var = tmp.idxvar
        if var == nothing
            return ""
        else
            return string(var)
        end
    end
    num_dims = length(j.indexsets)
    idxvars = Array(UTF8String, num_dims)
    dimidx = 1
    for i in 1:num_dims
        if j.indexexprs[i].idxvar == nothing
            while DIMS[dimidx] in locvars
                dimidx += 1
            end
            if dimidx > length(DIMS)
                error("Unexpectedly ran out of indices")
            end
            idxvars[i] = DIMS[dimidx]
            dimidx += 1
        else
            idxvars[i] = locvars[i]
        end
    end
    name_idx = string(j.name, ind_open, join(idxvars,","), ind_close)
    # 2. construct part with what we index over
    idx_sets = for_all*" "*join(map(dim->string(idxvars[dim], " ", in_set, " ", open_set,
                                cont_str_set(j.indexsets[dim], mid_set),
                                close_set), 1:num_dims), ", ")
    # 3. Handle any conditions
    if isa(j, JuMPDict) && j.condition != :()
       idx_sets *= " s.t. $(join(parse_conditions(j.condition), " and "))"
    end

    # 4. Bounds and category, if possible, and return final string
    a_var = first(j)[end]
    model = a_var.m
    var_cat = model.colCat[a_var.col]
    var_lb  = model.colLower[a_var.col]
    var_ub  = model.colUpper[a_var.col]
    # Variables may have different bounds, so we can't really print nicely
    # at this time (possibly ever, as they could have been changed post
    # creation, which we'd never be able to handle.
    all_same_lb = true
    all_same_ub = true
    for iter in j
        var = iter[end]
        all_same_lb &= model.colLower[var.col] == var_lb
        all_same_ub &= model.colUpper[var.col] == var_ub
    end
    str_lb = var_lb == -Inf ? "-$infty" : str_round(var_lb)
    str_ub = var_ub == +Inf ? infty     : str_round(var_ub)
    # Special case bounds printing based on the category
    if var_cat == :Bin  # x in {0,1}
        return "$name_idx $in_set $(open_set)0,1$close_set $idx_sets"
    elseif var_cat == :SemiInt  # x in union of 0 and {lb,...,ub}
        si_lb = all_same_lb ? str_lb : ".."
        si_ub = all_same_ub ? str_ub : ".."
        return "$name_idx $in_set $open_set$si_lb$mid_set$si_ub$close_set $union $(open_set)0$close_set $idx_sets"
    elseif var_cat == :SemiCont  # x in union of 0 and [lb,ub]
        si_lb = all_same_lb ? str_lb : ".."
        si_ub = all_same_ub ? str_ub : ".."
        return "$name_idx $in_set $open_rng$si_lb,$si_ub$close_rng $union $(open_set)0$close_set $idx_sets"
    elseif var_cat == :Fixed
        si_bnd = all_same_lb ? str_lb : ".."
        return "$name_idx = $si_bnd $idx_sets"
    end
    # Continuous and Integer
    idx_sets = var_cat == :Int ? ", $integer, $idx_sets" : " $idx_sets"
    if all_same_lb && all_same_ub
        # Free variable
        var_lb == -Inf && var_ub == +Inf && return "$name_idx free$idx_sets"
        # No lower bound
        var_lb == -Inf && return "$name_idx $leq $str_ub$idx_sets"
        # No upper bound
        var_ub == +Inf && return "$name_idx $geq $str_lb$idx_sets"
        # Range
        return "$str_lb $leq $name_idx $leq $str_ub$idx_sets"
    end
    if all_same_lb && !all_same_ub
        var_lb == -Inf && return "$name_idx $leq ..$idx_sets"
        return "$str_lb $leq $name_idx $leq ..$idx_sets"
    end
    if !all_same_lb && all_same_ub
        var_ub == +Inf && return "$name_idx $geq ..$idx_sets"
        return ".. $leq $name_idx $leq $str_ub$idx_sets"
    end
    return ".. $leq $name_idx $leq ..$idx_sets"
end
# UTILITY FUNCTIONS FOR cont_str
function cont_str_set(idxset::Union(Range,Array), mid_set)  # 2:2:20 -> {2,4..18,20}
    length(idxset) == 1 && return string(idxset[1])
    length(idxset) == 2 && return string(idxset[1],",",idxset[2])
    length(idxset) == 3 && return string(idxset[1],",",idxset[2],",",idxset[3])
    length(idxset) == 4 && return string(idxset[1],",",idxset[2],",",idxset[3],",",idxset[4])
    return string(idxset[1],",",idxset[2],mid_set,idxset[end-1],",",idxset[end])
end
cont_str_set(idxset, mid_set) = return ".." # Fallback
# parse_conditions
# Not exported. Traverses an expression and constructs an array with entries
# corresponding to each condition. More specifically, if the condition is
# a && (b || c) && (d && e), it returns [a, b || c, d, e].
parse_conditions(not_an_expr) = not_an_expr
function parse_conditions(expr::Expr)
    ret = Any[]
    if expr.head != :&&
        return push!(ret, expr)
    end
    recurse = map(parse_conditions, expr.args)
    vcat(ret, recurse...)
end

# Handlers to use correct symbols
cont_str(::Type{REPLMode}, j::JuMPContainer{Variable}; mathmode=false) =
    cont_str(REPLMode, j, repl_leq, repl_eq, repl_geq, repl_ind_open, repl_ind_close,
                repl_for_all, repl_in, repl_open_set, repl_mid_set, repl_close_set,
                repl_union, repl_infty, repl_open_rng, repl_close_rng, repl_integer)
cont_str(::Type{IJuliaMode}, j::JuMPContainer{Variable}; mathmode=true) =
    math(cont_str(IJuliaMode, j, ijulia_leq, ijulia_eq, ijulia_geq, ijulia_ind_open, ijulia_ind_close,
                ijulia_for_all, ijulia_in, ijulia_open_set, ijulia_mid_set, ijulia_close_set,
                ijulia_union, ijulia_infty, ijulia_open_rng, ijulia_close_rng, ijulia_integer), mathmode)

#------------------------------------------------------------------------
## JuMPContainer{Float64}
#------------------------------------------------------------------------
Base.print(io::IO, j::JuMPContainer{Float64}) = print(io, val_str(REPLMode,j))
Base.show( io::IO, j::JuMPContainer{Float64}) = print(io, val_str(REPLMode,j))
function val_str(mode, j::JuMPArray{Float64})
    dims = length(j.indexsets)
    out_str = "$(j.name): $dims dimensions:\n"

    function val_str_rec(depth, parent_index::Vector{Any}, parent_str::String)
        # Turn index set into strings
        indexset = j.indexsets[depth]
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
                            j[indexset[i]] :
                            j[parent_index...,indexset[i]]
                out_str *= indent * "[" * index_strs[i] * "] = $value\n"
            end
        else
            # At least one more layer to go
            for i = 1:length(indexset)
                index = indexset[i]
                # Print the ":" version of indices we will recurse over
                out_str *= indent * "[" * index_strs[i] * ",:"^(dims-depth) * "]\n"
                val_str_rec(depth+1,
                     length(parent_index) == 0 ? Any[index] : Any[parent_index...,index],
                    index_strs[i] * ",")
            end
        end
    end
    val_str_rec(1,Any[],"")
    return out_str
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
function val_str(mode, dict::JuMPDict{Float64})
    nelem = length(dict.tupledict)
    out_str  = "$(dict.name): $(length(dict.indexsets)) dimensions, $nelem "
    out_str *= nelem == 1 ? "entry" : "entries"
    isempty(dict) && return out_str
    out_str *= ":"

    sortedkeys = sort(collect(keys(dict.tupledict)), lt = _isless)

    ndim = length(first(keys(dict.tupledict)))

    key_strs = Array(String, length(dict), ndim)
    for (i, key) in enumerate(sortedkeys)
        for j in 1:ndim
            key_strs[i,j] = string(key[j])
        end
    end
    max_dim_lens = map(1:ndim) do i
        maximum(map(length,key_strs[:,i]))
    end
    key_str = map(1:length(dict)) do i
        join(map(1:ndim) do j
            lpad(key_strs[i,j], max_dim_lens[j])
        end, ",")
    end
    max_key_len = maximum(map(length,key_str))

    for (i,key) in enumerate(sortedkeys)
        val = dict[key...]
        out_str *= "\n" * lpad("[$(key_str[i])]", max_key_len+3)
        out_str *= " = $val"
    end
    return out_str
end


#------------------------------------------------------------------------
## AffExpr  (not GenericAffExpr)
#------------------------------------------------------------------------
Base.print(io::IO, a::AffExpr) = print(io, aff_str(REPLMode,a))
Base.show( io::IO, a::AffExpr) = print(io, aff_str(REPLMode,a))
Base.writemime(io::IO, ::MIME"text/latex", a::AffExpr) =
    print(io, math(aff_str(IJuliaMode,a),false))
# Generic string converter, called by mode-specific handlers
function aff_str(mode, a::AffExpr; show_constant=true)
    # If the expression is empty, return the constant (or 0)
    if length(a.vars) == 0
        return show_constant ? str_round(a.constant) : "0"
    end

    # Get reference to models included in this expression
    moddict = Dict{Model,IndexedVector}()
    for var in a.vars
        if !haskey(moddict, var.m)
            moddict[var.m] = IndexedVector(Float64,var.m.numCols)
        end
    end

    # Collect like terms
    for ind in 1:length(a.vars)
        addelt!(moddict[a.vars[ind].m], a.vars[ind].col, a.coeffs[ind])
    end

    elm = 1
    term_str = Array(UTF8String, 2*length(a.vars))
    # For each model
    for m in keys(moddict)
        indvec = moddict[m]
        # For each non-zero for this model
        for i in 1:indvec.nnz
            idx = indvec.nzidx[i]
            elt = indvec.elts[idx]
            abs(elt) < PRINT_ZERO_TOL && continue  # e.g. x - x

            pre = abs(abs(elt)-1) < PRINT_ZERO_TOL ? "" : str_round(abs(elt)) * " "
            var = var_str(mode,m,idx)

            term_str[2*elm-1] = elt < 0 ? " - " : " + "
            term_str[2*elm  ] = "$pre$var"
            elm += 1
        end
    end

    if elm == 1
        # Will happen with cancellation of all terms
        # We should just return the constant, if its desired
        return show_constant ? str_round(a.constant) : "0"
    else
        # Correction for very first term - don't want a " + "/" - "
        term_str[1] = (term_str[1] == " - ") ? "-" : ""
        ret = join(term_str[1:2*(elm-1)])
        if abs(a.constant) >= PRINT_ZERO_TOL && show_constant
            ret = string(ret, a.constant < 0 ? " - " : " + ", str_round(abs(a.constant)))
        end
        return ret
    end
end

# Backwards compatability shim
affToStr(a::AffExpr) = aff_str(REPLMode,a)


#------------------------------------------------------------------------
## GenericQuadExpr
#------------------------------------------------------------------------
Base.print(io::IO, q::GenericQuadExpr) = print(io, quad_str(REPLMode,q))
Base.show( io::IO, q::GenericQuadExpr) = print(io, quad_str(REPLMode,q))
Base.writemime(io::IO, ::MIME"text/latex", q::GenericQuadExpr) =
    print(io, quad_str(IJuliaMode,q,mathmode=false))
# Generic string converter, called by mode-specific handlers
function quad_str(mode, q::GenericQuadExpr, times::String, sq::String)
    length(q.qvars1) == 0 && return aff_str(mode,q.aff)

    # Canonicalize x_i * x_j so i <= j
    for ind in 1:length(q.qvars1)
        if q.qvars2[ind].col < q.qvars1[ind].col
            q.qvars1[ind],q.qvars2[ind] = q.qvars2[ind],q.qvars1[ind]
        end
    end
    # Merge duplicates
    Q = sparse([v.col for v in q.qvars1], [v.col for v in q.qvars2], q.qcoeffs)
    I,J,V = findnz(Q)
    Qnnz = length(V)

    # Odd terms are +/i, even terms are the variables/coeffs
    term_str = Array(UTF8String, 2*Qnnz)
    if Qnnz > 0
        for ind in 1:Qnnz
            val = abs(V[ind])
            pre = (val == 1.0 ? "" : str_round(val)*" ")

            x = var_str(mode,q.qvars1[ind].m,I[ind])
            y = var_str(mode,q.qvars1[ind].m,J[ind])

            term_str[2*ind-1] = V[ind] < 0 ? " - " : " + "
            term_str[2*ind  ] = "$pre$x" * (x == y ? sq : "$times$y")
        end
        # Correction for first term as there is no space
        # between - and variable coefficient/name
        term_str[1] = V[1] < 0 ? "-" : ""
    end
    ret = join(term_str)

    if q.aff.constant == 0 && length(q.aff.vars) == 0
        return ret
    else
        aff = aff_str(mode,q.aff)
        if aff[1] == '-'
            return string(ret, " - ", aff[2:end])
        else
            return string(ret, " + ", aff)
        end
    end
end

# Backwards compatability shim
quadToStr(q::GenericQuadExpr) = quad_str(REPLMode,q)
# Handlers to use correct symbols
quad_str(::Type{REPLMode}, q::GenericQuadExpr) =
    quad_str(REPLMode, q, repl_times, repl_sq)
quad_str(::Type{IJuliaMode}, q::GenericQuadExpr; mathmode=true) =
    math(quad_str(IJuliaMode, q, ijulia_times, ijulia_sq), mathmode)


#------------------------------------------------------------------------
## GenericRangeConstraint
#------------------------------------------------------------------------
Base.print(io::IO, c::GenericRangeConstraint) = print(io, con_str(REPLMode,c))
Base.show( io::IO, c::GenericRangeConstraint) = print(io, con_str(REPLMode,c))
Base.writemime(io::IO, ::MIME"text/latex", c::GenericRangeConstraint) =
    print(io, con_str(IJuliaMode,c,mathmode=false))
# Generic string converter, called by mode-specific handlers
function con_str(mode, c::GenericRangeConstraint, leq, eq, geq)
    s = sense(c)
    a = aff_str(mode,c.terms,show_constant=false)
    if s == :range
        out_str = "$(str_round(c.lb)) $leq $a $leq $(str_round(c.ub))"
    else
        rel = s == :<= ? leq : (s == :>= ? geq : eq)
        out_str = string(a," ",rel," ",str_round(rhs(c)))
    end
    out_str
end
# Backwards compatability shim
conToStr(c::GenericRangeConstraint) = con_str(REPLMode,c)
# Handlers to use correct symbols
con_str(::Type{REPLMode}, c::GenericRangeConstraint; args...) =
    con_str(REPLMode, c, repl_leq, repl_eq, repl_geq)
con_str(::Type{IJuliaMode}, c::GenericRangeConstraint; mathmode=true) =
    math(con_str(IJuliaMode, c, ijulia_leq, ijulia_eq, ijulia_geq), mathmode)


#------------------------------------------------------------------------
## QuadConstraint
#------------------------------------------------------------------------
Base.print(io::IO, c::QuadConstraint) = print(io, con_str(REPLMode,c))
Base.show( io::IO, c::QuadConstraint) = print(io, con_str(REPLMode,c))
Base.writemime(io::IO, ::MIME"text/latex", c::QuadConstraint) =
    print(io, con_str(IJuliaMode,c,mathmode=false))
# Generic string converter, called by mode-specific handlers
function con_str(mode, c::QuadConstraint, leq, eq, geq)
    s = c.sense
    r = (s == :<=) ? leq : (s == :>= ? geq : eq)
    "$(quad_str(mode,c.terms)) $r 0"
end
# Backwards compatability shim
conToStr(c::QuadConstraint) = con_str(REPLMode,c)
# Handlers to use correct symbols
con_str(::Type{REPLMode}, c::QuadConstraint; args...) =
    con_str(REPLMode, c, repl_leq, repl_eq, repl_geq)
con_str(::Type{IJuliaMode}, c::QuadConstraint; mathmode=true) =
    math(con_str(IJuliaMode, c, ijulia_leq, ijulia_eq, ijulia_geq), mathmode)


#------------------------------------------------------------------------
## SOSConstraint
#------------------------------------------------------------------------
Base.print(io::IO, c::SOSConstraint) = print(io, con_str(REPLMode,c))
Base.show( io::IO, c::SOSConstraint) = print(io, con_str(REPLMode,c))
Base.writemime(io::IO, ::MIME"text/latex", c::SOSConstraint) =
    print(io, con_str(IJuliaMode,c,mathmode=false))
# Generic string converter, called by mode-specific handlers
function con_str(mode, c::SOSConstraint, open_set, close_set)
    term_str = [string(str_round(c.weights[i]), " ", c.terms[i])
                    for i in 1:length(c.terms)]
    "$(c.sostype): $open_set$(join(term_str,", "))$close_set"
end
# Handlers to use correct symbols
con_str(::Type{REPLMode}, c::SOSConstraint; args...) =
    con_str(REPLMode, c, repl_open_set, repl_close_set)
con_str(::Type{IJuliaMode}, c::SOSConstraint; mathmode=true) =
    math(con_str(IJuliaMode, c, ijulia_open_set, ijulia_close_set), mathmode)


#------------------------------------------------------------------------
## ConstraintRef
#------------------------------------------------------------------------
Base.print(io::IO, c::ConstraintRef{LinearConstraint}) = print(io, con_str(REPLMode,c.m.linconstr[c.idx]))
Base.print(io::IO, c::ConstraintRef{QuadConstraint})   = print(io, con_str(REPLMode,c.m.quadconstr[c.idx]))
Base.print(io::IO, c::ConstraintRef{SOSConstraint})    = print(io, con_str(REPLMode,c.m.sosconstr[c.idx]))
Base.show( io::IO, c::ConstraintRef{LinearConstraint}) = print(io, con_str(REPLMode,c.m.linconstr[c.idx]))
Base.show( io::IO, c::ConstraintRef{QuadConstraint})   = print(io, con_str(REPLMode,c.m.quadconstr[c.idx]))
Base.show( io::IO, c::ConstraintRef{SOSConstraint})    = print(io, con_str(REPLMode,c.m.sosconstr[c.idx]))
Base.writemime(io::IO, ::MIME"text/latex", c::ConstraintRef{LinearConstraint}) =
    print(io, con_str(IJuliaMode,c.m.linconstr[c.idx],mathmode=false))
Base.writemime(io::IO, ::MIME"text/latex", c::ConstraintRef{QuadConstraint}) =
    print(io, con_str(IJuliaMode,c.m.quadconstr[c.idx],mathmode=false))
Base.writemime(io::IO, ::MIME"text/latex", c::ConstraintRef{SOSConstraint}) =
    print(io, con_str(IJuliaMode,c.m.sosconstr[c.idx],mathmode=false))

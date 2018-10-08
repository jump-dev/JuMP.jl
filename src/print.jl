#  Copyright 2017, Iain Dunning, Joey Huchette, Miles Lubin, and contributors
#  This Source Code Form is subject to the terms of the Mozilla Public
#  License, v. 2.0. If a copy of the MPL was not distributed with this
#  file, You can obtain one at http://mozilla.org/MPL/2.0/.
#############################################################################
# JuMP
# An algebraic modeling language for Julia
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
# - Base.print and Base.string both delegate to Base.show, if they are not
#   separately defined.
#############################################################################

# Used for dispatching
abstract type PrintMode end
abstract type REPLMode <: PrintMode end
abstract type IJuliaMode <: PrintMode end

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

# TODO: get rid of this! This is only a helper, and should be Base.values
# (and probably live there, as well)
_values(x::Array) = x
_values(x) = Base.values(x)

# REPL-specific symbols
# Anything here: https://en.wikipedia.org/wiki/Windows-1252
# should probably work fine on Windows
const repl = Dict{Symbol,String}(
    :leq        => (Sys.iswindows() ? "<=" : "≤"),
    :geq        => (Sys.iswindows() ? ">=" : "≥"),
    :eq         => (Sys.iswindows() ? "==" : "="),
    :times      => "*",
    :sq         => "²",
    :ind_open   => "[",
    :ind_close  => "]",
    :for_all    => (Sys.iswindows() ? "for all" : "∀"),
    :in         => (Sys.iswindows() ? "in" : "∈"),
    :open_set   => "{",
    :dots       => (Sys.iswindows() ? ".." : "…"),
    :close_set  => "}",
    :union      => (Sys.iswindows() ? "or" : "∪"),
    :infty      => (Sys.iswindows() ? "Inf" : "∞"),
    :open_rng   => "[",
    :close_rng  => "]",
    :integer    => "integer",
    :succeq0    => " is semidefinite",
    :Vert       => (Sys.iswindows() ? "||" : "‖"),
    :sub2       => (Sys.iswindows() ? "_2" : "₂"))

# IJulia-specific symbols
const ijulia = Dict{Symbol,String}(
    :leq        => "\\leq",
    :geq        => "\\geq",
    :eq         => "=",
    :times      => "\\times ",
    :sq         => "^2",
    :ind_open   => "_{",
    :ind_close  => "}",
    :for_all    => "\\quad\\forall",
    :in         => "\\in",
    :open_set   => "\\{",
    :dots       => "\\dots",
    :close_set  => "\\}",
    :union      => "\\cup",
    :infty      => "\\infty",
    :open_rng   => "\\[",
    :close_rng  => "\\]",
    :integer    => "\\in \\mathbb{Z}",
    :succeq0    => "\\succeq 0",
    :Vert       => "\\Vert",
    :sub2       => "_2")

const PrintSymbols = Dict{Symbol,String}

# If not already mathmode, then wrap in MathJax start/close tags
math(s,mathmode) = mathmode ? s : "\$\$ $s \$\$"

# helper to look up corresponding JuMPContainerData
printdata(v::JuMPContainer) = _getmodel(v).varData[v]
getname(x::JuMPContainer) = hasmeta(x, :model) ? printdata(x).name : "__anon__"
function printdata(v::Array{Variable})
    if isempty(v)
        error("Cannot locate printing data for an empty array")
    end
    m = first(v).m
    m.varData[v]
end

#------------------------------------------------------------------------
## Model
#------------------------------------------------------------------------
function Base.print(io::IO, m::Model; ignore_print_hook=(m.printhook==nothing))
    ignore_print_hook || return m.printhook(io, m)
    print(io, model_str(REPLMode,m))
end

function Base.show(io::IO, m::Model)
    plural(n) = (n==1 ? "" : "s")
    print(io, m.objSense == :Max ? "Maximization" : ((m.objSense == :Min && (!isempty(m.obj) || (m.nlpdata !== nothing && isa(m.nlpdata.nlobj, NonlinearExprData)))) ? "Minimization" : "Feasibility"))
    println(io, " problem with:")
    nlin = length(m.linconstr)
    println(io, " * $(nlin) linear constraint$(plural(nlin))")
    nquad = length(m.quadconstr)
    if nquad > 0
        println(io, " * $(nquad) quadratic constraint$(plural(nquad))")
    end
    nsos = length(m.sosconstr)
    if nsos > 0
        println(io, " * $(nsos) SOS constraint$(plural(nsos))")
    end
    nsoc = length(m.socconstr)
    if nsoc > 0
        println(io, " * $(nsoc) SOC constraint$(plural(nsoc))")
    end
    nsdp = length(m.sdpconstr)
    if nsdp > 0
        println(io, " * $(nsdp) semidefinite constraint$(plural(nsdp))")
    end
    nlp = m.nlpdata
    if nlp !== nothing && length(nlp.nlconstr) > 0
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
    print(io, "Solver is ")
    if isa(m.solver, UnsetSolver)
        print(io, "default solver")
    elseif VERSION < v"0.7-"
        print(io, split(split(string(m.solver), "Solver")[1], ".")[2])
    else
        print(io, split(string(m.solver), "Solver")[1])
    end
end
Base.show(io::IO, ::MIME"text/latex", m::Model) =
    print(io, model_str(IJuliaMode,m))
function model_str(mode, m::Model, sym::PrintSymbols)
    ijl = mode == IJuliaMode
    sep = ijl ? " & " : " "
    eol = ijl ? "\\\\\n" : "\n"
    nlp = m.nlpdata

    # Objective
    qobj_str = quad_str(mode, m.obj)
    obj_sense = ijl ? (m.objSense == :Max ? "\\max" : "\\min")*"\\quad" :
                      (m.objSense == :Max ? "Max" : "Min")
    str = obj_sense * sep
    if nlp !== nothing && nlp.nlobj !== nothing
        str *= (qobj_str=="0" ? "" : "$qobj_str + ") * expr_str(m, mode, nlp.nlobj)
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
    for c in m.socconstr
        str *= sep * con_str(mode,c,mathmode=true) * eol
    end
    if nlp !== nothing
        for c in nlp.nlconstr
            str *= sep * con_str(m,mode,c,mathmode=true) * eol
        end
    end

    # Display indexed variables
    in_dictlist = falses(m.numCols)
    for d in m.dictList
        # make sure that you haven't changed a variable type in the collection
        allsame = true
        if !isempty(d)
            firstval = first(_values(d))
            cat = getcategory(firstval)
            lb, ub = getlowerbound(firstval), getupperbound(firstval)
            for v in _values(d)
                if !(getcategory(v) == cat && getlowerbound(v) == lb && getupperbound(v) == ub)
                    allsame = false
                    break
                elseif v in m.customNames
                    allsame = false
                    break
                end
            end
        end
        if allsame
            for it in _values(d)  # Mark variables in JuMPContainer as printed
                in_dictlist[it.col] = true
            end
            str *= sep * cont_str(mode,d,mathmode=true)  * eol
        end
    end

    # Display non-indexed variables
    for i in 1:m.numCols
        in_dictlist[i] && continue
        var_name = var_str(mode,m,i)
        var_lb, var_ub = m.colLower[i], m.colUpper[i]
        str_lb = var_lb == -Inf ? "-"*sym[:infty] : str_round(var_lb)
        str_ub = var_ub == +Inf ?     sym[:infty] : str_round(var_ub)
        var_cat = m.colCat[i]
        if var_cat == :Bin  # x binary
            str *= string(sep, var_name,
                            " ", sym[:in],
                            " ", sym[:open_set],
                            "0,1", sym[:close_set])
        elseif var_cat == :SemiInt  # x in union of 0 and {lb,...,ub}
            str *= string(sep, var_name,
                            " ", sym[:in],
                            " ", sym[:open_set],
                            str_lb, ",", sym[:dots], ",", str_ub,
                            sym[:close_set],
                            " ", sym[:union], " ",
                            sym[:open_set], "0", sym[:close_set])
        elseif var_cat == :SemiCont  # x in union of 0 and [lb,ub]
            str *= string(sep, var_name,
                            " ", sym[:in],
                            " ", sym[:open_rng],
                            str_lb, ",", str_ub,
                            sym[:close_rng],
                            " ", sym[:union], " ",
                            sym[:open_set], "0", sym[:close_set])
        elseif var_cat == :Fixed
            str *= string(sep, var_name, " = ", str_lb)
        elseif var_lb == -Inf && var_ub == +Inf # Free variable
            str *= string(sep, var_name)
        elseif var_lb == -Inf  # No lower bound
            str *= string(sep, var_name, " ", sym[:leq], " ", str_ub)
        elseif var_ub == +Inf  # No upper bound
            str *= string(sep, var_name, " ", sym[:geq], " ", str_lb)
        else
            str *= string(sep, str_lb, " ", sym[:leq],
                            " ", var_name, " ",
                            sym[:leq], " ", str_ub)
        end
        if var_cat == :Int
            str *= string(", ", sym[:integer])
        end
        str *= eol
    end

    if nlp !== nothing
        for i in 1:length(nlp.nlexpr)
            if ijl
                str *= "subexpression_{$i} = \\quad &" * expr_str(m,mode,nlp.nlexpr[i]) * eol
            else
                str *= "subexpression[$i]: " * expr_str(m,mode,nlp.nlexpr[i]) * eol
            end
        end
    end

    ijl ? "\$\$ \\begin{alignat*}{1}"*str*"\\end{alignat*}\n \$\$" :
          str
end

# Handlers to use correct symbols
model_str(::Type{REPLMode}, m::Model) =
    model_str(REPLMode, m, repl)
model_str(::Type{IJuliaMode}, m::Model; mathmode=true) =
    math(model_str(IJuliaMode, m, ijulia), mathmode)


#------------------------------------------------------------------------
## Variable
#------------------------------------------------------------------------
Base.show(io::IO, v::Variable) = print(io, var_str(REPLMode,v))
Base.show(io::IO, ::MIME"text/latex", v::Variable) =
    print(io, var_str(IJuliaMode,v,mathmode=false))
function var_str(mode, m::Model, col::Int; mathmode=true)
    colNames = mode == REPLMode ? m.colNames : m.colNamesIJulia
    if colNames[col] === EMPTYSTRING
        for cont in m.dictList
            fill_var_names(mode, colNames, cont)
        end
    end
    return math(colNames[col] == "" ? "col_$col" : colNames[col], mathmode)
end
function fill_var_names(mode, colNames, v::JuMPArray{Variable,N}) where N
    data = printdata(v)
    idxsets = data.indexsets
    lengths = map(length, idxsets)
    name = data.name
    cprod = cumprod([lengths...])
    for (ind,var) in enumerate(v.innerArray)
        idx_strs = [string( idxsets[1][mod1(ind,lengths[1])] )]
        for i = 2:N
            push!(idx_strs, string(idxsets[i][Int(ceil(mod1(ind,cprod[i]) / cprod[i-1]))]))
        end
        if mode == IJuliaMode
            colNames[var.col] = string(name, "_{", join(idx_strs,",") , "}")
        else
            colNames[var.col] = string(name,  "[", join(idx_strs,",") , "]")
        end
    end
end
function fill_var_names(mode, colNames, v::JuMPDict{Variable})
    name = printdata(v).name
    for (ind,var) in zip(keys(v),values(v))
        if mode == IJuliaMode
            colNames[var.col] = string(name, "_{", join([string(i) for i in ind],","), "}")
        else
            colNames[var.col] = string(name,  "[", join([string(i) for i in ind],","), "]")
        end
    end
end
function fill_var_names(mode, colNames, v::Array{Variable})
    isempty(v) && return
    sizes = size(v)
    m = first(v).m
    if !haskey(m.varData, v)
        return
    end
    name = m.varData[v].name
    for (ii,var) in enumerate(v)
        @assert var.m === m
        ind = _ind2sub(sizes, ii)
        colNames[var.col] = if mode === IJuliaMode
            string(name, "_{", join(ind, ","), "}")
        else
            string(name,  "[", join(ind, ","), "]")
        end
    end
    return
end

# Handlers to use correct symbols
var_str(::Type{REPLMode}, v::Variable) =
    var_str(REPLMode, v.m, v.col)
var_str(::Type{IJuliaMode}, v::Variable; mathmode=true) =
    var_str(IJuliaMode, v.m, v.col, mathmode=mathmode)

#------------------------------------------------------------------------
## Norm
#------------------------------------------------------------------------
Base.show(io::IO, j::Norm) = print(io, norm_str(REPLMode,j))
Base.show(io::IO, ::MIME"text/latex", j::Norm) =
    print(io, norm_str(IJuliaMode,j))

function norm_str(mode, n::Norm, sym::PrintSymbols)
    string(sym[:Vert], "[",
            join(map(t->aff_str(mode,t),n.terms),","),
            "]", sym[:Vert], sym[:sub2])
end

# Handlers to use correct symbols
norm_str(::Type{REPLMode}, n::Norm) =
    norm_str(REPLMode, n, repl)
norm_str(::Type{IJuliaMode}, n::Norm; mathmode=true) =
    math(norm_str(IJuliaMode, n, ijulia), mathmode)

exprstr(n::Norm) = norm_str(REPLMode, n)

#------------------------------------------------------------------------
## JuMPContainer{Variable}
#------------------------------------------------------------------------
Base.show(io::IO, j::JuMPContainer{Variable}) = print(io, cont_str(REPLMode,j))
Base.show(io::IO, j::Array{Variable}) = print(io, cont_str(REPLMode,j))
Base.show(io::IO, ::MIME"text/latex", j::Union{JuMPContainer{Variable},Array{Variable}}) =
    print(io, cont_str(IJuliaMode,j,mathmode=false))
# Generic string converter, called by mode-specific handlers

# Assumes that !isempty(j)
_getmodel(j::Array{Variable}) = first(j).m
_getmodel(j::JuMPContainer) = getmeta(j, :model)

# 0.7 compat
const _print_array = if VERSION < v"0.7-"
    Base.showarray
else
    Base.print_array
end

function cont_str(mode, j, sym::PrintSymbols)
    # Check if anything in the container
    if isempty(j)
        name = isa(j, JuMPContainer) ? printdata(j).name : "Empty Array{Variable}"
        return "$name (no indices)"
    end

    m = _getmodel(j)

    # If this looks like a user-created Array, then defer to base printing
    if !haskey(m.varData, j)
        @assert isa(j, Array{Variable})
        if ndims(j) == 1
            return sprint((io,v) -> Base.show_vector(io, v, "[", "]"), j)
        else
            return sprint((io,X) -> _print_array(io, X), j)
        end
    end

    data = printdata(j)

    # 1. construct the part with variable name and indexing
    locvars = map(data.indexexprs) do tmp
        var = tmp.idxvar
        if var == nothing
            return ""
        else
            return string(var)
        end
    end
    num_dims = length(data.indexsets)
    idxvars = Array{String}(undef, num_dims)
    dimidx = 1
    for i in 1:num_dims
        if data.indexexprs[i].idxvar == nothing
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
    name_idx = string(data.name, sym[:ind_open], join(idxvars,","), sym[:ind_close])
    # 2. construct part with what we index over
    idx_sets = sym[:for_all]*" "*join(map(dim->string(idxvars[dim], " ", sym[:in],
                                " ", sym[:open_set],
                                cont_str_set(data.indexsets[dim],sym[:dots]),
                                sym[:close_set]), 1:num_dims), ", ")
    # 3. Handle any conditions
    if isa(j, JuMPDict) && data.condition != :()
       idx_sets *= string(" s.t. ",join(parse_conditions(data.condition), " and "))
    end

    # 4. Bounds and category, if possible, and return final string
    a_var = first(_values(j))
    model = a_var.m
    var_cat = model.colCat[a_var.col]
    var_lb  = model.colLower[a_var.col]
    var_ub  = model.colUpper[a_var.col]
    # Variables may have different bounds, so we can't really print nicely
    # at this time (possibly ever, as they could have been changed post
    # creation, which we'd never be able to handle.
    all_same_lb = true
    all_same_ub = true
    for var in _values(j)
        all_same_lb &= model.colLower[var.col] == var_lb
        all_same_ub &= model.colUpper[var.col] == var_ub
    end
    str_lb = var_lb == -Inf ? "-"*sym[:infty] : str_round(var_lb)
    str_ub = var_ub == +Inf ?     sym[:infty] : str_round(var_ub)

    # Special case bounds printing based on the category
    if var_cat == :Bin  # x in {0,1}
        return "$name_idx $(sym[:in]) $(sym[:open_set])0,1$(sym[:close_set]) $idx_sets"
    elseif var_cat == :SemiInt  # x in union of 0 and {lb,...,ub}
        si_lb = all_same_lb ? str_lb : sym[:dots]
        si_ub = all_same_ub ? str_ub : sym[:dots]
        return "$name_idx $(sym[:in]) $(sym[:open_set])$si_lb,$(sym[:dots]),$si_ub$(sym[:close_set]) $(sym[:union]) $(sym[:open_set])0$(sym[:close_set]) $idx_sets"
    elseif var_cat == :SemiCont  # x in union of 0 and [lb,ub]
        si_lb = all_same_lb ? str_lb : sym[:dots]
        si_ub = all_same_ub ? str_ub : sym[:dots]
        return "$name_idx $(sym[:in]) $(sym[:open_rng])$si_lb,$si_ub$(sym[:close_rng]) $(sym[:union]) $(sym[:open_set])0$(sym[:close_set]) $idx_sets"
    elseif var_cat == :Fixed
        si_bnd = all_same_lb ? str_lb : sym[:dots]
        return "$name_idx = $si_bnd $idx_sets"
    end
    # Continuous and Integer
    idx_sets = var_cat == :Int ? ", $(sym[:integer]), $idx_sets" : " $idx_sets"
    if all_same_lb && all_same_ub
        # Free variable
        var_lb == -Inf && var_ub == +Inf && return "$name_idx$idx_sets"
        # No lower bound
        var_lb == -Inf && return "$name_idx $(sym[:leq]) $str_ub$idx_sets"
        # No upper bound
        var_ub == +Inf && return "$name_idx $(sym[:geq]) $str_lb$idx_sets"
        # Range
        return "$str_lb $(sym[:leq]) $name_idx $(sym[:leq]) $str_ub$idx_sets"
    end
    if all_same_lb && !all_same_ub
        var_lb == -Inf && return "$name_idx $(sym[:leq]) $(sym[:dots])$idx_sets"
        return "$str_lb $(sym[:leq]) $name_idx $(sym[:leq]) $(sym[:dots])$idx_sets"
    end
    if !all_same_lb && all_same_ub
        var_ub == +Inf && return "$name_idx $(sym[:geq]) $(sym[:dots])$idx_sets"
        return "$(sym[:dots]) $(sym[:leq]) $name_idx $(sym[:leq]) $str_ub$idx_sets"
    end
    return "$(sym[:dots]) $(sym[:leq]) $name_idx $(sym[:leq]) $(sym[:dots])$idx_sets"
end

# UTILITY FUNCTIONS FOR cont_str
function cont_str_set(idxset::Union{AbstractRange,Array}, dots)  # 2:2:20 -> {2,4..18,20}
    length(idxset) == 0 && return dots
    length(idxset) == 1 && return string(idxset[1])
    length(idxset) == 2 && return string(idxset[1],",",idxset[2])
    length(idxset) == 3 && return string(idxset[1],",",idxset[2],",",idxset[3])
    length(idxset) == 4 && return string(idxset[1],",",idxset[2],",",idxset[3],",",idxset[4])
    length(idxset) == 5 && return string(idxset[1],",",idxset[2],",",idxset[3],",",idxset[4],",",idxset[5])
    return string(idxset[1],",",idxset[2],",",dots,",",idxset[end-1],",",idxset[end])
end
cont_str_set(idxset, dots) = return dots # Fallback
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
cont_str(::Type{REPLMode}, j; mathmode=false) =
    cont_str(REPLMode, j, repl)
cont_str(::Type{IJuliaMode}, j; mathmode=true) =
    math(cont_str(IJuliaMode, j, ijulia), mathmode)

#------------------------------------------------------------------------
## JuMPContainer{Float64}
#------------------------------------------------------------------------
Base.show(io::IO, j::JuMPContainer{Float64}) = print(io, val_str(REPLMode,j))
function val_str(mode, j::JuMPArray{Float64,N}) where N
    out_str = "$(getname(j)): $N dimensions:"
    if isempty(j)
        return out_str * "\n  (no entries)"
    end

    function val_str_rec(depth, parent_index::Vector{Any}, parent_str::AbstractString)
        # Turn index set into strings
        indexset = j.indexsets[depth]
        isa(indexset, AbstractArray) || (indexset = collect(indexset))
        index_strs = map(string, indexset)

        # Determine longest index so we can align columns
        max_index_len = 0
        for index_str in index_strs
            max_index_len = max(max_index_len, textwidth(index_str))
        end

        # If have recursed, we need to prepend the parent's index strings
        # accumulated, as well as white space so the alignment works.
        for i = 1:length(index_strs)
            index_strs[i] = parent_str * lpad(index_strs[i],max_index_len," ")
        end

        # Create a string for the number of spaces we need to indent
        indent = " "^(2*(depth-1))

        # Determine the need to recurse
        if depth == N
            # Deepest level
            for i = 1:length(indexset)
                value = length(parent_index) == 0 ?
                            j[indexset[i]] :
                            j[parent_index...,indexset[i]]
                out_str *= "\n" * indent * "[" * index_strs[i] * "] = $value"
            end
        else
            # At least one more layer to go
            for i = 1:length(indexset)
                index = indexset[i]
                # Print the ":" version of indices we will recurse over
                out_str *= "\n" * indent * "[" * index_strs[i] * ",:"^(N-depth) * "]"
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
function val_str(mode, dict::JuMPDict{Float64,N}) where N
    nelem = length(dict.tupledict)
    isempty(dict) && return ""
    out_str  = "$(getname(dict)): $N dimensions, $nelem "
    out_str *= nelem == 1 ? "entry" : "entries"
    out_str *= ":"

    sortedkeys = sort(collect(keys(dict.tupledict)), lt = _isless)

    ndim = length(first(keys(dict.tupledict)))

    key_strs = Array{String}(undef, length(dict), ndim)
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
Base.show(io::IO, a::AffExpr) = print(io, aff_str(REPLMode,a))
Base.show(io::IO, ::MIME"text/latex", a::AffExpr) =
    print(io, math(aff_str(IJuliaMode,a),false))
# Generic string converter, called by mode-specific handlers
function aff_str(mode, a::AffExpr, show_constant=true)
    # If the expression is empty, return the constant (or 0)
    if length(a.vars) == 0
        return show_constant ? str_round(a.constant) : "0"
    end

    # Get reference to models included in this expression
    moddict = Dict{Model,IndexedVector{Float64}}()
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
    term_str = Array{String}(undef, 2*length(a.vars))
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
# Precompile for faster boot times
Base.precompile(aff_str, (Type{JuMP.REPLMode}, AffExpr, Bool))
Base.precompile(aff_str, (Type{JuMP.IJuliaMode}, AffExpr, Bool))
Base.precompile(aff_str, (Type{JuMP.REPLMode}, AffExpr))
Base.precompile(aff_str, (Type{JuMP.IJuliaMode}, AffExpr))


#------------------------------------------------------------------------
## GenericQuadExpr
#------------------------------------------------------------------------
Base.show(io::IO, q::GenericQuadExpr) = print(io, quad_str(REPLMode,q))
Base.show(io::IO, ::MIME"text/latex", q::GenericQuadExpr) =
    print(io, quad_str(IJuliaMode,q,mathmode=false))
# Generic string converter, called by mode-specific handlers
function quad_str(mode, q::GenericQuadExpr, sym)
    length(q.qvars1) == 0 && return aff_str(mode,q.aff)

    # Canonicalize x_i * x_j so i <= j
    for ind in 1:length(q.qvars1)
        if q.qvars2[ind].col < q.qvars1[ind].col
            q.qvars1[ind],q.qvars2[ind] = q.qvars2[ind],q.qvars1[ind]
        end
    end
    # Merge duplicates
    Q = dropzeros(sparse([v.col for v in q.qvars1], [v.col for v in q.qvars2], q.qcoeffs))
    I,J,V = findnz(Q)
    Qnnz = length(V)

    # Odd terms are +/i, even terms are the variables/coeffs
    term_str = Array{String}(undef, 2*Qnnz)
    if Qnnz > 0
        for ind in 1:Qnnz
            val = abs(V[ind])
            pre = (val == 1.0 ? "" : str_round(val)*" ")

            x = var_str(mode,q.qvars1[ind].m,I[ind])
            y = var_str(mode,q.qvars1[ind].m,J[ind])

            term_str[2*ind-1] = V[ind] < 0 ? " - " : " + "
            term_str[2*ind  ] = "$pre$x" * (x == y ? sym[:sq] : "$(sym[:times])$y")
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
        elseif length(ret) == 0
            return string(aff)
        else
            return string(ret, " + ", aff)
        end
    end
end

# Handlers to use correct symbols
quad_str(::Type{REPLMode}, q::GenericQuadExpr) =
    quad_str(REPLMode, q, repl)
quad_str(::Type{IJuliaMode}, q::GenericQuadExpr; mathmode=true) =
    math(quad_str(IJuliaMode, q, ijulia), mathmode)

#------------------------------------------------------------------------
## SOCExpr
#------------------------------------------------------------------------
Base.show(io::IO, c::SOCExpr) = print(io, expr_str(REPLMode, c))
Base.show(io::IO, ::MIME"text/latex", c::SOCExpr) =
    print(io, expr_str(IJuliaMode, c))
function expr_str(mode, c::SOCExpr)
    coeff = (c.coeff == 1) ? "" : string(c.coeff, " ")
    aff   = aff_str(mode, c.aff)
    if aff[1] == '-'
        chain = " - "
        aff = aff[2:end]
    elseif aff == "0"
        aff = ""
        chain = ""
    else  # positive
        chain = " + "
    end
    string(coeff, norm_str(mode, c.norm), chain, aff)
end

#------------------------------------------------------------------------
## NonlinearExprData
#------------------------------------------------------------------------
#Base.show(io::IO, c::NonlinearExprData) = print(io, expr_str(REPLMode, c))
#Base.show(io::IO, ::MIME"text/latex", c::NonlinearExprData) =
#    print(io, expr_str(IJuliaMode, c))
function expr_str(m::Model, mode, c::NonlinearExprData)
    return string(tapeToExpr(m, 1, c.nd, adjmat(c.nd), c.const_values, [], [], m.nlpdata.user_operators, false, false, mode))
end

#------------------------------------------------------------------------
## GenericRangeConstraint
#------------------------------------------------------------------------
Base.show(io::IO, c::GenericRangeConstraint) = print(io, con_str(REPLMode,c))
Base.show(io::IO, ::MIME"text/latex", c::GenericRangeConstraint) =
    print(io, con_str(IJuliaMode,c,mathmode=false))
# Generic string converter, called by mode-specific handlers
function con_str(mode, c::GenericRangeConstraint, sym)
    s = sense(c)
    a = aff_str(mode,c.terms,false)
    if s == :range
        out_str = "$(str_round(c.lb)) $(sym[:leq]) $a $(sym[:leq]) $(str_round(c.ub))"
    else
        rel = s == :<= ? sym[:leq] : (s == :>= ? sym[:geq] : sym[:eq])
        out_str = string(a," ",rel," ",str_round(rhs(c)))
    end
    out_str
end
# Handlers to use correct symbols
con_str(::Type{REPLMode}, c::GenericRangeConstraint; args...) =
    con_str(REPLMode, c, repl)
con_str(::Type{IJuliaMode}, c::GenericRangeConstraint; mathmode=true) =
    math(con_str(IJuliaMode, c, ijulia), mathmode)


#------------------------------------------------------------------------
## QuadConstraint
#------------------------------------------------------------------------
Base.show(io::IO, c::QuadConstraint) = print(io, con_str(REPLMode,c))
Base.show(io::IO, ::MIME"text/latex", c::QuadConstraint) =
    print(io, con_str(IJuliaMode,c,mathmode=false))
# Generic string converter, called by mode-specific handlers
function con_str(mode, c::QuadConstraint, sym)
    s = c.sense
    r = (s == :<=) ? sym[:leq] : (s == :>= ? sym[:geq] : sym[:eq])
    "$(quad_str(mode,c.terms)) $r 0"
end
# Handlers to use correct symbols
con_str(::Type{REPLMode}, c::QuadConstraint; args...) =
    con_str(REPLMode, c, repl)
con_str(::Type{IJuliaMode}, c::QuadConstraint; mathmode=true) =
    math(con_str(IJuliaMode, c, ijulia), mathmode)

#------------------------------------------------------------------------
## SOCConstraint
#------------------------------------------------------------------------
Base.show(io::IO, c::SOCConstraint) = print(io, con_str(REPLMode,c))
Base.show(io::IO, ::MIME"text/latex", c::SOCConstraint) =
    print(io, con_str(IJuliaMode,c))
function con_str(mode, c::SOCConstraint, sym::PrintSymbols)
    ne = c.normexpr
    coeff = ne.coeff == 1 ? "" : string(ne.coeff, " ")
    nrm   = norm_str(mode, ne.norm)
    aff   = aff_str(mode, -ne.aff)
    string(coeff, nrm, " $(repl[:leq]) ", aff)
end
# Handlers to use correct symbols
con_str(::Type{REPLMode}, c::SOCConstraint; args...) =
    con_str(REPLMode, c, repl)
con_str(::Type{IJuliaMode}, c::SOCConstraint; mathmode=true) =
    math(con_str(IJuliaMode, c, ijulia), mathmode)

#------------------------------------------------------------------------
## SOSConstraint
#------------------------------------------------------------------------
Base.show(io::IO, c::SOSConstraint) = print(io, con_str(REPLMode,c))
Base.show(io::IO, ::MIME"text/latex", c::SOSConstraint) =
    print(io, con_str(IJuliaMode,c,mathmode=false))
# Generic string converter, called by mode-specific handlers
function con_str(mode, c::SOSConstraint, sym::PrintSymbols)
    term_str = [string(str_round(c.weights[i]), " ", c.terms[i])
                    for i in 1:length(c.terms)]
    "$(c.sostype): $(sym[:open_set])$(join(term_str,", "))$(sym[:close_set])"
end
# Handlers to use correct symbols
con_str(::Type{REPLMode}, c::SOSConstraint; args...) =
    con_str(REPLMode, c, repl)
con_str(::Type{IJuliaMode}, c::SOSConstraint; mathmode=true) =
    math(con_str(IJuliaMode, c, ijulia), mathmode)

#------------------------------------------------------------------------
## SDConstraint
#------------------------------------------------------------------------
Base.show(io::IO, c::SDConstraint) = print(io, con_str(REPLMode,c))
Base.show(io::IO, ::MIME"text/latex", c::SDConstraint) =
    print(io, con_str(IJuliaMode,c,mathmode=false))
# Generic string converter, called by mode-specific handlers
function con_str(mode, c::SDConstraint, succeq0)
    t = c.terms
    str = sprint(show, MIME"text/plain"(), t)
    splitted = split(str, "\n")[2:end]
    center = ceil(Int, length(splitted)/2)
    splitted[center] *= succeq0
    join(splitted, "\n")
end
# Handlers to use correct symbols
con_str(::Type{REPLMode}, c::SDConstraint; args...) =
    con_str(REPLMode, c, repl[:succeq0])
con_str(::Type{IJuliaMode}, c::SDConstraint; mathmode=true) =
    math(con_str(IJuliaMode, c, ijulia[:succeq0], mathmode))

#------------------------------------------------------------------------
## NonlinearConstraint
#------------------------------------------------------------------------
Base.show(io::IO, c::NonlinearConstraint) = print(io, con_str(REPLMode,c))
Base.show(io::IO, ::MIME"text/latex", c::NonlinearConstraint) =
    print(io, con_str(IJuliaMode,c,mathmode=false))
# Generic string converter, called by mode-specific handlers
function con_str(m::Model, mode, c::NonlinearConstraint, sym)
    s = sense(c)
    nl = expr_str(m, mode, c.terms)
    if s == :range
        out_str = "$(str_round(c.lb)) $(sym[:leq]) $nl $(sym[:leq]) $(str_round(c.ub))"
    else
        rel = s == :<= ? sym[:leq] : (s == :>= ? sym[:geq] : sym[:eq])
        out_str = string(nl," ",rel," ",str_round(rhs(c)))
    end
    out_str
end
# Handlers to use correct symbols
con_str(m::Model, ::Type{REPLMode}, c::GenericRangeConstraint; args...) =
    con_str(m, REPLMode, c, repl)
con_str(m::Model, ::Type{IJuliaMode}, c::GenericRangeConstraint; mathmode=true) =
    math(con_str(m, IJuliaMode, c, ijulia), mathmode)
#------------------------------------------------------------------------
## ConstraintRef
#------------------------------------------------------------------------
Base.show(io::IO, c::ConstraintRef{Model,LinearConstraint}) = print(io, con_str(REPLMode,c.m.linconstr[c.idx]))
Base.show(io::IO, c::ConstraintRef{Model,QuadConstraint})   = print(io, con_str(REPLMode,c.m.quadconstr[c.idx]))
Base.show(io::IO, c::ConstraintRef{Model,SOSConstraint})    = print(io, con_str(REPLMode,c.m.sosconstr[c.idx]))
Base.show(io::IO, c::ConstraintRef{Model,SOCConstraint})    = print(io, con_str(REPLMode,c.m.socconstr[c.idx]))
Base.show(io::IO, c::ConstraintRef{Model,SDConstraint})     = print(io, con_str(REPLMode,c.m.sdpconstr[c.idx]))
Base.show(io::IO, c::ConstraintRef{Model,NonlinearConstraint}) = print(io, con_str(c.m, REPLMode, c.m.nlpdata.nlconstr[c.idx]))

Base.show(io::IO, ::MIME"text/latex", c::ConstraintRef{Model,LinearConstraint}) =
    print(io, con_str(IJuliaMode,c.m.linconstr[c.idx],mathmode=false))
Base.show(io::IO, ::MIME"text/latex", c::ConstraintRef{Model,QuadConstraint}) =
    print(io, con_str(IJuliaMode,c.m.quadconstr[c.idx],mathmode=false))
Base.show(io::IO, ::MIME"text/latex", c::ConstraintRef{Model,SOSConstraint}) =
    print(io, con_str(IJuliaMode,c.m.sosconstr[c.idx],mathmode=false))
Base.show(io::IO, ::MIME"text/latex", c::ConstraintRef{Model,NonlinearConstraint}) = print(io, con_str(c.m, IJuliaMode, c.m.nlpdata.nlconstr[c.idx],mathmode=false))

#------------------------------------------------------------------------
## Nonlinear expression/parameter reference
#------------------------------------------------------------------------
Base.show(io::IO, ex::NonlinearExpression) = Base.show(io, "Reference to nonlinear expression #$(ex.index)")
Base.show(io::IO, p::NonlinearParameter) = Base.show(io, "Reference to nonlinear parameter #$(p.index)")

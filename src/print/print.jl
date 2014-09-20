#############################################################################
# JuMP
# An algebraic modelling langauge for Julia
# See http://github.com/JuliaOpt/JuMP.jl
#############################################################################
# print/print.jl
# All "pretty printers" for JuMP types.
# - Delegates to appropriate methods for REPL or IJulia as appropriate.
# - Provides generic conversion-to-string code for both.
# - To find printing code for a type, search for `## TypeName`
# - Code here does not need to be fast, in fact simplicity trumps speed
#   within reason as this code is thorny enough as it is.
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
    str = string(f)
    length(str) >= 2 && str[end-1:end] == ".0" ? str[1:end-2] : str
end

include("ijulia.jl")
include("repl.jl")

#########################################################################
# VARIABLES
#########################################################################
#------------------------------------------------------------------------
## Variable
#------------------------------------------------------------------------
Base.print(io::IO, v::Variable) = print(io, getName(v))
Base.show( io::IO, v::Variable) = show( io, getName(v))
Base.writemime(io::IO, ::MIME"text/latex", v::Variable) = 
    print(io, getName(v))
#------------------------------------------------------------------------
## JuMPContainer{Variable}
#------------------------------------------------------------------------
Base.print(io::IO, j::JuMPContainer{Variable}) = print(io, cont_str(REPLMode,j))
Base.show( io::IO, j::JuMPContainer{Variable}) = show( io, cont_str(REPLMode,j))
Base.writemime(io::IO, ::MIME"text/latex", j::JuMPContainer{Variable}) =
    print(io, cont_str(IJuliaMode,j,mathmode=false))
function cont_str(mode, j::JuMPContainer{Variable}, leq, eq, geq,
                            ind_open, ind_close, for_all, in_set,
                            open_set, mid_set, close_set, union, infty,
                            open_rng, close_rng, integer)
    # Check if anything in the container
    isempty(j) && return string(j.name, " (no indices)")

    # 1. construct the part with variable name and indexing
    num_dims = length(j.indexsets)
    name_idx = string(j.name, ind_open, join(DIMS[1:num_dims],","), ind_close)
    # 2. construct part with what we index over
    idx_sets = for_all*" "*join(map(dim->string(DIMS[dim], " ", in_set, " ", open_set,
                                cont_str_set(j.indexsets[dim], mid_set),
                                close_set), 1:num_dims), ", ")
    # 3. Handle any conditionals
    #if isa(dict, JuMPDict) && !isempty(dict.condition)
    #    tail_str *= " s.t. $(join(parse_conditions(j.condition[1]), " and "))"
    #end

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
    end
    # Continuous and Integer
    idx_sets = var_cat == :Int ? ", $integer, $idx_sets" : " $idx_sets"
    if all_same_lb && all_same_ub
        # Free variable
        var_lb == -Inf && var_ub == +Inf && return "$name_idx$idx_sets"
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
#=parse_conditions(not_an_expr) = not_an_expr
function parse_conditions(expr::Expr)
    ret = {}
    if expr.head != :&&
        return {expr}
    end
    recurse = map(parse_conditions, expr.args)
    vcat(ret, recurse...)
end=#




#########################################################################
# EXPRESSIONS
#########################################################################
#------------------------------------------------------------------------
## AffExpr  (not GenericAffExpr)
#------------------------------------------------------------------------
Base.print(io::IO, a::AffExpr) = print(io, aff_str(REPLMode,a))
Base.show( io::IO, a::AffExpr) = show( io, aff_str(REPLMode,a))
Base.writemime(io::IO, ::MIME"text/latex", a::AffExpr) =
    print(io, aff_str(IJuliaMode,a,mathmode=false))
# Generic string converter, called by mode-specific handlers
# Doesn't need to be passed `mode` as it has no need to pass it on.
function aff_str(a::AffExpr; show_constant=true)
    # If the expression is empty, return the constant (or 0)
    if length(a.vars) == 0
        return show_constant ? str_round(a.constant) : "0"
    end

    # Get reference to models included in this expression
    moddict = Dict{Model,IndexedVector}()
    for var in a.vars
        if !haskey(moddict, var.m)
            checkNameStatus(var.m)
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
            var = getName(m,idx)

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
Base.show( io::IO, q::GenericQuadExpr) = show( io, quad_str(REPLMode,q))
Base.writemime(io::IO, ::MIME"text/latex", q::GenericQuadExpr) =
    print(io, quad_str(IJuliaMode,q,mathmode=false))
# Generic string converter, called by mode-specific handlers
function quad_str(mode, q::GenericQuadExpr, times::String, sq::String)
    length(q.qvars1) == 0 && return aff_str(mode,q.aff)

    # Check model ownership
    for ind in 1:length(q.qvars1)
        if q.qvars1[ind].m != q.qvars2[ind].m
            error("You cannot have a quadratic term with variables from different models.")
        end
    end

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

            x = getName(Variable(q.qvars1[ind].m,I[ind]))
            y = getName(Variable(q.qvars1[ind].m,J[ind]))
            
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


#########################################################################
# CONSTRAINTS
#########################################################################
#------------------------------------------------------------------------
## GenericRangeConstraint
#------------------------------------------------------------------------
Base.print(io::IO, c::GenericRangeConstraint) = print(io, con_str(REPLMode,c))
Base.show( io::IO, c::GenericRangeConstraint) = show( io, con_str(REPLMode,c))
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
        out_str = "$a $rel $(str_round(rhs(c)))"
    end
    out_str
end
# Backwards compatability shim
conToStr(c::GenericRangeConstraint) = con_str(REPLMode,c)

#------------------------------------------------------------------------
## QuadConstraint
#------------------------------------------------------------------------
Base.print(io::IO, c::QuadConstraint) = print(io, con_str(REPLMode,c))
Base.show( io::IO, c::QuadConstraint) = show( io, con_str(REPLMode,c))
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

#------------------------------------------------------------------------
## SOSConstraint
#------------------------------------------------------------------------
Base.print(io::IO, c::SOSConstraint) = print(io, con_str(REPLMode,c))
Base.show( io::IO, c::SOSConstraint) = show( io, con_str(REPLMode,c))
Base.writemime(io::IO, ::MIME"text/latex", c::SOSConstraint) =
    print(io, con_str(IJuliaMode,c,mathmode=false))
# Generic string converter, called by mode-specific handlers
function con_str(mode, c::SOSConstraint, open_set, close_set)
    term_str = [string(str_round(c.weights[i]), " ", c.terms[i])
                    for i in 1:length(c.terms)]
    "$(c.sostype): $open_set$(join(term_str,", "))$close_set"
end

#------------------------------------------------------------------------
## ConstraintRef
#------------------------------------------------------------------------
Base.print(io::IO, c::ConstraintRef{LinearConstraint}) = print(io, con_str(REPLMode,c.m.linconstr[c.idx]))
Base.print(io::IO, c::ConstraintRef{QuadConstraint})   = print(io, con_str(REPLMode,c.m.quadconstr[c.idx]))
Base.print(io::IO, c::ConstraintRef{SOSConstraint})    = print(io, con_str(REPLMode,c.m.sosconstr[c.idx]))
Base.show( io::IO, c::ConstraintRef{LinearConstraint}) = show( io, con_str(REPLMode,c.m.linconstr[c.idx]))
Base.show( io::IO, c::ConstraintRef{QuadConstraint})   = show( io, con_str(REPLMode,c.m.quadconstr[c.idx]))
Base.show( io::IO, c::ConstraintRef{SOSConstraint})    = show( io, con_str(REPLMode,c.m.sosconstr[c.idx]))
Base.writemime(io::IO, ::MIME"text/latex", c::ConstraintRef{LinearConstraint}) =
    print(io, con_str(IJuliaMode,c.m.linconstr[c.idx],mathmode=false))
Base.writemime(io::IO, ::MIME"text/latex", c::ConstraintRef{QuadConstraint}) =
    print(io, con_str(IJuliaMode,c.m.quadconstr[c.idx],mathmode=false))
Base.writemime(io::IO, ::MIME"text/latex", c::ConstraintRef{SOSConstraint}) =
    print(io, con_str(IJuliaMode,c.m.sosconstr[c.idx],mathmode=false))
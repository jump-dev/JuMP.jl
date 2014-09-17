#############################################################################
# JuMP
# An algebraic modelling langauge for Julia
# See http://github.com/JuliaOpt/JuMP.jl
#############################################################################
# print/print.jl
# - Delegates to appropriate methods for REPL or IJulia as appropriate.
# - Provides generic conversion-to-string code for both.
#############################################################################

abstract PrintMode
abstract REPLMode <: PrintMode
abstract IJuliaMode <: PrintMode

const print_zero_tol = 1e-10

include("ijulia.jl")
include("repl.jl")

#########################################################################
# EXPRESSIONS
#########################################################################
#------------------------------------------------------------------------
# AffExpr  (not GenericAffExpr)
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
            abs(elt) < print_zero_tol && continue  # e.g. x - x

            pre = abs(abs(elt)-1) < print_zero_tol ? "" : str_round(abs(elt)) * " "
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
        if abs(a.constant) >= print_zero_tol && show_constant
            ret = string(ret, a.constant < 0 ? " - " : " + ", str_round(abs(a.constant)))
        end
        return ret
    end
end

# Backwards compatability shim
affToStr(a::AffExpr) = aff_str(REPLMode,a)

#------------------------------------------------------------------------
# GenericQuadExpr
#------------------------------------------------------------------------
Base.print(io::IO, q::GenericQuadExpr) = print(io, quad_str(REPLMode,q))
Base.show(io::IO,  q::GenericQuadExpr) = show(io,  quad_str(REPLMode,q))
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
# GenericRangeConstraint
#------------------------------------------------------------------------
Base.print(io::IO, c::GenericRangeConstraint) = print(io, con_str(REPLMode,c))
Base.show(io::IO,  c::GenericRangeConstraint) = show(io,  con_str(REPLMode,c))
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
# QuadConstraint
#------------------------------------------------------------------------
Base.print(io::IO, c::QuadConstraint) = print(io, con_str(c))
Base.show(io::IO,  c::QuadConstraint) = show(io,  con_str(c))
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
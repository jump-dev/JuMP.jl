#############################################################################
# JuMP
# An algebraic modelling langauge for Julia
# See http://github.com/JuliaOpt/JuMP.jl
#############################################################################
# print/print.jl
# Delegates to appropriate methods for REPL or IJulia as appropriate
#############################################################################

abstract PrintMode
abstract REPLMode <: PrintMode
abstract IJuliaMode <: PrintMode

include("ijulia.jl")
include("repl.jl")

# TODO
aff_str(mode,a::AffExpr,show_constant=true) = affToStr(a,show_constant)

#########################################################################
# EXPRESSIONS
#########################################################################
#------------------------------------------------------------------------
# GenericQuadExpr
#------------------------------------------------------------------------
Base.print(io::IO, q::GenericQuadExpr) =
    print(io, quad_str(REPLMode,q))
Base.show(io::IO, q::GenericQuadExpr) =
    show(io, quad_str(REPLMode,q))
Base.writemime(io::IO, ::MIME"text/latex", q::GenericQuadExpr) =
    print(io, quad_str(IJuliaMode,q,mathmode=false))

function quad_str(mode, q::GenericQuadExpr, times, sq)
    length(q.qvars1) == 0 && return aff_str(mode,q.aff)

    # Check model ownership
    for ind in 1:length(q.qvars1)
        if q.qvars1[ind].m != q.qvars2[ind].m
            error("You cannot have a quadratic term with variables from different models")
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
            term_str[2*ind-1] = V[ind] < 0 ? " - " : " + " 
            
            x = getName(Variable(q.qvars1[ind].m,I[ind]))
            v = abs(V[ind])
            if I[ind] == J[ind]
                # Squared term
                term_str[2*ind] = (v==1.0 ? "" : str_round(v)*" ") * "$x$sq"
            else
                # Normal term
                y = getName(Variable(q.qvars1[ind].m,J[ind]))
                term_str[2*ind] = (v==1.0 ? "" : str_round(v)*" ") * "$x$times$y"
            end
        end
        # Correction for first term
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
quadToStr(c::GenericQuadExpr) = quad_str(REPLMode,c)


#########################################################################
# CONSTRAINTS
#########################################################################
#------------------------------------------------------------------------
# GenericRangeConstraint
#------------------------------------------------------------------------
Base.print(io::IO, c::GenericRangeConstraint) =
    print(io, con_str(REPLMode,c))
Base.show(io::IO, c::GenericRangeConstraint) =
    show(io, con_str(REPLMode,c))
Base.writemime(io::IO, ::MIME"text/latex", c::GenericRangeConstraint) =
    print(io, con_str(IJuliaMode,c,mathmode=false))

# Generic constraint printer
function con_str(mode, c::GenericRangeConstraint, leq, eq, geq)
    s = sense(c)
    a = aff_str(mode,c.terms,false)
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
Base.print(io::IO, c::QuadConstraint) = 
    print(io, con_str(c))
Base.show(io::IO, c::QuadConstraint) =
    show(io, con_str(c))
Base.writemime(io::IO, ::MIME"text/latex", c::QuadConstraint) =
    print(io, con_str(IJuliaMode,c,mathmode=false))

# Generic constraint printer
function con_str(mode, c::QuadConstraint, leq, eq, geq)
    s = c.sense
    r = s == :<= ? leq : (s == :>= ? geq : eq)
    "$(quad_str(mode,c.terms)) $r 0"
end
# Backwards compatability shim
conToStr(c::QuadConstraint) = con_str(REPLMode,c)
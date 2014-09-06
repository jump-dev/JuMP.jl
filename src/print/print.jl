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
quad_str(mode,q::QuadExpr) = quadToStr(q)

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
function con_str(mode,c::GenericRangeConstraint, leq, eq, geq)
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
function con_str(mode,c::QuadConstraint, leq, eq, geq)
    s = c.sense
    r = s == :<= ? leq : (s == :>= ? geq : eq)
    "$(quad_str(mode,c.terms)) $r 0"
end
# Backwards compatability shim
conToStr(c::QuadConstraint) = con_str(REPLMode,c)
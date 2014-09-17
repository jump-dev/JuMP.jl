#############################################################################
# JuMP
# An algebraic modelling langauge for Julia
# See http://github.com/JuliaOpt/JuMP.jl
#############################################################################
# print/repl.jl
# Pretty printing for REPL
#############################################################################

const repl_leq = @windows? "<=" : "≤?"
const repl_geq = @windows? ">=" : "≥?"
const repl_eq  = @windows? "==" : "=?"
const repl_times = "*"
const repl_sq    = "\u00B2"  # Superscript 2

#########################################################################
# EXPRESSIONS
#########################################################################
#------------------------------------------------------------------------
# AffExpr  (not GenericAffExpr)
#------------------------------------------------------------------------
aff_str(::Type{REPLMode}, a::AffExpr; show_constant=true) =
    aff_str(a, show_constant=show_constant)
#------------------------------------------------------------------------
# GenericQuadExpr
#------------------------------------------------------------------------
quad_str(::Type{REPLMode}, q::GenericQuadExpr) = 
    quad_str(REPLMode, q, repl_times, repl_sq)

#########################################################################
# CONSTRAINTS
#########################################################################
#------------------------------------------------------------------------
# GenericRangeConstraint
#------------------------------------------------------------------------
con_str(::Type{REPLMode}, c::GenericRangeConstraint) =
    con_str(REPLMode, c, repl_leq, repl_eq, repl_geq)

#------------------------------------------------------------------------
# QuadConstraint
#------------------------------------------------------------------------
con_str(::Type{REPLMode}, c::QuadConstraint) =
    con_str(REPLMode, c, repl_leq, repl_eq, repl_geq)
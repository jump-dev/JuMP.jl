#############################################################################
# JuMP
# An algebraic modelling langauge for Julia
# See http://github.com/JuliaOpt/JuMP.jl
#############################################################################
# print/repl.jl
# Pretty printing for REPL
#############################################################################

const repl_leq = @windows? "<=" : "≤"
const repl_geq = @windows? ">=" : "≥"
const repl_eq  = @windows? "==" : "="

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
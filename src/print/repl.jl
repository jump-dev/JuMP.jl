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

#########################################################################
# VARIABLES
#########################################################################
#------------------------------------------------------------------------
## JuMPContainer{Variable}
#------------------------------------------------------------------------
cont_str(::Type{REPLMode}, j::JuMPContainer{Variable}) =
    cont_str(REPLMode, j, repl_leq, repl_eq, repl_geq,
                        repl_ind_open, repl_ind_close, repl_for_all, repl_in,
                        repl_open_set, repl_mid_set, repl_close_set,
                        repl_union, repl_infty, repl_open_rng, repl_close_rng,
                        repl_integer)

#########################################################################
# EXPRESSIONS
#########################################################################
#------------------------------------------------------------------------
## AffExpr  (not GenericAffExpr)
#------------------------------------------------------------------------
aff_str(::Type{REPLMode}, a::AffExpr; show_constant=true) =
    aff_str(a, show_constant=show_constant)
#------------------------------------------------------------------------
## GenericQuadExpr
#------------------------------------------------------------------------
quad_str(::Type{REPLMode}, q::GenericQuadExpr) = 
    quad_str(REPLMode, q, repl_times, repl_sq)

#########################################################################
# CONSTRAINTS
#########################################################################
#------------------------------------------------------------------------
## GenericRangeConstraint
#------------------------------------------------------------------------
con_str(::Type{REPLMode}, c::GenericRangeConstraint) =
    con_str(REPLMode, c, repl_leq, repl_eq, repl_geq)
#------------------------------------------------------------------------
## QuadConstraint
#------------------------------------------------------------------------
con_str(::Type{REPLMode}, c::QuadConstraint) =
    con_str(REPLMode, c, repl_leq, repl_eq, repl_geq)
#------------------------------------------------------------------------
## SOSConstraint
#------------------------------------------------------------------------
con_str(::Type{REPLMode}, c::SOSConstraint) =
    con_str(REPLMode, c, repl_open_set, repl_close_set)
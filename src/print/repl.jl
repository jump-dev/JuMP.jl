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

#------------------------------------------------------------------------
## Model
#------------------------------------------------------------------------
model_str(::Type{REPLMode}, m::Model) =
    model_str(REPLMode, m, repl_leq, repl_geq, repl_in,
                        repl_open_set, repl_mid_set, repl_close_set,
                        repl_union, repl_infty, repl_open_rng, repl_close_rng,
                        repl_integer)
#------------------------------------------------------------------------
## Variable
#------------------------------------------------------------------------
var_str(::Type{REPLMode}, v::Variable) = var_str(REPLMode, v.m, v.col)
var_str(::Type{REPLMode}, m::Model, col::Int) = 
    var_str(REPLMode, m, col, repl_ind_open, repl_ind_close)
#------------------------------------------------------------------------
## JuMPContainer{Variable}
#------------------------------------------------------------------------
cont_str(::Type{REPLMode}, j::JuMPContainer{Variable}; mathmode=false) =
    cont_str(REPLMode, j, repl_leq, repl_eq, repl_geq,
                        repl_ind_open, repl_ind_close, repl_for_all, repl_in,
                        repl_open_set, repl_mid_set, repl_close_set,
                        repl_union, repl_infty, repl_open_rng, repl_close_rng,
                        repl_integer)
#------------------------------------------------------------------------
## GenericQuadExpr
#------------------------------------------------------------------------
quad_str(::Type{REPLMode}, q::GenericQuadExpr) = 
    quad_str(REPLMode, q, repl_times, repl_sq)
#------------------------------------------------------------------------
## GenericRangeConstraint
#------------------------------------------------------------------------
con_str(::Type{REPLMode}, c::GenericRangeConstraint; args...) =
    con_str(REPLMode, c, repl_leq, repl_eq, repl_geq)
#------------------------------------------------------------------------
## QuadConstraint
#------------------------------------------------------------------------
con_str(::Type{REPLMode}, c::QuadConstraint; args...) =
    con_str(REPLMode, c, repl_leq, repl_eq, repl_geq)
#------------------------------------------------------------------------
## SOSConstraint
#------------------------------------------------------------------------
con_str(::Type{REPLMode}, c::SOSConstraint; args...) =
    con_str(REPLMode, c, repl_open_set, repl_close_set)
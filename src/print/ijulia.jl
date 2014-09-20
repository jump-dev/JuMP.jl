#############################################################################
# JuMP
# An algebraic modelling langauge for Julia
# See http://github.com/JuliaOpt/JuMP.jl
#############################################################################
# print/ijulia.jl
# Pretty printing for IJulia (MathJax)
#############################################################################

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
model_str(::Type{IJuliaMode}, m::Model; mathmode=true) =
    math(model_str(IJuliaMode, m, ijulia_leq, ijulia_geq, ijulia_in,
                        ijulia_open_set, ijulia_mid_set, ijulia_close_set, 
                        ijulia_union, ijulia_infty, ijulia_open_rng, ijulia_close_rng,
                        ijulia_integer), mathmode)
#------------------------------------------------------------------------
## JuMPContainer{Variable}
#------------------------------------------------------------------------
cont_str(::Type{IJuliaMode}, j::JuMPContainer{Variable}; mathmode=true) =
    math(cont_str(IJuliaMode, j, ijulia_leq, ijulia_eq, ijulia_geq,
                        ijulia_ind_open, ijulia_ind_close, ijulia_for_all, ijulia_in,
                        ijulia_open_set, ijulia_mid_set, ijulia_close_set, 
                        ijulia_union, ijulia_infty, ijulia_open_rng, ijulia_close_rng,
                        ijulia_integer), mathmode)
#------------------------------------------------------------------------
## AffExpr  (not GenericAffExpr)
#------------------------------------------------------------------------
aff_str(::Type{IJuliaMode}, a::AffExpr; mathmode=true, show_constant=true) =
    math(aff_str(a, show_constant=show_constant), mathmode)
#------------------------------------------------------------------------
## GenericQuadExpr
#------------------------------------------------------------------------
quad_str(::Type{IJuliaMode}, q::GenericQuadExpr; mathmode=true) =
    math(quad_str(IJuliaMode, q, ijulia_times, ijulia_sq), mathmode)
#------------------------------------------------------------------------
## GenericRangeConstraint
#------------------------------------------------------------------------
con_str(::Type{IJuliaMode}, c::GenericRangeConstraint; mathmode=true) =
    math(con_str(IJuliaMode, c, ijulia_leq, ijulia_eq, ijulia_geq), mathmode)
#------------------------------------------------------------------------
## QuadConstraint
#------------------------------------------------------------------------
con_str(::Type{IJuliaMode}, c::QuadConstraint; mathmode=true) =
    math(con_str(IJuliaMode, c, ijulia_leq, ijulia_eq, ijulia_geq), mathmode)
#------------------------------------------------------------------------
## SOSConstraint
#------------------------------------------------------------------------
con_str(::Type{IJuliaMode}, c::SOSConstraint; mathmode=true) =
    math(con_str(IJuliaMode, c, ijulia_open_set, ijulia_close_set), mathmode)
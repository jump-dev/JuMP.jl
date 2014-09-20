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
math(s) = "\$\$ $s \$\$"


#########################################################################
# VARIABLES
#########################################################################
#------------------------------------------------------------------------
## JuMPContainer{Variable}
#------------------------------------------------------------------------
function cont_str(::Type{IJuliaMode}, j::JuMPContainer{Variable}; mathmode=true)
    j_str = cont_str(IJuliaMode, j, ijulia_leq, ijulia_eq, ijulia_geq,
                        ijulia_ind_open, ijulia_ind_close, ijulia_for_all, ijulia_in,
                        ijulia_open_set, ijulia_mid_set, ijulia_close_set, 
                        ijulia_union, ijulia_infty, ijulia_open_rng, ijulia_close_rng,
                        ijulia_integer)
    return mathmode ? j_str : math(j_str)
end

#########################################################################
# EXPRESSIONS
#########################################################################
#------------------------------------------------------------------------
## AffExpr  (not GenericAffExpr)
#------------------------------------------------------------------------
function aff_str(::Type{IJuliaMode}, a::AffExpr; mathmode=true, show_constant=true)
    a_str = aff_str(a, show_constant=show_constant)
    mathmode ? a_str : math(a_str)
end
#------------------------------------------------------------------------
## GenericQuadExpr
#------------------------------------------------------------------------
function quad_str(::Type{IJuliaMode}, q::GenericQuadExpr; mathmode=true)
    q_str = quad_str(IJuliaMode, q, ijulia_times, ijulia_sq)
    mathmode ? q_str : math(q_str)
end

#########################################################################
# CONSTRAINTS
#########################################################################
#------------------------------------------------------------------------
## GenericRangeConstraint
#------------------------------------------------------------------------
function con_str(::Type{IJuliaMode}, c::GenericRangeConstraint; mathmode=true)
    c_str = con_str(IJuliaMode, c, ijulia_leq, ijulia_eq, ijulia_geq)
    mathmode ? c_str : math(c_str)
end
#------------------------------------------------------------------------
## QuadConstraint
#------------------------------------------------------------------------
function con_str(::Type{IJuliaMode}, c::QuadConstraint; mathmode=true)
    c_str = con_str(IJuliaMode, c, ijulia_leq, ijulia_eq, ijulia_geq)
    mathmode ? c_str : math(c_str)
end
#------------------------------------------------------------------------
## SOSConstraint
#------------------------------------------------------------------------
function con_str(::Type{IJuliaMode}, c::SOSConstraint; mathmode=true)
    c_str = con_str(IJuliaMode, c, ijulia_open_set, ijulia_close_set)
    mathmode ? c_str : math(c_str)
end
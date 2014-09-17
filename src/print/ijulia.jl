#############################################################################
# JuMP
# An algebraic modelling langauge for Julia
# See http://github.com/JuliaOpt/JuMP.jl
#############################################################################
# print/ijulia.jl
# Pretty printing for IJulia (MathJax)
#############################################################################

const ijulia_leq    = "\\leq"
const ijulia_geq    = "\\geq"
const ijulia_eq     = "="
const ijulia_times  = "\\times"
const ijulia_sq     = "^2"
math(s) = "\$\$ $s \$\$"

#########################################################################
# EXPRESSIONS
#########################################################################
#------------------------------------------------------------------------
# AffExpr  (not GenericAffExpr)
#------------------------------------------------------------------------
function aff_str(::Type{IJuliaMode}, a::AffExpr; mathmode=true, show_constant=true)
    a_str = aff_str(a, show_constant=show_constant)
    mathmode ? a_str : math(a_str)
end
#------------------------------------------------------------------------
# GenericQuadExpr
#------------------------------------------------------------------------
function quad_str(::Type{IJuliaMode}, q::GenericQuadExpr; mathmode=true)
    q_str = quad_str(IJuliaMode, q, ijulia_times, ijulia_sq)
    mathmode ? q_str : math(q_str)
end

#########################################################################
# CONSTRAINTS
#########################################################################
#------------------------------------------------------------------------
# GenericRangeConstraint
#------------------------------------------------------------------------
function con_str(::Type{IJuliaMode}, c::GenericRangeConstraint; mathmode=true)
    c_str = con_str(IJuliaMode, c, ijulia_leq, ijulia_eq, ijulia_geq)
    mathmode ? c_str : math(c_str)
end

#------------------------------------------------------------------------
# QuadConstraint
#------------------------------------------------------------------------
function con_str(::Type{IJuliaMode}, c::QuadConstraint; mathmode=true)
    c_str = con_str(IJuliaMode, c, ijulia_leq, ijulia_eq, ijulia_geq)
    mathmode ? c_str : math(c_str)
end
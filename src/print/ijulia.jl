#############################################################################
# JuMP
# An algebraic modelling langauge for Julia
# See http://github.com/JuliaOpt/JuMP.jl
#############################################################################
# print/ijulia.jl
# Pretty printing for IJulia (MathJax)
#############################################################################

math(s) = "\$\$ $s \$\$"


#########################################################################
# EXPRESSIONS
#########################################################################
#------------------------------------------------------------------------
# GenericQuadExpr
#------------------------------------------------------------------------
function quad_str(::Type{IJuliaMode}, q::GenericQuadExpr; mathmode=true)
    q_str = quad_str(IJuliaMode, q, "\\times", "^2")
    mathmode ? q_str : math(q_str)
end

#########################################################################
# CONSTRAINTS
#########################################################################
#------------------------------------------------------------------------
# GenericRangeConstraint
#------------------------------------------------------------------------
function con_str(::Type{IJuliaMode}, c::GenericRangeConstraint; mathmode=true)
    c_str = con_str(IJuliaMode, c, "\\leq", "=", "\\geq")
    mathmode ? c_str : math(c_str)
end

#------------------------------------------------------------------------
# QuadConstraint
#------------------------------------------------------------------------
function con_str(::Type{IJuliaMode}, c::QuadConstraint; mathmode=true)
    c_str = con_str(IJuliaMode, c, "\\leq", "=", "\\geq")
    mathmode ? c_str : math(c_str)
end
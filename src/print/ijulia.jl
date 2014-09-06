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
# CONSTRAINTS
#########################################################################
#------------------------------------------------------------------------
# GenericRangeConstraint
#------------------------------------------------------------------------
function con_str(::Type{IJuliaMode}, c::GenericRangeConstraint; mathmode=true)
    c = con_str(IJuliaMode, c, "\\leq", "=", "\\geq")
    mathmode ? c : math(c)
end

#------------------------------------------------------------------------
# QuadConstraint
#------------------------------------------------------------------------
function con_str(::Type{IJuliaMode}, c::QuadConstraint; mathmode=true)
    c = con_str(IJuliaMode, c, "\\leq", "=", "\\geq")
    mathmode ? c : math(c)
end
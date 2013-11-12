#############################################################################
# JuMP
# An algebraic modelling langauge for Julia
# See http://github.com/JuliaOpt/JuMP.jl
#############################################################################
# basic.jl
#
# Solves a simple LP:
# max 5x + 3y
#  st 1x + 5y <= 3
#      0 <= x <= 2
#      0 <= y <= 20
#############################################################################

using JuMP

m = Model()

@defVar(m, 0 <= x <= 2)
@defVar(m, 0 <= y <= 30)

@setObjective(m, Max, 5x + 3y)
@addConstraint(m, 1x + 5y <= 3.0)
    
print(m)
    
status = solve(m)
    
println("Objective value: ", getObjectiveValue(m))
println("x = ", getValue(x))
println("y = ", getValue(y))

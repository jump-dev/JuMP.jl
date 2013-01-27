########################################################################
# Jump 
# A MILP+QP modelling langauge for Julia
# Julia + Mathematical Programming = Jump
# By Iain Dunning and Miles Lubin
#
# knapsack.jl
# A simple knapsack-style problem.
# 
########################################################################


require("../src/Jump.jl")
using Jump

m = Model("max")

x = addVars(m, 0, 1, JUMP_CONTINUOUS, 5, "x")
y = addVars(m, 0, 1, JUMP_CONTINUOUS, (5,5), "y")

profit = [ 5, 3, 2, 7, 4 ]
weight = [ 2, 8, 4, 2, 5 ]
capacity = 10

@setObjective(m, sum{profit[i]*x[i], i=1:5} )

@addConstraint(m, sum{weight[i]*x[i], i=1:5} <= capacity)

status = solveClp(m)

println("Objective is: ",m.objVal)
println("Solution is:")
for i = 1:5
	println("x", i, " = ", getValue(x[i]))
end

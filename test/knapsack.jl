# Julp - LP for Julia
# knapsack.jl
# A simple example of modelling in Julia

require("../src/Julp.jl")
using Julp

m = Model("max")


x = [ Variable(m, 0, 1, 0, "x$i") for i=1:5 ]

profit = [ 5, 3, 2, 7, 4 ]
weight = [ 2, 8, 4, 2, 5 ]
capacity = 10

setObjective(m,  @sumExpr([ profit[i]*x[i] for i = 1:5 ]) )

addConstraint(m, @sumExpr([ weight[i]*x[i] for i = 1:5 ]) <= capacity)

status = solveClp(m)

println("Objective is: ",m.objVal)
println("Solution is:")
for i = 1:5
	println("x", i, " = ", getValue(x[i]))
end

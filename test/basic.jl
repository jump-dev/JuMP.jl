require("../src/MathProg.jl")
using MathProg

m = Model("max")

x = addVar(m, 0,  2, CONTINUOUS)
y = addVar(m, 0, 30, CONTINUOUS, "y")

setObjective(m, 5x + 3y)
addConstraint(m, 1x + 5y <= 3.0)
    
print(m)
    
status = solveClp(m)
    
println("Objective value: ", m.objVal)
println("x = ", getValue(x))
println("y = ", getValue(y))

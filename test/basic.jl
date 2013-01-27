require("../src/Jump.jl")
using Jump

m = Model(JUMP_MAX)

x = addVar(m, 0,  2, JUMP_CONTINUOUS)
y = addVar(m, 0, 30, JUMP_CONTINUOUS, "y")

setObjective(m, 5x + 3y)
addConstraint(m, 1x + 5y <= 3.0)
    
print(m)
    
status = solveClp(m)
    
println("Objective value: ", m.objVal)
println("x = ", getValue(x))
println("y = ", getValue(y))

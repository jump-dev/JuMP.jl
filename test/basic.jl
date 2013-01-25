require("../src/Jump.jl")
using Jump

m = Model(JUMP_MAX)

x = Variable(m, 0, 2, JUMP_CONTINUOUS)
y = Variable(m, 0, 30, JUMP_CONTINUOUS, "y")

setObjective(m, 5*x + 3*y)
addConstraint(m, 1*x + 5*y <= 3.0)
    
print(m)
    
status = solveClp(m)
    
println("Objective value: ", m.objVal)
println("x = ", getValue(x))
println("y = ", getValue(y))

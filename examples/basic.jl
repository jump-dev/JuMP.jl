using MathProg

m = Model(:Max)

@defVar(m, 0 <= x <= 2)
@defVar(m, 0 <= y <= 30)

@setObjective(m, 5x + 3y)
@addConstraint(m, 1x + 5y <= 3.0)
    
print(m)
    
status = solve(m)
    
println("Objective value: ", m.objVal)
println("x = ", getValue(x))
println("y = ", getValue(y))

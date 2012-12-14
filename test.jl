require("julp.jl")

m = Julp.Model("max")

x = Julp.Variable(m, 0., 3.,    0,"x")
y = Julp.Variable(m, 0., 1e10, 0)
z = Julp.Variable(m, 0., 33.,  0)

#Julp.lpSum([1*x,2*y+z])

Julp.setObjective(m, 5*x + 3*y)
Julp.addConstraint(m, 1.0*x + 5*y <= 3.0)
Julp.print(m)
Julp.writeLP(m,"test.lp")

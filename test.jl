require("julp.jl")

m = Julp.Model("max")

x = Julp.Variable(m, 0., 3.,   0,"x")
y = Julp.Variable(m, 0., Julp.Infinity, 0,"y")
#z = Julp.Variable(m, 0., Julp.Infinity,  0,"z")

#Julp.lpSum([1*x,2*y+z])

#qe = (z+2)*(x+3) + y^2
#Julp.print(qe)
#println("")


Julp.setObjective(m, 5*x + 3*y)
Julp.addConstraint(m, 1.0*x + 5*y <= 3.0)
Julp.print(m)
Julp.writeLP(m,"test.lp")
Julp.writeMPS(m,"test.mps")

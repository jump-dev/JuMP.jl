require("julp.jl")

m = Julp.Model("max")

x = Julp.Variable(m,"x",0,3)
y = Julp.Variable(m,0,1e10)

m.objective = 5*x + 3*y
Julp.AddConstraint(m, 1.0*x + 5*y <= 3.0)
Julp.PrintModel(m)
Julp.WriteLP(m,"test.lp")

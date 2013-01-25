Jump
====

Linear Programming, Quadratic Programming and Integer Programming 
modelling with Julia. The goal is to have the speed of AMPL embedded in
a fully functional language. Compare with PuLP/Pyomo/YALMIP/...

    Julia + Mathematical Programming = JuMP

# Installation

You can install Jump through the Julia package manager:
    Pkg.add("Jump")
    
Jump is currently dependent on the CLP package which provides a powerful
open-source LP solver. As the infrastructure and interfaces for solvers
are developed in Julia we will extend the functionality

# Simple Example

    using Jump

    m = Model(JUMP_MAX)

    x = addVar(m, 0,  2, JUMP_CONTINUOUS)
    y = addVar(m, 0, 30, JUMP_CONTINUOUS, "y")

    setObjective(m, 5*x + 3*y)
    addConstraint(m, 1*x + 5*y <= 3.0)
    
    print(m)
    
    status = solveClp(m)
    
    println("Objective value: ", m.objVal)
    println("x = ", getValue(x))
    println("y = ", getValue(y))
    
# More advanced

You can see more examples in the ```test/``` folder.



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
    
# Getting Speed

Using Julia's powerful metaprogramming features, we can turn easy-to-read
statements into the a sparse internal representation very quickly. To 
invoke this, use the @sumExpr command. Here are some examples:

    addConstraint(m, @sumExpr(1x[i] + -1s[i]) <= 0)
    
    setObjective(m, @sumExpr([1.0*x[i] for i = 1:numLocation]) )
    
There are some restrictions on what can go inside the @sumExpr(...)
 * You can build complicated expressions, but they must always simplify
   down to Coefficient*Variable
 * Variables must always have a coefficient, even if it is just 1
 * No standalone constants
 
# Full function listing

 * ```Model(sense)```
   * Construct a Model with the objective sense provided. There are two
     constants exported by Jump, JUMP_MAXIMIZE and JUMP_MINIMIZE, that
     you should use.
 
 * print(Model)
   * Displays the Model to the screen in a user-readable way





Julp
====

MILP modelling with Julia


    Julia
    +  LP
    -----
    =Julp.

# Usage Example

Example code (so far)

    m = Model("max")

    x = Variable(m,"x",0,1e10)
    y = Variable(m,"y",0,1e10)

    m.objective = 5*x + 3*y
    AddConstraint(m, 1.0*x + 5*y <= 3.0)
    PrintModel(m)



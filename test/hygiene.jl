# hygiene.jl
# Make sure that our macros have good hygiene
using Base.Test

require("MathProg")

mymod = MathProg.Model(:Min)
MathProg.@defVar(mymod, x >= 0)
r = 3:5
MathProg.@defVar(mymod, y[i=r] <= i)
MathProg.@addConstraint(mymod, x + sum{ j*y[j], j=r } <= 1)
MathProg.@addConstraint(mymod, sum{ y[j], j=r ; j == 4} <= 1)
MathProg.@addConstraint(mymod, -1 <= x + y[3] <= 1)
MathProg.@setObjective(mymod, y[4])

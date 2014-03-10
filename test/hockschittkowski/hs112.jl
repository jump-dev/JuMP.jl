#############################################################################
# JuMP
# An algebraic modelling langauge for Julia
# See http://github.com/JuliaOpt/JuMP.jl
#############################################################################
# Hock-Schittkowski Nonlinear Test Suite  -  HS112
# This file is JuMP implementation of the model described in 
#  W. Hock, K. Schittkowski, Test Examples for Nonlinear Programming
#  Codes, Lecture Notes in Economics and Mathematical Systems, 
#  Springer, No, 187, 1981 
# More information, including original model description, at
# http://www.ai7.uni-bayreuth.de/downloads.htm
#
# This problem has an objective with the logarithm of a fraction where both
# the nominator and denominator have variables in them. Constraints linear.
#############################################################################

using JuMP
using Base.Test

let

c = [-6.089, -17.164, -34.054, -5.914, -24.721, -14.986, -24.100, -10.708, -26.662, -22.179]

m = Model()
@defVar(m, x[1:10] >= 1e-6)
for i = 1:10
    setValue(x[i], 0.1)
end

@setNLObjective(m, Min, sum{x[j]*(c[j] + log(x[j]/sum{x[k],k=1:10})), j=1:10})

@addNLConstraint(m, x[1] + 2*x[2] + 2*x[3] + x[6] + x[10] == 2)
@addNLConstraint(m, x[4] + 2*x[5] + x[6] + x[7] == 1)
@addNLConstraint(m, x[3] + x[7] + x[8] + 2*x[9] + x[10] == 1)

solve(m)

@test_approx_eq_eps getObjectiveValue(m) -47.76109026 1e-5

end

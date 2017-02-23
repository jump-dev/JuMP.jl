#  Copyright 2017, Iain Dunning, Joey Huchette, Miles Lubin, and contributors
#  This Source Code Form is subject to the terms of the Mozilla Public
#  License, v. 2.0. If a copy of the MPL was not distributed with this
#  file, You can obtain one at http://mozilla.org/MPL/2.0/.
#############################################################################
# JuMP
# An algebraic modelling langauge for Julia
# See http://github.com/JuliaOpt/JuMP.jl
#############################################################################
# Hock-Schittkowski Nonlinear Test Suite  -  HS114
# This file is JuMP implementation of the model described in
#  W. Hock, K. Schittkowski, Test Examples for Nonlinear Programming
#  Codes, Lecture Notes in Economics and Mathematical Systems,
#  Springer, No, 187, 1981
# More information, including original model description, at
# http://www.ai7.uni-bayreuth.de/downloads.htm
#
# This problem has a quadratic objective with quadratic constraints as well
# as a couple constraints that have fractions of variables.
#############################################################################

@testset "HS114" begin

n = 10
a = 0.99
b = 0.9

lower = [0.0001, 0.0001, 0.0001, 0.0001, 0.0001, 85, 90, 3, 1.2, 145]
upper = [2000, 16000, 120, 5000, 2000, 93, 95, 12, 4, 162]
start = [1745, 12000, 110, 3048, 1974, 89.2, 92.8, 8, 3.6, 145]

m = Model(solver=nlp_solver)
@variable(m, lower[i] <= x[i=1:n] <= upper[i], start = start[i])

@NLobjective(m, Min, 5.04*x[1] + .035*x[2] + 10*x[3] + 3.36*x[5] - .063*x[4]*x[7])

@NLconstraint(m, 35.82 - .222*x[10] - b*x[9] >= 0)
@NLconstraint(m, -133 + 3*x[7] - a*x[10] >= 0)
@NLconstraint(m, -(35.82 - .222*x[10] - b*x[9]) + x[9]*(1/b - b) >= 0)
@NLconstraint(m, -(-133 + 3*x[7] - a*x[10]) + (1/a - a)*x[10] >= 0)
@NLconstraint(m, 1.12*x[1] + .13167*x[1]*x[8] - .00667*x[1]*x[8]^2 - a*x[4] >= 0)
@NLconstraint(m, 57.425 + 1.098*x[8] - .038*x[8]^2 + .325*x[6] - a*x[7] >= 0)
@NLconstraint(m, -(1.12*x[1] + .13167*x[1]*x[8] - .00667*x[1]*x[8]^2 - a*x[4]) + (1/a - a)*x[4] >= 0)
@NLconstraint(m, -(57.425 + 1.098*x[8] - .038*x[8]^2 + .325*x[6] - a*x[7]) + (1/a - a)*x[7] >= 0)
@NLconstraint(m, 1.22*x[4] - x[1] - x[5] == 0)
@NLconstraint(m, 98000*x[3]/(x[4]*x[9] + 1000*x[3]) - x[6] == 0)
@NLconstraint(m, (x[2] + x[5])/x[1] - x[8] == 0)

solve(m)

@test isapprox(getobjectivevalue(m), -1768.80696, atol=1e-3)

end

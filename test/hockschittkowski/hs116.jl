#  Copyright 2017, Iain Dunning, Joey Huchette, Miles Lubin, and contributors
#  This Source Code Form is subject to the terms of the Mozilla Public
#  License, v. 2.0. If a copy of the MPL was not distributed with this
#  file, You can obtain one at http://mozilla.org/MPL/2.0/.
#############################################################################
# JuMP
# An algebraic modelling langauge for Julia
# See http://github.com/JuliaOpt/JuMP.jl
#############################################################################
# Hock-Schittkowski Nonlinear Test Suite  -  HS116
# This file is JuMP implementation of the model described in
#  W. Hock, K. Schittkowski, Test Examples for Nonlinear Programming
#  Codes, Lecture Notes in Economics and Mathematical Systems,
#  Springer, No, 187, 1981
# More information, including original model description, at
# http://www.ai7.uni-bayreuth.de/downloads.htm
#
# This problem has a linear objective with quadratic constraints.
#############################################################################

@testset "HS116" begin

N = 13
a = 0.002
b = 1.262626
c = 1.231059
d = 0.03475
e = 0.975
f = 0.00975

lower = [0.1, 0.1, 0.1, 0.0001, 0.1, 0.1, 0.1, 0.1, 500, 0.1, 1.0, 0.0001, 0.0001, 0.0, 0.0, 0.0]
upper = [1.0, 1.0, 1.0, 0.1, 0.9, 0.9, 1000, 1000, 1000, 500, 150, 150, 150, Inf, Inf, Inf]
start = [0.5  2 0.8  3 0.9  4 0.1  5 0.14  6 0.5  7 489  8 80  9 650 0.5  2 0.8  3 0.9  4 0.1  5 0.14  6 0.5  7 489  8 80  9 650]

m = Model(solver=nlp_solver)
@variable(m, lower[i] <= x[i=1:N] <= upper[i], start = start[i])
@NLobjective(m, Min, x[11] + x[12] + x[13])

@NLconstraint(m, x[3] - x[2] >= 0)
@NLconstraint(m, x[2] - x[1] >= 0)
@NLconstraint(m, 1 - a * x[7] + a * x[8] >= 0)
@NLconstraint(m, x[11] + x[12] + x[13] >= 50)
@NLconstraint(m, x[13] - b * x[10] + c * x[3] * x[10] >= 0)
@NLconstraint(m, x[5] - d * x[2] - e * x[2] * x[5] + f * x[2]^2 >= 0)
@NLconstraint(m, x[6] - d * x[3] - e * x[3] * x[6] + f * x[3]^2 >= 0)
@NLconstraint(m, x[4] - d * x[1] - e * x[1] * x[4] + f * x[1]^2 >= 0)
@NLconstraint(m, x[12] - b * x[9] + c * x[2] * x[9] >= 0)
@NLconstraint(m, x[11] - b * x[8] + c * x[1] * x[8] >= 0)
@NLconstraint(m, x[5] * x[7] - x[1] * x[8] - x[4] * x[7] + x[4] * x[8] >= 0)
@NLconstraint(m, 1 - a * (x[2] * x[9] + x[5] * x[8] - x[1] * x[8] - x[6] * x[9]) -
          x[5] - x[6] >= 0)
@NLconstraint(m, x[2] * x[9] - x[3] * x[10] - x[6] * x[9] - 500 * x[2] +
          500 * x[6] + x[2] * x[10] >= 0)
@NLconstraint(m, x[2] - 0.9 - a * (x[2] * x[10] - x[3] * x[10]) >= 0)
@NLconstraint(m, x[11] + x[12] + x[13] <= 250)


solve(m)

@test isapprox(getobjectivevalue(m), 97.588409, atol=1e-3)

end

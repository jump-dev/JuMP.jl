#  Copyright 2017, Iain Dunning, Joey Huchette, Miles Lubin, and contributors
#  This Source Code Form is subject to the terms of the Mozilla Public
#  License, v. 2.0. If a copy of the MPL was not distributed with this
#  file, You can obtain one at http://mozilla.org/MPL/2.0/.
#############################################################################
# JuMP
# An algebraic modelling langauge for Julia
# See http://github.com/JuliaOpt/JuMP.jl
#############################################################################
# Hock-Schittkowski Nonlinear Test Suite  -  HS109
# This file is JuMP implementation of the model described in
#  W. Hock, K. Schittkowski, Test Examples for Nonlinear Programming
#  Codes, Lecture Notes in Economics and Mathematical Systems,
#  Springer, No, 187, 1981
# More information, including original model description, at
# http://www.ai7.uni-bayreuth.de/downloads.htm
#
# Nonconvex objective with cubics, and constraints contain trignometry and
# products of variables.
#############################################################################

@testset "HS109" begin

a  = 50.176
b1 = 0.25
b = sin(b1)
c = cos(b1)

L = [0.0, 0.0, -0.55, -0.55, 196, 196, 196, -400, -400]
U = [Inf, Inf,  0.55,  0.55, 252, 252, 252,  800,  800]

m = Model(solver=nlp_solver)
@variable(m, L[i] <= x[i=1:9] <= U[i], start = 0.0)

@NLobjective(m, Min, 3 * x[1] + 1e-6 * x[1]^3 + 2 * x[2] + .522074e-6 * x[2]^3)

@NLconstraint(m, x[4] - x[3] + 0.55 >= 0)
@NLconstraint(m, x[3] - x[4] + 0.55 >= 0)
@NLconstraint(m, 2250000 - x[1]^2 - x[8]^2 >= 0)
@NLconstraint(m, 2250000 - x[2]^2 - x[9]^2 >= 0)
@NLconstraint(m,
    x[5] * x[6] * sin(-x[3] - .25) + x[5] * x[7] * sin(-x[4] - .25) +
    2 * b * x[5]^2 - a * x[1] + 400 * a == 0)
@NLconstraint(m,
    x[5] * x[6] * sin(x[3] - .25) + x[6] * x[7] * sin(x[3] - x[4] - .25) +
    2 * b * x[6]^2 - a * x[2] + 400 * a == 0)
@NLconstraint(m,
    x[5] * x[7] * sin(x[4] - .25) + x[6] * x[7] * sin(x[4] - x[3] - .25) +
    2 * b * x[7]^2 + 881.779 * a == 0)
@NLconstraint(m,
    a * x[8] + x[5] * x[6] * cos(-x[3] - .25) +
    x[5] * x[7] * cos(-x[4] - .25) - 200 * a - 2 * c * x[5]^2 +
    .7533e-3 * a * x[5]^2 == 0)
@NLconstraint(m,
    a * x[9] + x[5] * x[6] * cos(x[3] - .25) +
    x[6] * x[7] * cos(x[3] - x[4] - .25) - 2 * c * x[6]^2 +
    .7533e-3 * a * x[6]^2 - 200 * a == 0)
@NLconstraint(m,
    x[5] * x[7] * cos(x[4] - .25) + x[6] * x[7] * cos(x[4] - x[3] - .25) -
    2 * c * x[7]^2 + 22.938 * a + .7533e-3 * a * x[7] ^2 == 0)


solve(m)

@test isapprox(getobjectivevalue(m),5326.851310161077, atol=1e-5)

end

#  Copyright 2017, Iain Dunning, Joey Huchette, Miles Lubin, and contributors
#  This Source Code Form is subject to the terms of the Mozilla Public
#  License, v. 2.0. If a copy of the MPL was not distributed with this
#  file, You can obtain one at http://mozilla.org/MPL/2.0/.
#############################################################################
# JuMP
# An algebraic modelling langauge for Julia
# See http://github.com/JuliaOpt/JuMP.jl
#############################################################################
# Hock-Schittkowski Nonlinear Test Suite  -  HS111
# This file is JuMP implementation of the model described in
#  W. Hock, K. Schittkowski, Test Examples for Nonlinear Programming
#  Codes, Lecture Notes in Economics and Mathematical Systems,
#  Springer, No, 187, 1981
# More information, including original model description, at
# http://www.ai7.uni-bayreuth.de/downloads.htm
#
# This problem has an objective and constraints with both logarithms and
# exponents. Related to HS112.
#############################################################################

@testset "HS111" begin

c = [-6.089, -17.164, -34.054, -5.914, -24.721, -14.986, -24.100, -10.708, -26.662, -22.179]

m = Model(solver=nlp_solver)
@variable(m, -100 <= x[1:10] <= 100, start = -2.3)

@NLobjective(m, Min,
    sum(exp(x[j]) * (c[j] + x[j] - log( sum(exp(x[k]) for k=1:10) ) ) for j=1:10))

@NLconstraint(m, exp(x[1]) + 2*exp(x[2]) + 2*exp(x[3]) +   exp(x[6]) + exp(x[10]) == 2)
@NLconstraint(m, exp(x[4]) + 2*exp(x[5]) +   exp(x[6]) +   exp(x[7])              == 1)
@NLconstraint(m, exp(x[3]) +   exp(x[7]) +   exp(x[8]) + 2*exp(x[9]) + exp(x[10]) == 1)

solve(m)

@test isapprox(getobjectivevalue(m), -47.76109026, atol=1e-5)

end

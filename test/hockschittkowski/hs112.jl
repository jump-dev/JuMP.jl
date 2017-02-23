#  Copyright 2017, Iain Dunning, Joey Huchette, Miles Lubin, and contributors
#  This Source Code Form is subject to the terms of the Mozilla Public
#  License, v. 2.0. If a copy of the MPL was not distributed with this
#  file, You can obtain one at http://mozilla.org/MPL/2.0/.
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

@testset "HS112" begin

c = [-6.089, -17.164, -34.054, -5.914, -24.721, -14.986, -24.100, -10.708, -26.662, -22.179]

m = Model(solver=nlp_solver)
@variable(m, x[1:10] >= 1e-6, start = 0.1)

@NLobjective(m, Min, sum(x[j]*(c[j] + log(x[j]/sum(x[k] for k=1:10))) for j=1:10))

@NLconstraint(m, x[1] + 2*x[2] + 2*x[3] + x[6] + x[10] == 2)
@NLconstraint(m, x[4] + 2*x[5] + x[6] + x[7] == 1)
@NLconstraint(m, x[3] + x[7] + x[8] + 2*x[9] + x[10] == 1)

solve(m)

@test isapprox(getobjectivevalue(m), -47.76109026, atol=1e-5)

end

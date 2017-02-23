#  Copyright 2017, Iain Dunning, Joey Huchette, Miles Lubin, and contributors
#  This Source Code Form is subject to the terms of the Mozilla Public
#  License, v. 2.0. If a copy of the MPL was not distributed with this
#  file, You can obtain one at http://mozilla.org/MPL/2.0/.
#############################################################################
# JuMP
# An algebraic modelling langauge for Julia
# See http://github.com/JuliaOpt/JuMP.jl
#############################################################################
# Hock-Schittkowski Nonlinear Test Suite  -  HS110
# This file is JuMP implementation of the model described in
#  W. Hock, K. Schittkowski, Test Examples for Nonlinear Programming
#  Codes, Lecture Notes in Economics and Mathematical Systems,
#  Springer, No, 187, 1981
# More information, including original model description, at
# http://www.ai7.uni-bayreuth.de/downloads.htm
#
# This problem has an objective with squared logarithms and a product
# of variables squared.
#############################################################################

@testset "HS110" begin

m = Model(solver=nlp_solver)
@variable(m, -2.001 <= x[1:10] <= 9.999, start = 9)

@NLobjective(m, Min,
    sum( log(x[j] - 2)^2 + log(10 - x[j])^2 for j=1:10) -
    prod( x[i] for i=1:10) ^ 0.2
)

solve(m)

@test isapprox(getobjectivevalue(m), -45.77846971, atol=1e-5)

end

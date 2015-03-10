#  Copyright 2015, Iain Dunning, Joey Huchette, Miles Lubin, and contributors
#  This Source Code Form is subject to the terms of the Mozilla Public
#  License, v. 2.0. If a copy of the MPL was not distributed with this
#  file, You can obtain one at http://mozilla.org/MPL/2.0/.
#############################################################################
# JuMP
# An algebraic modelling langauge for Julia
# See http://github.com/JuliaOpt/JuMP.jl
#############################################################################
# Hock-Schittkowski Nonlinear Test Suite  -  HS118
# This file is JuMP implementation of the model described in
#  W. Hock, K. Schittkowski, Test Examples for Nonlinear Programming
#  Codes, Lecture Notes in Economics and Mathematical Systems,
#  Springer, No, 187, 1981
# More information, including original model description, at
# http://www.ai7.uni-bayreuth.de/downloads.htm
#
# This problem has a quadratic objective with linear constraints.
#############################################################################

using JuMP
using Base.Test

let

m = Model()

L = zeros(15)
L[1] =  8.0
L[2] = 43.0
L[3] =  3.0

U = zeros(15)
U[1] = 21.0
U[2] = 57.0
U[3] = 16.0
for k in 1:4
    U[3*k+1] =  90.0
    U[3*k+2] = 120.0
    U[3*k+3] =  60.0
end

@defVar(m, L[i] <= x[i=1:15] <= U[i])

@setNLObjective(m, Min,
    sum{2.3     * x[3*k+1]   +
        0.0001  * x[3*k+1]^2 +
        1.7     * x[3*k+2]   +
        0.0001  * x[3*k+2]^2 +
        2.2     * x[3*k+3] +
        0.00015 * x[3*k+3]^2, k=0:4})

# setObjective(m, :Min,
#   sum([(2.3     * x[3*k+1]   +
#         0.0001  * x[3*k+1]^2 +
#         1.7     * x[3*k+2]   +
#         0.0001  * x[3*k+2]^2 +
#         2.2     * x[3*k+3]^2 +
#         0.00015 * x[3*k+3]^2) for k in 0:4 ]))


# constr1
for j in 1:4
    @addConstraint(m, x[3*j+1] - x[3*j-2] + 7 <= 13)
    @addConstraint(m, x[3*j+1] - x[3*j-2] + 7 >=  0)
end

# constr2
for j in 1:4
    @addConstraint(m, x[3*j+2] - x[3*j-1] + 7 <= 14)
    @addConstraint(m, x[3*j+2] - x[3*j-1] + 7 >=  0)
end

# constr3
for j in 1:4
    @addConstraint(m, x[3*j+3] - x[3*j  ] + 7 <= 13)
    @addConstraint(m, x[3*j+3] - x[3*j  ] + 7 >=  0)
end

@addConstraint(m, x[1] + x[2] + x[3]    >= 60)
@addConstraint(m, x[4] + x[5] + x[6]    >= 50)
@addConstraint(m, x[7] + x[8] + x[9]    >= 70)
@addConstraint(m, x[10] + x[11] + x[12] >= 85)
@addConstraint(m, x[13] + x[14] + x[15] >= 100)

# Initial solution (could also use 'start' keyword in @defVar)
setValue(x[1], 20.0)
setValue(x[2], 55.0)
setValue(x[3], 15.0)
setValue(x[4], 20.0)
setValue(x[5], 60.0)
setValue(x[6], 20.0)
setValue(x[7], 20.0)
setValue(x[8], 60.0)
setValue(x[9], 20.0)
setValue(x[10], 20.0)
setValue(x[11], 60.0)
setValue(x[12], 20.0)
setValue(x[13], 20.0)
setValue(x[14], 60.0)
setValue(x[15], 20.0)

solve(m)

#println(getValue(x))

@test_approx_eq_eps getValue(x[1]) 8.0  1e-5
@test_approx_eq_eps getValue(x[2]) 49.0 1e-5
@test_approx_eq_eps getValue(x[3]) 3.0  1e-5
@test_approx_eq_eps getValue(x[4]) 1.0  1e-5
@test_approx_eq_eps getObjectiveValue(m) 664.82045 1e-5

end

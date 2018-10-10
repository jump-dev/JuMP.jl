#  Copyright 2017, Iain Dunning, Joey Huchette, Miles Lubin, and contributors
#  This Source Code Form is subject to the terms of the Mozilla Public
#  License, v. 2.0. If a copy of the MPL was not distributed with this
#  file, You can obtain one at http://mozilla.org/MPL/2.0/.
#############################################################################
# JuMP
# An algebraic modeling langauge for Julia
# See http://github.com/JuliaOpt/JuMP.jl
#############################################################################
# basic.jl
#
# Solves a simple LP:
# max 5x + 3y
#  st 1x + 5y <= 3
#      0 <= x <= 2
#      0 <= y <= 30
#############################################################################

using JuMP, Clp

m = Model(with_optimizer(Clp.Optimizer))

@variable(m, 0 <= x <= 2)
@variable(m, 0 <= y <= 30)

@objective(m, Max, 5x + 3y)
@constraint(m, 1x + 5y <= 3.0)

print(m)

JuMP.optimize!(m)
status = JuMP.termination_status(m)

println("Objective value: ", JuMP.objective_value(m))
println("x = ", JuMP.result_value(x))
println("y = ", JuMP.result_value(y))

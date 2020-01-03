#  Copyright 2017, Iain Dunning, Joey Huchette, Miles Lubin, and contributors
#  This Source Code Form is subject to the terms of the Mozilla Public
#  License, v. 2.0. If a copy of the MPL was not distributed with this
#  file, You can obtain one at http://mozilla.org/MPL/2.0/.

# hygiene.jl
# Make sure that our macros have good hygiene

module M
using LinearAlgebra
using Test
import JuMP

model = JuMP.Model()
sense = JuMP.MathOptInterface.MIN_SENSE
JuMP.@variable(model, x >= 0)
r = 3:5
JuMP.@variable(model, y[i=r] <= i)
JuMP.@variable(model, z[i=1:2, j=1:2], Symmetric)
@test z isa LinearAlgebra.Symmetric

JuMP.@constraint(model, x + sum( j*y[j] for j=r ) <= 1)
JuMP.@constraint(model, sum( y[j] for j=r if j == 4) <= 1)
JuMP.@constraint(model, -1 <= x + y[3] <= 1)
JuMP.@constraints(model, begin
    x + sum( j*y[j] for j=r ) <= 1
    sum( y[j] for j=r if j == 4) <= 1
end)

JuMP.@SDconstraint(model, [x x; -x x] ⪰ 0)

JuMP.@objective(model, sense, y[4])
JuMP.@objective(model, Min, x + sum( j*y[j] for j=r ))
JuMP.@objective(model, Max, sum( y[j] for j=r if j == 4))

JuMP.@NLconstraint(model, y[3] == 1)
JuMP.@NLconstraint(model, x + sum( j*y[j] for j=r ) <= 1)
JuMP.@NLconstraint(model, sum( y[j] for j=r if j == 4) <= 1)
JuMP.@NLconstraints(model, begin
    x + sum( j*y[j] for j=r ) <= 1
    sum( y[j] for j=r if j == 4) <= 1
end)

JuMP.@NLobjective(model, sense, y[4])
JuMP.@NLobjective(model, Min, x + sum( j*y[j] for j=r ))
JuMP.@NLobjective(model, Max, sum( y[j] for j=r if j == 4))

# TODO: Add tests for the content of the model.

model = JuMP.Model()
i = 10
j = 10
JuMP.@expression(model, ex[j = 2:3], sum(i for i in 1:j))
@test ex[2] == 3
@test ex[3] == 6
@test i == 10
@test j == 10

end

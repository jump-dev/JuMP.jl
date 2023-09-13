#  Copyright 2017, Iain Dunning, Joey Huchette, Miles Lubin, and contributors
#  This Source Code Form is subject to the terms of the Mozilla Public
#  License, v. 2.0. If a copy of the MPL was not distributed with this
#  file, You can obtain one at https://mozilla.org/MPL/2.0/.

module TestHygiene

import LinearAlgebra
import JuMP
import Test

# Check that the non-prefixed `Model` is not available in the current scope
Test.@test_throws UndefVarError Model

model = JuMP.Model()
sense = JuMP.MOI.MIN_SENSE
JuMP.@variable(model, x >= 0)
r = 3:5
JuMP.@variable(model, y[i = r] <= i)
JuMP.@variable(model, z[i = 1:2, j = 1:2], Symmetric)
Test.@test z isa LinearAlgebra.Symmetric

JuMP.@constraint(model, x + sum(j * y[j] for j in r) <= 1)
JuMP.@constraint(model, sum(y[j] for j in r if j == 4) <= 1)
JuMP.@constraint(model, -1 <= x + y[3] <= 1)
JuMP.@constraints(model, begin
    x + sum(j * y[j] for j in r) <= 1
    sum(y[j] for j in r if j == 4) <= 1
end)

JuMP.@constraint(model, [x x; -x x] >= 0, JuMP.PSDCone())

JuMP.@objective(model, sense, y[4])
JuMP.@objective(model, Min, x + sum(j * y[j] for j in r))
JuMP.@objective(model, Max, sum(y[j] for j in r if j == 4))

JuMP.@NLconstraint(model, y[3] == 1)
JuMP.@NLconstraint(model, x + sum(j * y[j] for j in r) <= 1)
JuMP.@NLconstraint(model, sum(y[j] for j in r if j == 4) <= 1)
JuMP.@NLconstraints(model, begin
    x + sum(j * y[j] for j in r) <= 1
    sum(y[j] for j in r if j == 4) <= 1
end)

JuMP.@NLobjective(model, sense, y[4])
JuMP.@NLobjective(model, Min, x + sum(j * y[j] for j in r))
JuMP.@NLobjective(model, Max, sum(y[j] for j in r if j == 4))

# TODO: Add tests for the content of the model.

model = JuMP.Model()
i = 10
j = 10
JuMP.@expression(model, ex[j = 2:3], sum(i for i in 1:j))
Test.@test ex[2] == 3
Test.@test ex[3] == 6
Test.@test i == 10
Test.@test j == 10

# Test that `model` is inferred correctly inside macros, despite being a
# non-const global.
m = JuMP.Model()
JuMP.@variable(m, x[1:0])
Test.@test x == JuMP.VariableRef[]

end

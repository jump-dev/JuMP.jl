#  Copyright 2017, Iain Dunning, Joey Huchette, Miles Lubin, and contributors
#  This Source Code Form is subject to the terms of the Mozilla Public
#  License, v. 2.0. If a copy of the MPL was not distributed with this
#  file, You can obtain one at http://mozilla.org/MPL/2.0/.

# hygiene.jl
# Make sure that our macros have good hygiene

module M
import JuMP

model = JuMP.Model()
sense = :Min
JuMP.@variable(model, x >= 0)
r = 3:5
JuMP.@variable(model, y[i=r] <= i)
JuMP.@constraint(model, x + sum( j*y[j] for j=r ) <= 1)
JuMP.@constraint(model, sum( y[j] for j=r if j == 4) <= 1)
JuMP.@constraint(model, -1 <= x + y[3] <= 1)
JuMP.@objective(model, sense, y[4])
JuMP.@NLconstraint(model, y[3] == 1)
        
# TODO: Add tests for the content of the model.

end

#  Copyright 2016, Iain Dunning, Joey Huchette, Miles Lubin, and contributors
#  This Source Code Form is subject to the terms of the Mozilla Public
#  License, v. 2.0. If a copy of the MPL was not distributed with this
#  file, You can obtain one at http://mozilla.org/MPL/2.0/.

# hygiene.jl
# Make sure that our macros have good hygiene
using Base.Test

module M
import JuMP

mymod = JuMP.Model()
mysense = :Min
JuMP.@defVar(mymod, x >= 0)
r = 3:5
JuMP.@defVar(mymod, y[i=r] <= i)
JuMP.@addConstraint(mymod, x + sum{ j*y[j], j=r } <= 1)
JuMP.@addConstraint(mymod, sum{ y[j], j=r ; j == 4} <= 1)
JuMP.@addConstraint(mymod, -1 <= x + y[3] <= 1)
JuMP.@setObjective(mymod, mysense, y[4])
JuMP.@addNLConstraint(mymod, y[3] == 1)

end

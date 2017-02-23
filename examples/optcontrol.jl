#  Copyright 2017, Iain Dunning, Joey Huchette, Miles Lubin, and contributors
#  This Source Code Form is subject to the terms of the Mozilla Public
#  License, v. 2.0. If a copy of the MPL was not distributed with this
#  file, You can obtain one at http://mozilla.org/MPL/2.0/.
using JuMP, Ipopt

solver = IpoptSolver()

# clnlbeam
# Based on AMPL model
# Copyright (C) 2001 Princeton University
# All Rights Reserved
#
# Permission to use, copy, modify, and distribute this software and
# its documentation for any purpose and without fee is hereby
# granted, provided that the above copyright notice appear in all
# copies and that the copyright notice and this
# permission notice appear in all supporting documentation.

#   Source:
#   H. Maurer and H.D. Mittelman,
#   "The non-linear beam via optimal control with bound state variables",
#   Optimal Control Applications and Methods 12, pp. 19-31, 1991.
let
    N     = 1000
    ni    = N
    h     = 1/ni
    alpha = 350

    m = Model(solver=solver)

    @variable(m, -1 <= t[1:(ni+1)] <= 1)
    @variable(m, -0.05 <= x[1:(ni+1)] <= 0.05)
    @variable(m, u[1:(ni+1)])

    @NLobjective(m, Min, sum( 0.5*h*(u[i+1]^2 + u[i]^2) + 0.5*alpha*h*(cos(t[i+1]) + cos(t[i])) for i = 1:ni))

    # cons1
    for i in 1:ni
        @NLconstraint(m, x[i+1] - x[i] - (0.5h)*(sin(t[i+1])+sin(t[i])) == 0)
    end
    # cons2
    for i in 1:ni
        @constraint(m, t[i+1] - t[i] - (0.5h)*u[i+1] - (0.5h)*u[i] == 0)
    end

    solve(m)

end

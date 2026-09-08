# Copyright 2017, Iain Dunning, Joey Huchette, Miles Lubin, and contributors    #src
# This Source Code Form is subject to the terms of the Mozilla Public License   #src
# v.2.0. If a copy of the MPL was not distributed with this file, You can       #src
# obtain one at https://mozilla.org/MPL/2.0/.                                   #src

# # The clnlbeam problem

# Based on an AMPL model by Hande Y. Benson
#
# Copyright (C) 2001 Princeton University
# All Rights Reserved

# Permission to use, copy, modify, and distribute this software and its
# documentation for any purpose and without fee is hereby granted, provided that
# the above copyright notice appear in all copies and that the copyright notice
# and this permission notice appear in all supporting documentation.

# Source: H. Maurer and H.D. Mittelman, "The non-linear beam via optimal control
# with bound state variables", Optimal Control Applications and Methods 12, pp.
# 19-31, 1991.

using JuMP
import Ipopt
import Test  #src

function example_clnlbeam()
    N = 1000
    h = 1 / N
    alpha = 350
    model = Model(Ipopt.Optimizer)
    @variables(model, begin
        -1 <= t[1:(N+1)] <= 1
        -0.05 <= x[1:(N+1)] <= 0.05
        u[1:(N+1)]
    end)
    @NLobjective(
        model,
        Min,
        sum(
            0.5 * h * (u[i+1]^2 + u[i]^2) +
            0.5 * alpha * h * (cos(t[i+1]) + cos(t[i])) for i in 1:N
        ),
    )
    @NLconstraint(
        model,
        [i = 1:N],
        x[i+1] - x[i] - 0.5 * h * (sin(t[i+1]) + sin(t[i])) == 0,
    )
    @constraint(
        model,
        [i = 1:N],
        t[i+1] - t[i] - 0.5 * h * u[i+1] - 0.5 * h * u[i] == 0,
    )
    optimize!(model)
    println("""
    termination_status = $(termination_status(model))
    primal_status      = $(primal_status(model))
    objective_value    = $(objective_value(model))
    """)
    Test.@test termination_status(model) == MOI.LOCALLY_SOLVED  #src
    Test.@test primal_status(model) == MOI.FEASIBLE_POINT  #src
    Test.@test objective_value(model) ≈ 350.0  #src
    return
end

example_clnlbeam()

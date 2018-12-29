#  Copyright 2017, Iain Dunning, Joey Huchette, Miles Lubin, and contributors
#  This Source Code Form is subject to the terms of the Mozilla Public
#  License, v. 2.0. If a copy of the MPL was not distributed with this
#  file, You can obtain one at http://mozilla.org/MPL/2.0/.

using JuMP, Ipopt, Test

function example_rosekbrock(; verbose = true)
    model = Model(with_optimizer(Ipopt.Optimizer, print_level=0))
    @variable(model, x)
    @variable(model, y)
    @NLobjective(model, Min, (1 - x)^2 + 100 * (y - x^2)^2)
    JuMP.optimize!(model)
    if verbose
        println("x = ", JuMP.value(x), " y = ", JuMP.value(y))
    end
    @test JuMP.termination_status(model) == MOI.LOCALLY_SOLVED
    @test JuMP.primal_status(model) == MOI.FEASIBLE_POINT
    @test JuMP.objective_value(model) ≈ 0.0 atol = 1e-10
    @test JuMP.value(x) ≈ 1.0
    @test JuMP.value(y) ≈ 1.0
end

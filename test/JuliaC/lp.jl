#  Copyright 2026, Iain Dunning, Joey Huchette, Miles Lubin, and contributors
#  This Source Code Form is subject to the terms of the Mozilla Public
#  License, v. 2.0. If a copy of the MPL was not distributed with this
#  file, You can obtain one at https://mozilla.org/MPL/2.0/.
#############################################################################
# JuMP
# An algebraic modeling language for Julia
# See https://github.com/jump-dev/JuMP.jl
#############################################################################

using JuMP
import HiGHS

function @main(args::Vector{String})
    optimizer = HiGHS.Optimizer()
    model = JuMP.concrete_direct_model(optimizer)
    @assert backend(model) === optimizer
    set_silent(model)
    @variable(model, 0 <= x <= 1)
    @variable(model, y >= 0)
    @constraint(model, x + y <= 3)
    @objective(model, Max, 2x + y)
    optimize!(model)
    @assert termination_status(model) == MOI.OPTIMAL
    @assert primal_status(model) == MOI.FEASIBLE_POINT
    @assert isapprox(value(x), 1.0; atol = 1e-8)
    @assert isapprox(value(y), 2.0; atol = 1e-8)
    @assert isapprox(objective_value(model), 4.0; atol = 1e-8)
    return 0
end

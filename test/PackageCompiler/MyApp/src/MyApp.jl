#  Copyright 2017, Iain Dunning, Joey Huchette, Miles Lubin, and contributors
#  This Source Code Form is subject to the terms of the Mozilla Public
#  License, v. 2.0. If a copy of the MPL was not distributed with this
#  file, You can obtain one at https://mozilla.org/MPL/2.0/.
#############################################################################
# JuMP
# An algebraic modeling language for Julia
# See https://github.com/jump-dev/JuMP.jl
#############################################################################

module MyApp

using JuMP
import HiGHS
import Ipopt

function run_highs_example()
    capacity = 10.0
    profit = [5.0, 3.0, 2.0, 7.0, 4.0]
    weight = [2.0, 8.0, 4.0, 2.0, 5.0]
    model = Model(HiGHS.Optimizer)
    @variable(model, x[1:length(weight)], Bin)
    @constraint(model, weight' * x <= capacity)
    @objective(model, Max, profit' * x)
    optimize!(model)
    assert_is_solved_and_feasible(model)
    print(solution_summary(model))
    return
end

function run_ipopt_example()
    model = Model(Ipopt.Optimizer)
    @variable(model, x)
    @variable(model, y)
    @objective(model, Min, (1 - x)^2 + 100 * (y - x^2)^2)
    optimize!(model)
    assert_is_solved_and_feasible(model)
    print(solution_summary(model))
    return
end

function julia_main()::Cint
    solver = only(match(r"--solver=(.+)", only(ARGS)))
    if solver == "ipopt"
        run_ipopt_example()
    else
        @assert solver == "highs"
        run_highs_example()
    end
    return 0
end

end

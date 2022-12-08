#  Copyright 2017, Iain Dunning, Joey Huchette, Miles Lubin, and contributors
#  This Source Code Form is subject to the terms of the Mozilla Public
#  License, v. 2.0. If a copy of the MPL was not distributed with this
#  file, You can obtain one at https://mozilla.org/MPL/2.0/.

module TestSolutionSummary

using JuMP
using Test

function runtests()
    for name in names(@__MODULE__; all = true)
        if startswith("$(name)", "test_")
            @testset "$(name)" begin
                getfield(@__MODULE__, name)()
            end
        end
    end
    return
end

function test_empty_model()
    model = Model()
    @test sprint(show, solution_summary(model)) == """
* Solver : No optimizer attached.

* Status
  Result             : 1 of 0
  Termination status : OPTIMIZE_NOT_CALLED
  Primal status      : NO_SOLUTION
  Dual status        : NO_SOLUTION
  Message from the solver:
  "optimize not called"

* Candidate solution

* Work counters
"""
    return
end

function test_solution_summary()
    model = Model()
    @variable(model, x <= 2.0)
    @variable(model, y >= 0.0)
    @objective(model, Min, -x)
    c = @constraint(model, x + y <= 1) # anonymous constraint
    set_optimizer(
        model,
        () -> MOI.Utilities.MockOptimizer(
            MOI.Utilities.Model{Float64}();
            eval_objective_value = false,
        ),
    )
    optimize!(model)
    mock = JuMP.unsafe_backend(model)
    MOI.set(mock, MOI.TerminationStatus(), MOI.OPTIMAL)
    MOI.set(mock, MOI.RawStatusString(), "solver specific string")
    MOI.set(mock, MOI.ObjectiveValue(1), -1.0)
    MOI.set(mock, MOI.ObjectiveValue(2), -0.0)
    MOI.set(mock, MOI.ObjectiveBound(), -3.0)
    MOI.set(mock, MOI.RelativeGap(), abs((-3 - -1) / -3))
    MOI.set(mock, MOI.ResultCount(), 2)
    MOI.set(mock, MOI.PrimalStatus(1), MOI.FEASIBLE_POINT)
    MOI.set(mock, MOI.DualStatus(1), MOI.FEASIBLE_POINT)
    MOI.set(mock, MOI.VariablePrimal(1), JuMP.optimizer_index(x), 1.0)
    MOI.set(mock, MOI.VariablePrimal(1), JuMP.optimizer_index(y), 0.0)
    MOI.set(mock, MOI.ConstraintDual(1), JuMP.optimizer_index(c), -1.0)
    MOI.set(
        mock,
        MOI.ConstraintDual(),
        JuMP.optimizer_index(JuMP.UpperBoundRef(x)),
        0.0,
    )
    MOI.set(
        mock,
        MOI.ConstraintDual(),
        JuMP.optimizer_index(JuMP.LowerBoundRef(y)),
        1.0,
    )
    MOI.set(mock, MOI.PrimalStatus(2), MOI.FEASIBLE_POINT)
    MOI.set(mock, MOI.DualStatus(2), MOI.NO_SOLUTION)
    MOI.set(mock, MOI.VariablePrimal(2), JuMP.optimizer_index(x), 0.0)
    MOI.set(mock, MOI.VariablePrimal(2), JuMP.optimizer_index(y), 0.0)
    MOI.set(mock, MOI.SimplexIterations(), Int64(3))
    MOI.set(mock, MOI.BarrierIterations(), Int64(2))
    MOI.set(mock, MOI.NodeCount(), Int64(1))
    MOI.set(mock, MOI.SolveTimeSec(), 5.0)
    @test sprint(show, solution_summary(model)) == """
* Solver : Mock

* Status
  Result             : 1 of 2
  Termination status : OPTIMAL
  Primal status      : FEASIBLE_POINT
  Dual status        : FEASIBLE_POINT
  Message from the solver:
  "solver specific string"

* Candidate solution
  Objective value      : -1.00000e+00
  Objective bound      : -3.00000e+00
  Relative gap         : 6.66667e-01
  Dual objective value : -1.00000e+00

* Work counters
  Solve time (sec)   : 5.00000e+00
  Simplex iterations : 3
  Barrier iterations : 2
  Node count         : 1
"""

    summary = solution_summary(model; verbose = true)
    @test sprint(show, summary) == """
* Solver : Mock

* Status
  Result             : 1 of 2
  Termination status : OPTIMAL
  Primal status      : FEASIBLE_POINT
  Dual status        : FEASIBLE_POINT
  Message from the solver:
  "solver specific string"

* Candidate solution
  Objective value      : -1.00000e+00
  Objective bound      : -3.00000e+00
  Relative gap         : 6.66667e-01
  Dual objective value : -1.00000e+00
  Primal solution :
    x : 1.00000e+00
    y : 0.00000e+00
  Dual solution :

* Work counters
  Solve time (sec)   : 5.00000e+00
  Simplex iterations : 3
  Barrier iterations : 2
  Node count         : 1
"""

    summary = solution_summary(model; result = 2)
    @test sprint(show, summary) == """
* Solver : Mock

* Status
  Result             : 2 of 2
  Termination status : OPTIMAL
  Primal status      : FEASIBLE_POINT
  Dual status        : NO_SOLUTION

* Candidate solution
  Objective value      : -0.00000e+00

"""

    summary = solution_summary(model; result = 2, verbose = true)
    @test sprint(show, summary) == """
* Solver : Mock

* Status
  Result             : 2 of 2
  Termination status : OPTIMAL
  Primal status      : FEASIBLE_POINT
  Dual status        : NO_SOLUTION

* Candidate solution
  Objective value      : -0.00000e+00
  Primal solution :
    x : 0.00000e+00
    y : 0.00000e+00

"""

    summary = solution_summary(model; result = 3)
    @test sprint(show, summary) == """
* Solver : Mock

* Status
  Result             : 3 of 2
  Termination status : OPTIMAL
  Primal status      : NO_SOLUTION
  Dual status        : NO_SOLUTION

* Candidate solution

"""
    return
end

end

TestSolutionSummary.runtests()

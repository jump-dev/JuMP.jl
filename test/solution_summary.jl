#  Copyright 2017, Iain Dunning, Joey Huchette, Miles Lubin, and contributors
#  This Source Code Form is subject to the terms of the Mozilla Public
#  License, v. 2.0. If a copy of the MPL was not distributed with this
#  file, You can obtain one at https://mozilla.org/MPL/2.0/.

using JuMP
using Test

@testset "Print solution summary" begin
    model = Model()
    @variable(model, x <= 2.0)
    @variable(model, y >= 0.0)
    @objective(model, Min, -x)
    c = @constraint(model, x + y <= 1) # anonymous constraint
    set_optimizer(
        model,
        () -> MOI.Utilities.MockOptimizer(
            MOI.Utilities.Model{Float64}(),
            eval_objective_value = false,
        ),
    )
    optimize!(model)
    mockoptimizer = JuMP.unsafe_backend(model)
    MOI.set(mockoptimizer, MOI.TerminationStatus(), MOI.OPTIMAL)
    MOI.set(mockoptimizer, MOI.RawStatusString(), "solver specific string")
    MOI.set(mockoptimizer, MOI.ObjectiveValue(), -1.0)
    MOI.set(mockoptimizer, MOI.ObjectiveBound(), 3.0)
    MOI.set(mockoptimizer, MOI.ResultCount(), 1)
    MOI.set(mockoptimizer, MOI.PrimalStatus(), MOI.FEASIBLE_POINT)
    MOI.set(mockoptimizer, MOI.DualStatus(), MOI.FEASIBLE_POINT)
    MOI.set(mockoptimizer, MOI.VariablePrimal(), JuMP.optimizer_index(x), 1.0)
    MOI.set(mockoptimizer, MOI.VariablePrimal(), JuMP.optimizer_index(y), 0.0)
    MOI.set(mockoptimizer, MOI.ConstraintDual(), JuMP.optimizer_index(c), -1.0)
    MOI.set(
        mockoptimizer,
        MOI.ConstraintDual(),
        JuMP.optimizer_index(JuMP.UpperBoundRef(x)),
        0.0,
    )
    MOI.set(
        mockoptimizer,
        MOI.ConstraintDual(),
        JuMP.optimizer_index(JuMP.LowerBoundRef(y)),
        1.0,
    )
    MOI.set(mockoptimizer, MOI.SimplexIterations(), Int64(3))
    MOI.set(mockoptimizer, MOI.BarrierIterations(), Int64(2))
    MOI.set(mockoptimizer, MOI.NodeCount(), Int64(1))
    MOI.set(mockoptimizer, MOI.SolveTimeSec(), 5.0)
    @test sprint(show, solution_summary(model)) == """
* Solver : Mock

* Status
  Termination status : OPTIMAL
  Primal status      : FEASIBLE_POINT
  Dual status        : FEASIBLE_POINT
  Message from the solver:
  "solver specific string"

* Candidate solution
  Objective value      : -1.0
  Objective bound      : 3.0
  Dual objective value : -1.0

* Work counters
  Solve time (sec)   : 5.00000
  Simplex iterations : 3
  Barrier iterations : 2
  Node count         : 1
"""
    @test sprint(
        (io, model) -> show(io, solution_summary(model, verbose = true)),
        model,
    ) == """
* Solver : Mock

* Status
  Termination status : OPTIMAL
  Primal status      : FEASIBLE_POINT
  Dual status        : FEASIBLE_POINT
  Result count       : 1
  Has duals          : true
  Message from the solver:
  "solver specific string"

* Candidate solution
  Objective value      : -1.0
  Objective bound      : 3.0
  Dual objective value : -1.0
  Primal solution :
    x : 1.0
    y : 0.0
  Dual solution :

* Work counters
  Solve time (sec)   : 5.00000
  Simplex iterations : 3
  Barrier iterations : 2
  Node count         : 1
"""
end

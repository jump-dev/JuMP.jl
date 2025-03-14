#  Copyright 2017, Iain Dunning, Joey Huchette, Miles Lubin, and contributors
#  This Source Code Form is subject to the terms of the Mozilla Public
#  License, v. 2.0. If a copy of the MPL was not distributed with this
#  file, You can obtain one at https://mozilla.org/MPL/2.0/.

module TestSolutionSummary

using JuMP
using Test

function test_empty_model()
    model = Model()
    @test sprint(show, solution_summary(model)) == """
    Solution summary
    ├ solver_name          : No optimizer attached.
    ├ Solution quality
    │ ├ termination_status : OPTIMIZE_NOT_CALLED
    │ ├ result_count       : 0
    │ └ raw_status         : optimize not called
    └ Solution (; result = 1)
      ├ primal_status        : NO_SOLUTION
      └ dual_status          : NO_SOLUTION"""
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
    mock = unsafe_backend(model)
    MOI.set(mock, MOI.TerminationStatus(), MOI.OPTIMAL)
    MOI.set(mock, MOI.RawStatusString(), "solver specific string")
    MOI.set(mock, MOI.ObjectiveValue(1), -1.0)
    MOI.set(mock, MOI.ObjectiveValue(2), -0.0)
    MOI.set(mock, MOI.ObjectiveBound(), -3.0)
    MOI.set(mock, MOI.RelativeGap(), abs((-3 - -1) / -3))
    MOI.set(mock, MOI.ResultCount(), 2)
    MOI.set(mock, MOI.PrimalStatus(1), MOI.FEASIBLE_POINT)
    MOI.set(mock, MOI.DualStatus(1), MOI.FEASIBLE_POINT)
    MOI.set(mock, MOI.VariablePrimal(1), optimizer_index(x), 1.0)
    MOI.set(mock, MOI.VariablePrimal(1), optimizer_index(y), 0.0)
    MOI.set(mock, MOI.ConstraintDual(1), optimizer_index(c), -1.0)
    MOI.set(mock, MOI.ConstraintDual(), optimizer_index(UpperBoundRef(x)), 0.0)
    MOI.set(mock, MOI.ConstraintDual(), optimizer_index(LowerBoundRef(y)), 1.0)
    MOI.set(mock, MOI.PrimalStatus(2), MOI.FEASIBLE_POINT)
    MOI.set(mock, MOI.DualStatus(2), MOI.NO_SOLUTION)
    MOI.set(mock, MOI.VariablePrimal(2), optimizer_index(x), 0.0)
    MOI.set(mock, MOI.VariablePrimal(2), optimizer_index(y), 0.0)
    MOI.set(mock, MOI.SimplexIterations(), Int64(3))
    MOI.set(mock, MOI.BarrierIterations(), Int64(2))
    MOI.set(mock, MOI.NodeCount(), Int64(1))
    MOI.set(mock, MOI.SolveTimeSec(), 5.0)
    @test sprint(show, solution_summary(model)) == """
    Solution summary
    ├ solver_name          : Mock
    ├ Solution quality
    │ ├ termination_status : OPTIMAL
    │ ├ result_count       : 2
    │ ├ raw_status         : solver specific string
    │ ├ objective_bound    : -3.00000e+00
    │ └ relative_gap       : 6.66667e-01
    ├ Work counters
    │ ├ solve_time (sec)   : 5.00000e+00
    │ ├ simplex_iterations : 3
    │ ├ barrier_iterations : 2
    │ └ node_count         : 1
    └ Solution (; result = 1)
      ├ primal_status        : FEASIBLE_POINT
      ├ dual_status          : FEASIBLE_POINT
      ├ objective_value      : -1.00000e+00
      └ dual_objective_value : -1.00000e+00"""
    @test sprint(show, solution_summary(model; verbose = true)) == """
    Solution summary
    ├ solver_name          : Mock
    ├ Solution quality
    │ ├ termination_status : OPTIMAL
    │ ├ result_count       : 2
    │ ├ raw_status         : solver specific string
    │ ├ objective_bound    : -3.00000e+00
    │ └ relative_gap       : 6.66667e-01
    ├ Work counters
    │ ├ solve_time (sec)   : 5.00000e+00
    │ ├ simplex_iterations : 3
    │ ├ barrier_iterations : 2
    │ └ node_count         : 1
    └ Solution (; result = 1)
      ├ primal_status        : FEASIBLE_POINT
      ├ dual_status          : FEASIBLE_POINT
      ├ objective_value      : -1.00000e+00
      ├ dual_objective_value : -1.00000e+00
      └ value
        ├ x : 1.00000e+00
        └ y : 0.00000e+00"""
    @test sprint(show, solution_summary(model; result = 2)) == """
    Solution summary
    ├ solver_name          : Mock
    ├ Solution quality
    │ ├ termination_status : OPTIMAL
    │ ├ result_count       : 2
    │ ├ raw_status         : solver specific string
    │ ├ objective_bound    : -3.00000e+00
    │ └ relative_gap       : 6.66667e-01
    ├ Work counters
    │ ├ solve_time (sec)   : 5.00000e+00
    │ ├ simplex_iterations : 3
    │ ├ barrier_iterations : 2
    │ └ node_count         : 1
    └ Solution (; result = 2)
      ├ primal_status        : FEASIBLE_POINT
      ├ dual_status          : NO_SOLUTION
      └ objective_value      : -0.00000e+00"""
    suummary = solution_summary(model; result = 2, verbose = true)
    @test sprint(show, summary) == """
    Solution summary
    ├ solver_name          : Mock
    ├ Solution quality
    │ ├ termination_status : OPTIMAL
    │ ├ result_count       : 2
    │ ├ raw_status         : solver specific string
    │ ├ objective_bound    : -3.00000e+00
    │ └ relative_gap       : 6.66667e-01
    ├ Work counters
    │ ├ solve_time (sec)   : 5.00000e+00
    │ ├ simplex_iterations : 3
    │ ├ barrier_iterations : 2
    │ └ node_count         : 1
    └ Solution (; result = 2)
      ├ primal_status        : FEASIBLE_POINT
      ├ dual_status          : NO_SOLUTION
      ├ objective_value      : -0.00000e+00
      └ value
        ├ x : 0.00000e+00
        └ y : 0.00000e+00"""
    @test sprint(show, solution_summary(model; result = 3)) == """
    Solution summary
    ├ solver_name          : Mock
    ├ Solution quality
    │ ├ termination_status : OPTIMAL
    │ ├ result_count       : 2
    │ ├ raw_status         : solver specific string
    │ ├ objective_bound    : -3.00000e+00
    │ └ relative_gap       : 6.66667e-01
    ├ Work counters
    │ ├ solve_time (sec)   : 5.00000e+00
    │ ├ simplex_iterations : 3
    │ ├ barrier_iterations : 2
    │ └ node_count         : 1
    └ Solution (; result = 3)
      ├ primal_status        : NO_SOLUTION
      └ dual_status          : NO_SOLUTION"""
    return
end

function test_solution_summary_vector_dual()
    model = Model() do
        return MOI.Utilities.MockOptimizer(MOI.Utilities.Model{Float64}())
    end
    @variable(model, x[1:2])
    @constraint(model, c, x >= 0)
    optimize!(model)
    mock = unsafe_backend(model)
    MOI.set(mock, MOI.TerminationStatus(), MOI.OPTIMAL)
    MOI.set(mock, MOI.RawStatusString(), "solver specific string")
    MOI.set(mock, MOI.ResultCount(), 1)
    MOI.set(mock, MOI.PrimalStatus(), MOI.FEASIBLE_POINT)
    MOI.set(mock, MOI.DualStatus(), MOI.FEASIBLE_POINT)
    MOI.set(mock, MOI.VariablePrimal(), optimizer_index.(x), [1.0, 2.0])
    MOI.set(mock, MOI.ConstraintDual(), optimizer_index.(c), [3.0, 4.0])
    @test sprint(show, solution_summary(model; verbose = true)) == """
    Solution summary
    ├ solver_name          : Mock
    ├ Solution quality
    │ ├ termination_status : OPTIMAL
    │ ├ result_count       : 1
    │ └ raw_status         : solver specific string
    └ Solution (; result = 1)
      ├ primal_status        : FEASIBLE_POINT
      ├ dual_status          : FEASIBLE_POINT
      ├ objective_value      : 0.00000e+00
      ├ dual_objective_value : 0.00000e+00
      ├ value
      │ ├ x[1] : 1.00000e+00
      │ └ x[2] : 2.00000e+00
      └ dual
        └ c : [3.00000e+00,4.00000e+00]"""
    return
end

function test_solution_summary_same_names()
    model = Model() do
        return MOI.Utilities.MockOptimizer(MOI.Utilities.Model{Float64}())
    end
    @variable(model, x[1:2])
    @variable(model, y)
    z = @variable(model)
    set_name.(x, "x")
    @constraint(model, c, x .>= 0)
    @constraint(model, d, 2x[1] <= 1)
    e = @constraint(model, 2x[2] == 1)
    optimize!(model)
    mock = unsafe_backend(model)
    MOI.set(mock, MOI.TerminationStatus(), MOI.OPTIMAL)
    MOI.set(mock, MOI.RawStatusString(), "solver specific string")
    MOI.set(mock, MOI.ResultCount(), 1)
    MOI.set(mock, MOI.PrimalStatus(), MOI.FEASIBLE_POINT)
    MOI.set(mock, MOI.DualStatus(), MOI.FEASIBLE_POINT)
    MOI.set(mock, MOI.VariablePrimal(), optimizer_index.(x), [1.0, 2.0])
    MOI.set(mock, MOI.VariablePrimal(), optimizer_index(y), 3.0)
    MOI.set(mock, MOI.VariablePrimal(), optimizer_index(z), 4.0)
    MOI.set(mock, MOI.ConstraintDual(), optimizer_index.(c), [3.0, 4.0])
    MOI.set(mock, MOI.ConstraintDual(), optimizer_index(d), 5.0)
    MOI.set(mock, MOI.ConstraintDual(), optimizer_index(e), 6.0)
    @test sprint(show, solution_summary(model; verbose = true)) == """
    Solution summary
    ├ solver_name          : Mock
    ├ Solution quality
    │ ├ termination_status : OPTIMAL
    │ ├ result_count       : 1
    │ └ raw_status         : solver specific string
    └ Solution (; result = 1)
      ├ primal_status        : FEASIBLE_POINT
      ├ dual_status          : FEASIBLE_POINT
      ├ objective_value      : 0.00000e+00
      ├ dual_objective_value : 1.10000e+01
      ├ value
      │ ├ x : multiple variables with the same name
      │ └ y : 3.00000e+00
      └ dual
        ├ c : multiple constraints with the same name
        └ d : 5.00000e+00"""
    return
end

end  # module

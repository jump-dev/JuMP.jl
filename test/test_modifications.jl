#  Copyright 2017, Iain Dunning, Joey Huchette, Miles Lubin, and contributors
#  This Source Code Form is subject to the terms of the Mozilla Public
#  License, v. 2.0. If a copy of the MPL was not distributed with this
#  file, You can obtain one at https://mozilla.org/MPL/2.0/.

"""
    TestModifications

This module tests the modify-then-query behavior of JuMP models.

In order to simplify the testing, it frequently hacks the `model.is_model_dirty`
field of a JuMP model.
"""
module TestModifications

using JuMP
using Test

function test_variable_add_delete()
    model = Model()
    @test !model.is_model_dirty
    @variable(model, x)
    @test model.is_model_dirty
    model.is_model_dirty = false  # Hack!
    delete(model, x)
    @test model.is_model_dirty
    return
end

function test_variable_integer()
    model = Model()
    @variable(model, x)
    model.is_model_dirty = false  # Hack!
    set_integer(x)
    @test model.is_model_dirty
    model.is_model_dirty = false  # Hack!
    unset_integer(x)
    @test model.is_model_dirty
    return
end

function test_variable_binary()
    model = Model()
    @variable(model, x)
    model.is_model_dirty = false  # Hack!
    set_binary(x)
    @test model.is_model_dirty
    model.is_model_dirty = false  # Hack!
    unset_binary(x)
    @test model.is_model_dirty
    return
end

function test_variable_lower()
    model = Model()
    @variable(model, x)
    model.is_model_dirty = false  # Hack!
    set_lower_bound(x, 0)
    @test model.is_model_dirty
    model.is_model_dirty = false  # Hack!
    delete_lower_bound(x)
    @test model.is_model_dirty
    return
end

function test_variable_upper()
    model = Model()
    @variable(model, x)
    model.is_model_dirty = false  # Hack!
    set_upper_bound(x, 0)
    @test model.is_model_dirty
    model.is_model_dirty = false  # Hack!
    delete_upper_bound(x)
    @test model.is_model_dirty
    return
end

function test_variable_fix()
    model = Model()
    @variable(model, x)
    model.is_model_dirty = false  # Hack!
    fix(x, 0)
    @test model.is_model_dirty
    model.is_model_dirty = false  # Hack!
    unfix(x)
    @test model.is_model_dirty
    return
end

function test_constraint_add_delete()
    model = Model()
    @variable(model, x)
    model.is_model_dirty = false  # Hack!
    c = @constraint(model, x <= 1)
    @test model.is_model_dirty
    model.is_model_dirty = false  # Hack!
    delete(model, c)
    @test model.is_model_dirty
    return
end

function test_constraint_rhs()
    model = Model()
    @variable(model, x)
    c = @constraint(model, x <= 1)
    model.is_model_dirty = false  # Hack!
    set_normalized_rhs(c, 2.0)
    @test model.is_model_dirty
    return
end

function test_constraint_coef()
    model = Model()
    @variable(model, x)
    c = @constraint(model, x <= 1)
    model.is_model_dirty = false  # Hack!
    set_normalized_coefficient(c, x, 2.0)
    @test model.is_model_dirty
    return
end

function test_objective()
    model = Model()
    @variable(model, x)
    model.is_model_dirty = false  # Hack!
    @objective(model, Min, x)
    @test model.is_model_dirty
    model.is_model_dirty = false  # Hack!
    @objective(model, Max, -x)
    @test model.is_model_dirty
    return
end

function test_objective_sense()
    model = Model()
    @variable(model, x)
    @objective(model, Min, x)
    model.is_model_dirty = false  # Hack!
    set_objective_sense(model, MIN_SENSE)
    @test model.is_model_dirty
    return
end

function test_objective_function()
    model = Model()
    @variable(model, x)
    @objective(model, Min, x)
    model.is_model_dirty = false  # Hack!
    set_objective_function(model, 2 * x)
    @test model.is_model_dirty
    return
end

function test_objective_modify()
    model = Model()
    @variable(model, x)
    @objective(model, Min, x)
    model.is_model_dirty = false  # Hack!
    set_objective_coefficient(model, x, 2)
    @test model.is_model_dirty
    return
end

function test_status_caching()
    model = Model()
    @variable(model, x)
    @constraint(model, c, x == 1.2)
    @objective(model, Min, x - 0.1)
    set_optimizer(
        model,
        () -> MOI.Utilities.MockOptimizer(MOI.Utilities.Model{Float64}()),
    )
    optimize!(model)
    mock = unsafe_backend(model)
    MOI.set(mock, MOI.TerminationStatus(), MOI.OPTIMAL)
    MOI.set(mock, MOI.PrimalStatus(), MOI.FEASIBLE_POINT)
    MOI.set(mock, MOI.DualStatus(), MOI.FEASIBLE_POINT)
    MOI.set(mock, MOI.ResultCount(), 1)
    xis = MOI.get(mock, MOI.ListOfVariableIndices())
    MOI.set(mock, MOI.VariablePrimal(), xis[1], 1.2)
    cis = MOI.get(
        mock,
        MOI.ListOfConstraintIndices{
            MOI.ScalarAffineFunction{Float64},
            MOI.EqualTo{Float64},
        }(),
    )
    MOI.set(mock, MOI.ConstraintDual(), cis[1], 1.3)
    @test termination_status(model) == MOI.OPTIMAL
    @test dual_status(model) == MOI.FEASIBLE_POINT
    @test primal_status(model) == MOI.FEASIBLE_POINT
    @test objective_value(model) ≈ 1.1
    @test value(x) == 1.2
    @test dual(c) == 1.3
    model.is_model_dirty = true  # Hack!
    @test termination_status(model) == MOI.OPTIMIZE_NOT_CALLED
    @test dual_status(model) == MOI.NO_SOLUTION
    @test primal_status(model) == MOI.NO_SOLUTION
    @test_throws OptimizeNotCalled() value(x)
    @test_throws OptimizeNotCalled() dual(c)
    @test_throws OptimizeNotCalled() objective_value(model)
    return
end

function test_status_direct()
    mock = MOI.Utilities.MockOptimizer(MOI.Utilities.Model{Float64}())
    model = direct_model(mock)
    @variable(model, x)
    @constraint(model, c, x == 1.2)
    @objective(model, Min, x - 0.1)
    optimize!(model)
    MOI.set(mock, MOI.TerminationStatus(), MOI.OPTIMAL)
    MOI.set(mock, MOI.PrimalStatus(), MOI.FEASIBLE_POINT)
    MOI.set(mock, MOI.DualStatus(), MOI.FEASIBLE_POINT)
    MOI.set(mock, MOI.ResultCount(), 1)
    xis = MOI.get(mock, MOI.ListOfVariableIndices())
    MOI.set(mock, MOI.VariablePrimal(), xis[1], 1.2)
    cis = MOI.get(
        mock,
        MOI.ListOfConstraintIndices{
            MOI.ScalarAffineFunction{Float64},
            MOI.EqualTo{Float64},
        }(),
    )
    MOI.set(mock, MOI.ConstraintDual(), cis[1], 1.3)
    @test termination_status(model) == MOI.OPTIMAL
    @test dual_status(model) == MOI.FEASIBLE_POINT
    @test primal_status(model) == MOI.FEASIBLE_POINT
    @test objective_value(model) ≈ 1.1
    @test value(x) == 1.2
    @test dual(c) == 1.3
    model.is_model_dirty = true  # Hack!
    @test termination_status(model) == MOI.OPTIMAL
    @test dual_status(model) == MOI.FEASIBLE_POINT
    @test primal_status(model) == MOI.FEASIBLE_POINT
    @test objective_value(model) ≈ 1.1
    @test value(x) == 1.2
    @test dual(c) == 1.3
    return
end

end

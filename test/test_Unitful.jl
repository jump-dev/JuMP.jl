#  Copyright 2017, Iain Dunning, Joey Huchette, Miles Lubin, and contributors
#  This Source Code Form is subject to the terms of the Mozilla Public
#  License, v. 2.0. If a copy of the MPL was not distributed with this
#  file, You can obtain one at https://mozilla.org/MPL/2.0/.

module TestUnitful

using Test

using JuMP
using Unitful

include(joinpath(@__DIR__, "utilities.jl"))

function test_unitful_variable()
    model = Model()
    @variable(model, x, u"m/s")
    @test x isa GenericAffExpr{typeof(1.0u"m/s"),VariableRef}
    @test x == 1.0u"m/s" * only(all_variables(model))
    @test name(only(all_variables(model))) == "x"
    @test sprint(show, x) == "(x) m s^-1"
    return
end

function test_unitful_variable_bounds()
    model = Model()
    @variable(model, 0u"m/s" <= x <= 36u"km/hr", u"m/s", start = 18u"km/hr")
    v = only(all_variables(model))
    @test lower_bound(v) == 0
    @test upper_bound(v) ≈ 10
    @test start_value(v) ≈ 5
    @variable(model, y == 2u"km", u"m")
    @test fix_value(all_variables(model)[2]) ≈ 2000
    @test_throws Unitful.DimensionError @variable(model, z >= 1u"s", u"m")
    return
end

function test_unitful_variable_bounds_without_units()
    model = Model()
    @test_throws_runtime(
        ErrorException(
            """
            In `@variable(model, x >= 0, u"m")`:

            Unable to use `0` as a bound or start value because it has no units.

            Use a quantity with the same dimension as `m`.
            """,
        ),
        @variable(model, x >= 0, u"m"),
    )
    @test_throws ErrorException @variable(model, x <= 1, u"m")
    @test_throws ErrorException @variable(model, x == 1, u"m")
    @test_throws ErrorException @variable(model, x, u"m", start = 1)
    @test_throws ErrorException @variable(model, x >= 1, 10u"m")
    @test num_variables(model) == 0
    return
end

function test_unitful_variable_scale()
    model = Model()
    @variable(model, 0u"m/s" <= x <= 36u"km/hr", 10u"m/s")
    v = only(all_variables(model))
    @test x isa GenericAffExpr{typeof(1.0u"m/s"),VariableRef}
    @test x == 10.0u"m/s" * v
    @test sprint(show, x) == "(10 x) m s^-1"
    @test lower_bound(v) == 0
    @test upper_bound(v) ≈ 1
    @test value(_ -> 0.5, x) ≈ 5u"m/s"
    return
end

function test_unitful_variable_container()
    model = Model()
    @variable(model, x[1:2] >= 0u"m", u"m")
    @test x isa Vector{GenericAffExpr{typeof(1.0u"m"),VariableRef}}
    @test name.(all_variables(model)) == ["x[1]", "x[2]"]
    return
end

function test_unitful_variable_value_type()
    model = GenericModel{Float32}()
    @variable(model, x, u"m")
    @test x isa GenericAffExpr{typeof(1.0f0u"m"),GenericVariableRef{Float32}}
    @test 2u"cm" * only(all_variables(model)) isa
          GenericAffExpr{typeof(1.0f0u"cm"),GenericVariableRef{Float32}}
    return
end

function test_unitful_expression()
    model = Model()
    @variable(model, x, u"m")
    @variable(model, y, u"ft")
    @expression(model, aff, x + y + 1u"cm")
    @test aff isa GenericAffExpr{typeof(1.0u"m"),VariableRef}
    @test sprint(show, aff) == "(x + 0.3048 y + 0.01) m"
    @test value(v -> 2.0, aff) ≈ 2.6196u"m"
    @expression(model, quad, x * y)
    @test sprint(show, quad) == "(x*y) ft m"
    @test value(v -> 2.0, quad) ≈ 1.2192u"m^2"
    return
end

function test_unitful_constraint()
    model = Model()
    @variable(model, x, u"m/s")
    @variable(model, y, u"ft/s")
    v = all_variables(model)
    c = @constraint(model, x + y <= 36u"km/hr")
    @test isequal_canonical(constraint_object(c).func, v[1] + 0.3048 * v[2])
    @test constraint_object(c).set == MOI.LessThan(10.0)
    c = @constraint(model, x == y)
    @test isequal_canonical(constraint_object(c).func, v[1] - 0.3048 * v[2])
    @test constraint_object(c).set == MOI.EqualTo(0.0)
    c = @constraint(model, x * x >= 1u"m^2/s^2")
    @test isequal_canonical(constraint_object(c).func, v[1]^2)
    @test constraint_object(c).set == MOI.GreaterThan(1.0)
    @test_throws Unitful.DimensionError @constraint(model, x <= 1u"m")
    return
end

function test_unitful_constraint_with_unit()
    model = Model()
    @variable(model, x, u"m/s")
    v = only(all_variables(model))
    c = @constraint(model, x <= 10u"m/s", u"km/hr")
    @test isequal_canonical(constraint_object(c).func, 3.6 * v)
    @test constraint_object(c).set == MOI.LessThan(36.0)
    return
end

function test_unitful_interval()
    model = Model()
    @variable(model, x, u"m")
    v = only(all_variables(model))
    c = @constraint(model, 1u"cm" <= x <= 2u"m")
    @test isequal_canonical(constraint_object(c).func, 1.0 * v)
    @test constraint_object(c).set == MOI.Interval(0.01, 2.0)
    c = @constraint(model, 1u"cm" <= x <= 2u"m", u"cm")
    @test isequal_canonical(constraint_object(c).func, 100.0 * v)
    @test constraint_object(c).set == MOI.Interval(1.0, 200.0)
    return
end

function test_unitful_constraint_string()
    model = Model()
    @variable(model, x, u"m/s")
    @variable(model, y, u"ft/s")
    @constraint(model, c1, x + y <= 36u"km/hr")
    @test sprint(show, c1) == "c1 : (x + 0.3048 y) m s^-1 ≤ 10 m s^-1"
    @constraint(model, c2, x <= 36u"km/hr", u"km/hr")
    @test sprint(show, c2) == "c2 : (3.6 x) km hr^-1 ≤ 36 km hr^-1"
    @constraint(model, c3, 1u"cm/s" <= x <= 2u"m/s")
    @test sprint(show, c3) == "c3 : (x) m s^-1 ∈ [0.01, 2] m s^-1"
    return
end

function test_unitful_constraint_value_and_dual()
    model = Model() do
        return MOI.Utilities.MockOptimizer(MOI.Utilities.Model{Float64}())
    end
    @variable(model, x, u"m/s")
    @constraint(model, c1, x <= 36u"km/hr")
    @constraint(model, c2, x <= 36u"km/hr", u"km/hr")
    @constraint(model, c3, 0u"m/s" <= x <= 36u"km/hr", u"km/hr")
    MOI.Utilities.attach_optimizer(model)
    mock = unsafe_backend(model)
    MOI.set(mock, MOI.TerminationStatus(), MOI.OPTIMAL)
    MOI.set(mock, MOI.PrimalStatus(), MOI.FEASIBLE_POINT)
    MOI.set(mock, MOI.DualStatus(), MOI.FEASIBLE_POINT)
    v = optimizer_index(only(all_variables(model)))
    MOI.set(mock, MOI.VariablePrimal(), v, 10.0)
    MOI.set(mock, MOI.ConstraintDual(), optimizer_index(c1), -2.0)
    MOI.set(mock, MOI.ConstraintDual(), optimizer_index(c2), -2.0)
    MOI.set(mock, MOI.ConstraintDual(), optimizer_index(c3), -2.0)
    optimize!(model)
    @test value(c1) ≈ 10u"m/s"
    @test value(c1) isa typeof(1.0u"m/s")
    @test value(c2) ≈ 36u"km/hr"
    @test value(c2) isa typeof(1.0u"km/hr")
    @test value(c3) ≈ 36u"km/hr"
    @test value(_ -> 1.0, c1) ≈ 1u"m/s"
    @test dual(c1) ≈ -2u"s/m"
    @test dual(c1) isa typeof(1.0u"s/m")
    @test dual(c2) ≈ -2u"hr/km"
    @test dual(c2) isa typeof(1.0u"hr/km")
    @test dual(c3) ≈ -2u"hr/km"
    return
end

function test_unitful_constraint_start_value()
    model = Model()
    @variable(model, x, u"m/s")
    @constraint(model, c, x <= 36u"km/hr")
    @test start_value(c) === nothing
    @test dual_start_value(c) === nothing
    set_start_value(c, 18u"km/hr")
    @test start_value(c) ≈ 5u"m/s"
    @test MOI.get(model, MOI.ConstraintPrimalStart(), c) ≈ 5
    set_dual_start_value(c, -2u"s/m")
    @test dual_start_value(c) ≈ -2u"s/m"
    @test MOI.get(model, MOI.ConstraintDualStart(), c) ≈ -2
    set_start_value(c, nothing)
    @test start_value(c) === nothing
    set_dual_start_value(c, nothing)
    @test dual_start_value(c) === nothing
    @test_throws Unitful.DimensionError set_start_value(c, 1u"m")
    return
end

function test_unitful_constraint_ref()
    model = Model()
    @variable(model, x, u"m")
    @variable(model, y)
    @constraint(model, c1, x <= 1u"m")
    @constraint(model, c2, 2y <= 1)
    @test shape(constraint_object(c1)) == c1.shape
    @test c2.shape == ScalarShape()
    @test constraint_by_name(model, "c1") == c1
    @test typeof(constraint_by_name(model, "c1")) == typeof(c1)
    @test typeof(constraint_by_name(model, "c2")) == typeof(c2)
    @test all_constraints(model, AffExpr, MOI.LessThan{Float64}) == [c1, c2]
    new_model, ref_map = copy_model(model)
    new_c1 = constraint_by_name(new_model, "c1")
    @test new_c1.shape == c1.shape
    @test ref_map[c1] == new_c1
    return
end

function test_unitful_objective()
    model = Model()
    @variable(model, x, u"m")
    @variable(model, y, u"cm")
    v = all_variables(model)
    @objective(model, Min, x + y)
    @test isequal_canonical(objective_function(model), v[1] + 0.01 * v[2])
    @objective(model, Min, x * y)
    @test isequal_canonical(objective_function(model), v[1] * v[2])
    return
end

end

#  Copyright 2017, Iain Dunning, Joey Huchette, Miles Lubin, and contributors
#  This Source Code Form is subject to the terms of the Mozilla Public
#  License, v. 2.0. If a copy of the MPL was not distributed with this
#  file, You can obtain one at https://mozilla.org/MPL/2.0/.

module TestUnitful

using Test

using JuMP
using Unitful

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
    @variable(model, 0 <= x <= 36u"km/hr", u"m/s", start = 18u"km/hr")
    v = only(all_variables(model))
    @test lower_bound(v) == 0
    @test upper_bound(v) ≈ 10
    @test start_value(v) ≈ 5
    @test_throws Unitful.DimensionError @variable(model, y >= 1u"s", u"m")
    return
end

function test_unitful_variable_container()
    model = Model()
    @variable(model, x[1:2] >= 0, u"m")
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

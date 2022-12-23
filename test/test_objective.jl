#  Copyright 2017, Iain Dunning, Joey Huchette, Miles Lubin, and contributors
#  This Source Code Form is subject to the terms of the Mozilla Public
#  License, v. 2.0. If a copy of the MPL was not distributed with this
#  file, You can obtain one at https://mozilla.org/MPL/2.0/.
#############################################################################
# JuMP
# An algebraic modeling language for Julia
# See https://github.com/jump-dev/JuMP.jl
#############################################################################

module TestObjective

using JuMP
using Test

struct DummyOptimizer <: MOI.AbstractOptimizer end

MOI.is_empty(::DummyOptimizer) = true

function test_unsupported_objective_function()
    model = Model(DummyOptimizer)
    func = MOI.VariableIndex(1)
    @test_throws ErrorException set_objective_function(model, func)
    return
end

function test_unsupported_function_in_macro()
    model = Model()
    @variable(model, x[1:2, 1:2])
    @test_throws(
        ErrorException("The objective function `$x` is not supported by JuMP."),
        @objective(model, Min, x),
    )
    return
end

function test_objective_coef_update_linear_objective_changes()
    model = Model()
    @variable(model, x)
    @objective(model, Max, x)
    set_objective_coefficient(model, x, 4.0)
    @test isequal_canonical(objective_function(model), 4x)
    @variable(model, y)
    @objective(model, Max, x + y)
    set_objective_coefficient(model, x, 4.0)
    @test isequal_canonical(objective_function(model), 4x + y)
    @objective(model, Min, x)
    set_objective_coefficient(model, y, 2.0)
    @test isequal_canonical(objective_function(model), x + 2.0 * y)
    return
end

function test_objective_coef_update_quadratic_objective_changes()
    model = Model()
    @variable(model, x)
    @objective(model, Max, x^2 + x)
    set_objective_coefficient(model, x, 4.0)
    @test isequal_canonical(objective_function(model), x^2 + 4x)
    return
end

function test_extension_objective_sense_get_set(
    ModelType = Model,
    VariableType = VariableRef,
)
    model = ModelType()
    set_objective_sense(model, FEASIBILITY_SENSE)
    @test FEASIBILITY_SENSE == @inferred objective_sense(model)
    return
end

function test_extension_objective_variable(
    ModelType = Model,
    VariableType = VariableRef,
)
    model = ModelType()
    @variable(model, x)
    @objective(model, Min, x)
    @test MIN_SENSE == @inferred objective_sense(model)
    @test objective_function_type(model) == VariableType
    @test objective_function(model) == x
    @test x == @inferred objective_function(model, VariableType)
    @objective(model, Max, x)
    @test MAX_SENSE == @inferred objective_sense(model)
    @test objective_function_type(model) == VariableType
    @test objective_function(model) == x
    @test x == @inferred objective_function(model, VariableType)
    return
end

function test_extension_objective_affine(
    ModelType = Model,
    VariableType = VariableRef,
)
    T = value_type(ModelType)
    model = ModelType()
    @variable(model, x)
    @objective(model, Min, 2x)
    @test MIN_SENSE == @inferred objective_sense(model)
    @test objective_function_type(model) == GenericAffExpr{T,VariableType}
    @test isequal_canonical(objective_function(model), 2x)
    @test isequal_canonical(
        2x,
        @inferred objective_function(model, GenericAffExpr{T,VariableType})
    )
    @objective(model, Max, x + 3x + 1)
    @test MAX_SENSE == @inferred objective_sense(model)
    @test objective_function_type(model) == GenericAffExpr{T,VariableType}
    @test isequal_canonical(objective_function(model), 4x + 1)
    @test isequal_canonical(
        4x + 1,
        @inferred objective_function(model, GenericAffExpr{T,VariableType})
    )
    return
end

function test_extension_objective_quadratic(
    ModelType = Model,
    VariableType = VariableRef,
)
    T = value_type(ModelType)
    model = ModelType()
    @variable(model, x)
    @objective(model, Min, x^2 + 2x)
    @test MIN_SENSE == @inferred objective_sense(model)
    @test objective_function_type(model) == GenericQuadExpr{T,VariableType}
    @test isequal_canonical(objective_function(model), x^2 + 2x)
    @test isequal_canonical(
        x^2 + 2x,
        @inferred objective_function(model, GenericQuadExpr{T,VariableType})
    )
    @test_throws InexactError objective_function(
        model,
        GenericAffExpr{T,VariableType},
    )
    return
end

function test_extension_objective_sense_as_symbol(
    ModelType = Model,
    VariableType = VariableRef,
)
    model = ModelType()
    @variable(model, x)
    @test_throws ErrorException @objective(model, :Min, 2x)
    return
end

function test_extension_objective_sense_as_binnding(
    ModelType = Model,
    VariableType = VariableRef,
)
    T = value_type(ModelType)
    model = ModelType()
    @variable(model, x)
    sense = MIN_SENSE
    @objective(model, sense, 2x)
    @test MIN_SENSE == @inferred objective_sense(model)
    @test isequal_canonical(
        2x,
        @inferred objective_function(model, GenericAffExpr{T,VariableType})
    )
    sense = :Min
    @test_throws ErrorException @objective(model, sense, 2x)
    return
end

function test_extension_objective_constant(
    ModelType = Model,
    VariableType = VariableRef,
)
    model = ModelType()
    @objective(model, Min, 3)
    @test objective_sense(model) == MIN_SENSE
    @test isequal_canonical(
        GenericAffExpr{Float64,VariableType}(3.0),
        objective_function(model, GenericAffExpr{Float64,VariableType}),
    )
    return
end

function test_extension_objective_vector_of_variables(
    ModelType = Model,
    VariableType = VariableRef,
)
    model = ModelType()
    @variable(model, x[1:2])
    @objective(model, Min, x)
    @test isequal_canonical(objective_function(model), x)
    @test isequal_canonical(objective_function(model, typeof(x)), x)
    return
end

function test_extension_objective_vector_affine_function(
    ModelType = Model,
    VariableType = VariableRef,
)
    model = ModelType()
    @variable(model, x[1:2])
    f = 1.0 .* x .+ [2.0, 3.0]
    @objective(model, Min, f)
    @test isequal_canonical(objective_function(model), f)
    @test isequal_canonical(objective_function(model, typeof(f)), f)
    return
end

function test_extension_objective_vector_quadratic_function(
    ModelType = Model,
    VariableType = VariableRef,
)
    model = ModelType()
    @variable(model, x[1:2])
    f = 1.0 .* x .* x .+ [2.0, 3.0] .* x + [4.0, 5.0]
    @objective(model, Min, f)
    @test isequal_canonical(objective_function(model), f)
    @test isequal_canonical(objective_function(model, typeof(f)), f)
    return
end

end  # module

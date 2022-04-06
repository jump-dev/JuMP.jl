#  Copyright 2017, Iain Dunning, Joey Huchette, Miles Lubin, and contributors
#  This Source Code Form is subject to the terms of the Mozilla Public
#  License, v. 2.0. If a copy of the MPL was not distributed with this
#  file, You can obtain one at https://mozilla.org/MPL/2.0/.
#############################################################################
# JuMP
# An algebraic modeling language for Julia
# See https://github.com/jump-dev/JuMP.jl
#############################################################################
# This file contains objective-related functions

"""
    relative_gap(model::Model)

Return the final relative optimality gap after a call to `optimize!(model)`.
Exact value depends upon implementation of MathOptInterface.RelativeGap()
by the particular solver used for optimization.
"""
function relative_gap(model::Model)::Float64
    return MOI.get(model, MOI.RelativeGap())
end

"""
    objective_bound(model::Model)

Return the best known bound on the optimal objective value after a call to
`optimize!(model)`.
"""
function objective_bound(model::Model)::Float64
    return MOI.get(model, MOI.ObjectiveBound())
end

"""
    objective_value(model::Model; result::Int = 1)

Return the objective value associated with result index `result` of the
most-recent solution returned by the solver.

See also: [`result_count`](@ref).
"""
function objective_value(model::Model; result::Int = 1)::Float64
    return MOI.get(model, MOI.ObjectiveValue(result))
end

"""
    dual_objective_value(model::Model; result::Int = 1)

Return the value of the objective of the dual problem associated with result
index `result` of the most-recent solution returned by the solver.

Throws `MOI.UnsupportedAttribute{MOI.DualObjectiveValue}` if the solver does
not support this attribute.

See also: [`result_count`](@ref).
"""
function dual_objective_value(model::Model; result::Int = 1)::Float64
    return MOI.get(model, MOI.DualObjectiveValue(result))
end

"""
    objective_sense(model::Model)::MOI.OptimizationSense

Return the objective sense.
"""
function objective_sense(model::Model)
    return MOI.get(model, MOI.ObjectiveSense())::MOI.OptimizationSense
end

"""
    set_objective_sense(model::Model, sense::MOI.OptimizationSense)

Sets the objective sense of the model to the given sense. See
[`set_objective_function`](@ref) to set the objective function. These are
low-level functions; the recommended way to set the objective is with the
[`@objective`](@ref) macro.
"""
function set_objective_sense(model::Model, sense::MOI.OptimizationSense)
    return MOI.set(model, MOI.ObjectiveSense(), sense)
end

"""
    set_objective_function(
        model::Model,
        func::Union{AbstractJuMPScalar, MathOptInterface.AbstractScalarFunction})

Sets the objective function of the model to the given function. See
[`set_objective_sense`](@ref) to set the objective sense. These are low-level
functions; the recommended way to set the objective is with the
[`@objective`](@ref) macro.
"""
function set_objective_function end

function set_objective_function(model::Model, func::MOI.AbstractScalarFunction)
    attr = MOI.ObjectiveFunction{typeof(func)}()
    if !MOI.supports(backend(model), attr)
        error(
            "The solver does not support an objective function of type ",
            typeof(func),
            ".",
        )
    end
    MOI.set(model, attr, func)
    # Nonlinear objectives override regular objectives, so if there was a
    # nonlinear objective set, we must clear it.
    if model.nlp_data !== nothing
        model.nlp_data.objective = nothing
    end
    return
end

function set_objective_function(model::Model, func::AbstractJuMPScalar)
    check_belongs_to_model(func, model)
    return set_objective_function(model, moi_function(func))
end

function set_objective_function(model::Model, func::Real)
    return set_objective_function(
        model,
        MOI.ScalarAffineFunction(
            MOI.ScalarAffineTerm{Float64}[],
            Float64(func),
        ),
    )
end

function set_objective_function(model::AbstractModel, func)
    return error("The objective function `$(func)` is not supported by JuMP.")
end

"""
    set_objective(model::AbstractModel, sense::MOI.OptimizationSense, func)

The functional equivalent of the [`@objective`](@ref) macro.

Sets the objective sense and objective function simultaneously, and is
equivalent to:
```julia
set_objective_sense(model, sense)
set_objective_function(model, func)
```

## Examples

```jldoctest; setup=:(using JuMP)
model = Model()
@variable(model, x)
set_objective(model, MIN_SENSE, x)
```
"""
function set_objective(model::AbstractModel, sense::MOI.OptimizationSense, func)
    set_objective_sense(model, sense)
    return set_objective_function(model, func)
end

"""
    objective_function_type(model::Model)::AbstractJuMPScalar

Return the type of the objective function.
"""
function objective_function_type(model::Model)
    return jump_function_type(
        model,
        MOI.get(backend(model), MOI.ObjectiveFunctionType()),
    )
end

"""
    objective_function(model::Model,
                   T::Type{<:AbstractJuMPScalar}=objective_function_type(model))

Return an object of type `T` representing the objective function.
Error if the objective is not convertible to type `T`.

## Examples

```jldoctest objective_function; setup = :(using JuMP)
julia> model = Model()
A JuMP Model
Feasibility problem with:
Variables: 0
Model mode: AUTOMATIC
CachingOptimizer state: NO_OPTIMIZER
Solver name: No optimizer attached.

julia> @variable(model, x)
x

julia> @objective(model, Min, 2x + 1)
2 x + 1

julia> objective_function(model, AffExpr)
2 x + 1

julia> objective_function(model, QuadExpr)
2 x + 1

julia> typeof(objective_function(model, QuadExpr))
GenericQuadExpr{Float64,VariableRef}
```
We see with the last two commands that even if the objective function is affine,
as it is convertible to a quadratic function, it can be queried as a quadratic
function and the result is quadratic.

However, it is not convertible to a variable.
```jldoctest objective_function; filter = r"MathOptInterface\\."s
julia> objective_function(model, VariableRef)
ERROR: InexactError: convert(MathOptInterface.VariableIndex, MathOptInterface.ScalarAffineFunction{Float64}(MathOptInterface.ScalarAffineTerm{Float64}[MathOptInterface.ScalarAffineTerm{Float64}(2.0, MathOptInterface.VariableIndex(1))], 1.0))
[...]
```
"""
function objective_function(
    model::Model,
    FunType::Type{<:AbstractJuMPScalar} = objective_function_type(model),
)
    MOIFunType = moi_function_type(FunType)
    func =
        MOI.get(backend(model), MOI.ObjectiveFunction{MOIFunType}())::MOIFunType
    return jump_function(model, func)
end

"""
    set_objective_coefficient(model::Model, variable::VariableRef, coefficient::Real)

Set the linear objective coefficient associated with `Variable` to `coefficient`.

Note: this function will throw an error if a nonlinear objective is set.
"""
function set_objective_coefficient(
    model::Model,
    variable::VariableRef,
    coeff::Real,
)
    if model.nlp_data !== nothing && _nlp_objective_function(model) !== nothing
        error("A nonlinear objective is already set in the model")
    end
    coeff = convert(Float64, coeff)::Float64
    obj_fct_type = objective_function_type(model)
    if obj_fct_type == VariableRef
        # Promote the objective function to be an affine expression.
        current_obj = objective_function(model)
        if index(current_obj) == index(variable)
            set_objective_function(model, coeff * variable)
        else
            set_objective_function(
                model,
                add_to_expression!(coeff * variable, current_obj),
            )
        end
    elseif obj_fct_type == AffExpr || obj_fct_type == QuadExpr
        MOI.modify(
            backend(model),
            MOI.ObjectiveFunction{moi_function_type(obj_fct_type)}(),
            MOI.ScalarCoefficientChange(index(variable), coeff),
        )
    else
        error("Objective function type not supported: $(obj_fct_type)")
    end
    model.is_model_dirty = true
    return
end

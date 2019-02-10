#  Copyright 2017, Iain Dunning, Joey Huchette, Miles Lubin, and contributors
#  This Source Code Form is subject to the terms of the Mozilla Public
#  License, v. 2.0. If a copy of the MPL was not distributed with this
#  file, You can obtain one at http://mozilla.org/MPL/2.0/.
#############################################################################
# JuMP
# An algebraic modeling language for Julia
# See http://github.com/JuliaOpt/JuMP.jl
#############################################################################
# This file contains objective-related functions

"""
    objective_bound(model::Model)

Return the best known bound on the optimal objective value after a call to
`optimize!(model)`.
"""
objective_bound(model::Model)::Float64 = MOI.get(model, MOI.ObjectiveBound())

"""
    objective_value(model::Model)

Return the objective value after a call to `optimize!(model)`.
"""
objective_value(model::Model)::Float64 = MOI.get(model, MOI.ObjectiveValue())

"""
    objective_sense(model::Model)::MathOptInterface.OptimizationSense

Return the objective sense.
"""
function objective_sense(model::Model)
    return MOI.get(model, MOI.ObjectiveSense())::MOI.OptimizationSense
end

"""
    set_objective_sense(model::Model, sense::MathOptInterface.OptimizationSense)

Sets the objective sense of the model to the given sense. See
[`set_objective_function`](@ref) to set the objective function. These are
low-level functions; the recommended way to set the objective is with the
[`@objective`](@ref) macro.
"""
function set_objective_sense(model::Model, sense::MOI.OptimizationSense)
    MOI.set(model, MOI.ObjectiveSense(), sense)
end

"""
    set_objective_function(model::Model,
                           func::Union{AbstractJuMPScalar,
                                       MathOptInterface.AbstractScalarFunction})

Sets the objective function of the model to the given function. See
[`set_objective_sense`](@ref) to set the objective sense. These are low-level
functions; the recommended way to set the objective is with the
[`@objective`](@ref) macro.
"""
function set_objective_function end

function set_objective_function(model::Model, func::MOI.AbstractScalarFunction)
    attr = MOI.ObjectiveFunction{typeof(func)}()
    if !MOI.supports(backend(model), attr)
        error("The solver does not support an objective function of type ",
              typeof(func), ".")
    end
    MOI.set(model, attr, func)
    # Nonlinear objectives override regular objectives, so if there was a
    # nonlinear objective set, we must clear it.
    if model.nlp_data !== nothing
        model.nlp_data.nlobj = nothing
    end
    return
end

function set_objective_function(model::Model, func::AbstractJuMPScalar)
    check_belongs_to_model(func, model)
    set_objective_function(model, moi_function(func))
end

function set_objective(model::Model, sense::MOI.OptimizationSense,
                       func::AbstractJuMPScalar)
    set_objective_sense(model, sense)
    set_objective_function(model, func)
end

function set_objective(model::Model, sense::MOI.OptimizationSense, func)
    error("The objective function `$(func)` is not supported by JuMP.")
end

"""
    objective_function_type(model::Model)::AbstractJuMPScalar

Return the type of the objective function.
"""
function objective_function_type(model::Model)
    jump_function_type(model,
                       MOI.get(backend(model), MOI.ObjectiveFunctionType()))
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
```jldoctest objective_function; filter = r"Stacktrace:.*"s
julia> objective_function(model, VariableRef)
ERROR: InexactError: convert(MathOptInterface.SingleVariable, MathOptInterface.ScalarAffineFunction{Float64}(MathOptInterface.ScalarAffineTerm{Float64}[ScalarAffineTerm{Float64}(2.0, VariableIndex(1))], 1.0))
Stacktrace:
 [1] convert at /home/blegat/.julia/dev/MathOptInterface/src/functions.jl:398 [inlined]
 [2] get(::JuMP.JuMPMOIModel{Float64}, ::MathOptInterface.ObjectiveFunction{MathOptInterface.SingleVariable}) at /home/blegat/.julia/dev/MathOptInterface/src/Utilities/model.jl:290
 [3] get at /home/blegat/.julia/dev/MathOptInterface/src/Utilities/universalfallback.jl:114 [inlined]
 [4] get at /home/blegat/.julia/dev/MathOptInterface/src/Utilities/cachingoptimizer.jl:439 [inlined]
 [5] get(::MathOptInterface.Bridges.LazyBridgeOptimizer{MathOptInterface.Utilities.CachingOptimizer{MathOptInterface.AbstractOptimizer,MathOptInterface.Utilities.UniversalFallback{JuMP.JuMPMOIModel{Float64}}},MathOptInterface.Bridges.AllBridgedConstraints{Float64}}, ::MathOptInterface.ObjectiveFunction{MathOptInterface.SingleVariable}) at /home/blegat/.julia/dev/MathOptInterface/src/Bridges/bridgeoptimizer.jl:172
 [6] objective_function(::Model, ::Type{VariableRef}) at /home/blegat/.julia/dev/JuMP/src/objective.jl:129
 [7] top-level scope at none:0
```
"""
function objective_function(model::Model,
             FunType::Type{<:AbstractJuMPScalar}=objective_function_type(model))
    MOIFunType = moi_function_type(FunType)
    func = MOI.get(backend(model),
                   MOI.ObjectiveFunction{MOIFunType}())::MOIFunType
    return jump_function(model, func)
end

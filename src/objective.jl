#  Copyright 2017, Iain Dunning, Joey Huchette, Miles Lubin, and contributors
#  This Source Code Form is subject to the terms of the Mozilla Public
#  License, v. 2.0. If a copy of the MPL was not distributed with this
#  file, You can obtain one at https://mozilla.org/MPL/2.0/.
#############################################################################
# JuMP
# An algebraic modeling language for Julia
# See https://github.com/jump-dev/JuMP.jl
#############################################################################

"""
    relative_gap(model::GenericModel)

Return the final relative optimality gap after a call to `optimize!(model)`.

Exact value depends upon implementation of [`MOI.RelativeGap`](@ref) by the
particular solver used for optimization.

This function is equivalent to querying the [`MOI.RelativeGap`](@ref) attribute.

## Example

```jldoctest
julia> import HiGHS

julia> model = Model(HiGHS.Optimizer);

julia> set_silent(model)

julia> @variable(model, x >= 1, Int);

julia> @objective(model, Min, 2 * x + 1);

julia> optimize!(model)

julia> relative_gap(model)
0.0
```
"""
function relative_gap(model::GenericModel{T})::T where {T}
    return MOI.get(model, MOI.RelativeGap())
end

"""
    objective_bound(model::GenericModel)

Return the best known bound on the optimal objective value after a call to
`optimize!(model)`.

For scalar-valued objectives, this function returns a `Float64`. For
vector-valued objectives, it returns a `Vector{Float64}`.

In the case of a vector-valued objective, this returns the _ideal point_, that
is, the point obtained if each objective was optimized independently.

This function is equivalent to querying the [`MOI.ObjectiveBound`](@ref) attribute.

## Example

```jldoctest
julia> import HiGHS

julia> model = Model(HiGHS.Optimizer);

julia> set_silent(model)

julia> @variable(model, x >= 1, Int);

julia> @objective(model, Min, 2 * x + 1);

julia> optimize!(model)

julia> objective_bound(model)
3.0
```
"""
function objective_bound(model::GenericModel{T})::Union{T,Vector{T}} where {T}
    return MOI.get(model, MOI.ObjectiveBound())
end

"""
    objective_value(model::GenericModel; result::Int = 1)

Return the objective value associated with result index `result` of the
most-recent solution returned by the solver.

For scalar-valued objectives, this function returns a `Float64`. For
vector-valued objectives, it returns a `Vector{Float64}`.

This function is equivalent to querying the [`MOI.ObjectiveValue`](@ref) attribute.

See also: [`result_count`](@ref).

## Example

```jldoctest
julia> import HiGHS

julia> model = Model(HiGHS.Optimizer);

julia> set_silent(model)

julia> @variable(model, x >= 1);

julia> @objective(model, Min, 2 * x + 1);

julia> optimize!(model)

julia> objective_value(model)
3.0

julia> objective_value(model; result = 2)
ERROR: Result index of attribute MathOptInterface.ObjectiveValue(2) out of bounds. There are currently 1 solution(s) in the model.
Stacktrace:
[...]
```
"""
function objective_value(
    model::GenericModel{T};
    result::Int = 1,
)::Union{T,Vector{T}} where {T}
    return MOI.get(model, MOI.ObjectiveValue(result))
end

"""
    dual_objective_value(model::GenericModel; result::Int = 1)

Return the value of the objective of the dual problem associated with result
index `result` of the most-recent solution returned by the solver.

Throws `MOI.UnsupportedAttribute{MOI.DualObjectiveValue}` if the solver does
not support this attribute.

This function is equivalent to querying the [`MOI.DualObjectiveValue`](@ref)
attribute.

See also: [`result_count`](@ref).

## Example

```jldoctest
julia> import HiGHS

julia> model = Model(HiGHS.Optimizer);

julia> set_silent(model)

julia> @variable(model, x >= 1);

julia> @objective(model, Min, 2 * x + 1);

julia> optimize!(model)

julia> dual_objective_value(model)
3.0

julia> dual_objective_value(model; result = 2)
ERROR: Result index of attribute MathOptInterface.DualObjectiveValue(2) out of bounds. There are currently 1 solution(s) in the model.
Stacktrace:
[...]
```
"""
function dual_objective_value(
    model::GenericModel{T};
    result::Int = 1,
)::T where {T}
    return MOI.get(model, MOI.DualObjectiveValue(result))
end

"""
    objective_sense(model::GenericModel)::MOI.OptimizationSense

Return the objective sense.

This function is equivalent to querying the [`MOI.ObjectiveSense`](@ref) attribute.

## Example

```jldoctest
julia> model = Model();

julia> objective_sense(model)
FEASIBILITY_SENSE::OptimizationSense = 2

julia> @variable(model, x);

julia> @objective(model, Max, x)
x

julia> objective_sense(model)
MAX_SENSE::OptimizationSense = 1
```
"""
function objective_sense(model::GenericModel)
    return MOI.get(model, MOI.ObjectiveSense())::MOI.OptimizationSense
end

"""
    set_objective_sense(model::GenericModel, sense::MOI.OptimizationSense)

Sets the objective sense of the model to the given sense.

See [`set_objective_function`](@ref) to set the objective function.

These are low-level functions; the recommended way to set the objective is with
the [`@objective`](@ref) macro.

## Example

```jldoctest
julia> model = Model();

julia> objective_sense(model)
FEASIBILITY_SENSE::OptimizationSense = 2

julia> set_objective_sense(model, MOI.MAX_SENSE)

julia> objective_sense(model)
MAX_SENSE::OptimizationSense = 1
```
"""
function set_objective_sense(model::GenericModel, sense::MOI.OptimizationSense)
    return MOI.set(model, MOI.ObjectiveSense(), sense)
end

"""
    set_objective_function(model::GenericModel, func::MOI.AbstractFunction)
    set_objective_function(model::GenericModel, func::AbstractJuMPScalar)
    set_objective_function(model::GenericModel, func::Real)
    set_objective_function(model::GenericModel, func::Vector{<:AbstractJuMPScalar})

Sets the objective function of the model to the given function.

See [`set_objective_sense`](@ref) to set the objective sense.

These are low-level functions; the recommended way to set the objective is with
the [`@objective`](@ref) macro.

## Example

```jldoctest
julia> model = Model();

julia> @variable(model, x);

julia> @objective(model, Min, x);

julia> objective_function(model)
x

julia> set_objective_function(model, 2 * x + 1)

julia> objective_function(model)
2 x + 1
```
"""
function set_objective_function end

function set_objective_function(model::GenericModel, func::MOI.AbstractFunction)
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
    nlp = nonlinear_model(model)
    if nlp !== nothing
        MOI.Nonlinear.set_objective(nlp, nothing)
    end
    return
end

function set_objective_function(model::GenericModel, func::AbstractJuMPScalar)
    check_belongs_to_model(func, model)
    return set_objective_function(model, moi_function(func))
end

function set_objective_function(model::GenericModel{T}, func::Real) where {T}
    return set_objective_function(
        model,
        MOI.ScalarAffineFunction(MOI.ScalarAffineTerm{T}[], convert(T, func)),
    )
end

function set_objective_function(
    model::GenericModel,
    func::AbstractVector{<:AbstractJuMPScalar},
)
    for f in func
        check_belongs_to_model(f, model)
    end
    return set_objective_function(model, moi_function(func))
end

function set_objective_function(model::AbstractModel, func)
    return error("The objective function `$(func)` is not supported by JuMP.")
end

"""
    set_objective(model::AbstractModel, sense::MOI.OptimizationSense, func)

The functional equivalent of the [`@objective`](@ref) macro.

Sets the objective sense and objective function simultaneously, and is
equivalent to calling [`set_objective_sense`](@ref) and
[`set_objective_function`](@ref) separately.

## Example

```jldoctest
julia> model = Model();

julia> @variable(model, x)
x

julia> set_objective(model, MIN_SENSE, x)
```
"""
function set_objective(model::AbstractModel, sense::MOI.OptimizationSense, func)
    set_objective_sense(model, sense)
    return set_objective_function(model, func)
end

"""
    objective_function_type(model::GenericModel)::AbstractJuMPScalar

Return the type of the objective function.

This function is equivalent to querying the [`MOI.ObjectiveFunctionType`](@ref)
attribute.

## Example

```jldoctest
julia> model = Model();

julia> @variable(model, x);

julia> @objective(model, Min, 2 * x + 1);

julia> objective_function_type(model)
AffExpr (alias for GenericAffExpr{Float64, GenericVariableRef{Float64}})
```
"""
function objective_function_type(model::GenericModel)
    return jump_function_type(
        model,
        MOI.get(backend(model), MOI.ObjectiveFunctionType()),
    )
end

"""
    objective_function(
        model::GenericModel,
        ::Type{F} = objective_function_type(model),
    ) where {F}

Return an object of type `F` representing the objective function.

Errors if the objective is not convertible to type `F`.

This function is equivalent to querying the [`MOI.ObjectiveFunction{F}`](@ref)
attribute.

## Example

```jldoctest objective_function
julia> model = Model();

julia> @variable(model, x)
x

julia> @objective(model, Min, 2x + 1)
2 x + 1

julia> objective_function(model, AffExpr)
2 x + 1

julia> objective_function(model, QuadExpr)
2 x + 1

julia> typeof(objective_function(model, QuadExpr))
QuadExpr (alias for GenericQuadExpr{Float64, GenericVariableRef{Float64}})
```

We see with the last two commands that even if the objective function is affine,
as it is convertible to a quadratic function, it can be queried as a quadratic
function and the result is quadratic.

However, it is not convertible to a variable:

```jldoctest objective_function; filter = r"MathOptInterface\\."s
julia> objective_function(model, VariableRef)
ERROR: InexactError: convert(MathOptInterface.VariableIndex, 1.0 + 2.0 MOI.VariableIndex(1))
[...]
```
"""
function objective_function(
    model::GenericModel,
    ::Type{F},
) where {F<:MOI.AbstractFunction}
    func = MOI.get(backend(model), MOI.ObjectiveFunction{F}())::F
    return jump_function(model, func)
end

function objective_function(model::GenericModel, ::Type{T}) where {T}
    return objective_function(model, moi_function_type(T))
end

function objective_function(model::GenericModel)
    F = MOI.get(backend(model), MOI.ObjectiveFunctionType())
    return objective_function(model, F)
end

"""
    set_objective_coefficient(
        model::GenericModel,
        variable::GenericVariableRef,
        coefficient::Real,
    )

Set the linear objective coefficient associated with `variable` to `coefficient`.

Note: this function will throw an error if a nonlinear objective is set.

## Example

```jldoctest
julia> model = Model();

julia> @variable(model, x);

julia> @objective(model, Min, 2x + 1)
2 x + 1

julia> set_objective_coefficient(model, x, 3)

julia> objective_function(model)
3 x + 1
```
"""
function set_objective_coefficient(
    model::GenericModel{T},
    variable::GenericVariableRef{T},
    coeff::Real,
) where {T}
    if _nlp_objective_function(model) !== nothing
        error(
            "A nonlinear objective created by the legacy `@NLobjective` is " *
            "set in the model. This does not support modification.",
        )
    end
    coeff_t = convert(T, coeff)::T
    F = objective_function_type(model)
    _set_objective_coefficient(model, variable, coeff_t, F)
    model.is_model_dirty = true
    return
end

function _set_objective_coefficient(
    model::GenericModel{T},
    variable::GenericVariableRef{T},
    coeff::T,
    ::Type{GenericVariableRef{T}},
) where {T}
    current_obj = objective_function(model)
    if index(current_obj) == index(variable)
        set_objective_function(model, coeff * variable)
    else
        set_objective_function(
            model,
            add_to_expression!(coeff * variable, current_obj),
        )
    end
    return
end

function _set_objective_coefficient(
    model::GenericModel{T},
    variable::GenericVariableRef{T},
    coeff::T,
    ::Type{F},
) where {T,F}
    MOI.modify(
        backend(model),
        MOI.ObjectiveFunction{moi_function_type(F)}(),
        MOI.ScalarCoefficientChange(index(variable), coeff),
    )
    return
end

"""
    set_objective_coefficient(
        model::GenericModel,
        variables::Vector{<:GenericVariableRef},
        coefficients::Vector{<:Real},
    )

Set multiple linear objective coefficients associated with `variables` to
`coefficients`, in a single call.

Note: this function will throw an error if a nonlinear objective is set.

## Example

```jldoctest
julia> model = Model();

julia> @variable(model, x);

julia> @variable(model, y);

julia> @objective(model, Min, 3x + 2y + 1)
3 x + 2 y + 1

julia> set_objective_coefficient(model, [x, y], [5, 4])

julia> objective_function(model)
5 x + 4 y + 1
```
"""
function set_objective_coefficient(
    model::GenericModel{T},
    variables::AbstractVector{<:GenericVariableRef{T}},
    coeffs::AbstractVector{<:Real},
) where {T}
    if _nlp_objective_function(model) !== nothing
        error(
            "A nonlinear objective created by the legacy `@NLobjective` is " *
            "set in the model. This does not support modification.",
        )
    end
    n, m = length(variables), length(coeffs)
    if !(n == m)
        msg = "The number of variables ($n) and coefficients ($m) must match"
        throw(DimensionMismatch(msg))
    end
    F = objective_function_type(model)
    _set_objective_coefficient(model, variables, convert.(T, coeffs), F)
    model.is_model_dirty = true
    return
end

function _set_objective_coefficient(
    model::GenericModel{T},
    variables::AbstractVector{<:GenericVariableRef{T}},
    coeffs::AbstractVector{<:T},
    ::Type{GenericVariableRef{T}},
) where {T}
    new_objective = LinearAlgebra.dot(coeffs, variables)
    current_obj = objective_function(model)::GenericVariableRef{T}
    if !(current_obj in variables)
        add_to_expression!(new_objective, current_obj)
    end
    set_objective_function(model, new_objective)
    return
end

function _set_objective_coefficient(
    model::GenericModel{T},
    variables::AbstractVector{<:GenericVariableRef{T}},
    coeffs::AbstractVector{<:T},
    ::Type{F},
) where {T,F}
    MOI.modify(
        backend(model),
        MOI.ObjectiveFunction{moi_function_type(F)}(),
        MOI.ScalarCoefficientChange.(index.(variables), coeffs),
    )
    return
end

"""
    set_objective_coefficient(
        model::GenericModel{T},
        variable_1::GenericVariableRef{T},
        variable_2::GenericVariableRef{T},
        coefficient::Real,
    ) where {T}

Set the quadratic objective coefficient associated with `variable_1` and
`variable_2` to `coefficient`.

Note: this function will throw an error if a nonlinear objective is set.

## Example

```jldoctest
julia> model = Model();

julia> @variable(model, x[1:2]);

julia> @objective(model, Min, x[1]^2 + x[1] * x[2])
x[1]² + x[1]*x[2]

julia> set_objective_coefficient(model, x[1], x[1], 2)

julia> set_objective_coefficient(model, x[1], x[2], 3)

julia> objective_function(model)
2 x[1]² + 3 x[1]*x[2]
```
"""
function set_objective_coefficient(
    model::GenericModel{T},
    variable_1::GenericVariableRef{T},
    variable_2::GenericVariableRef{T},
    coeff::Real,
) where {T}
    if _nlp_objective_function(model) !== nothing
        error(
            "A nonlinear objective created by the legacy `@NLobjective` is " *
            "set in the model. This does not support modification.",
        )
    end
    coeff_t = convert(T, coeff)::T
    F = moi_function_type(objective_function_type(model))
    _set_objective_coefficient(model, variable_1, variable_2, coeff_t, F)
    model.is_model_dirty = true
    return
end

function _set_objective_coefficient(
    model::GenericModel{T},
    variable_1::GenericVariableRef{T},
    variable_2::GenericVariableRef{T},
    coeff::T,
    ::Type{F},
) where {T,F}
    current_obj = objective_function(model)
    new_obj = add_to_expression!(coeff * variable_1 * variable_2, current_obj)
    set_objective_function(model, new_obj)
    return
end

function _set_objective_coefficient(
    model::GenericModel{T},
    variable_1::GenericVariableRef{T},
    variable_2::GenericVariableRef{T},
    coeff::T,
    ::Type{MOI.ScalarQuadraticFunction{T}},
) where {T}
    if variable_1 == variable_2
        coeff *= T(2)
    end
    MOI.modify(
        backend(model),
        MOI.ObjectiveFunction{MOI.ScalarQuadraticFunction{T}}(),
        MOI.ScalarQuadraticCoefficientChange(
            index(variable_1),
            index(variable_2),
            coeff,
        ),
    )
    return
end

"""
    set_objective_coefficient(
        model::GenericModel{T},
        variables_1::AbstractVector{<:GenericVariableRef{T}},
        variables_2::AbstractVector{<:GenericVariableRef{T}},
        coefficients::AbstractVector{<:Real},
    ) where {T}

Set multiple quadratic objective coefficients associated with `variables_1` and
`variables_2` to `coefficients`, in a single call.

Note: this function will throw an error if a nonlinear objective is set.

## Example

```jldoctest
julia> model = Model();

julia> @variable(model, x[1:2]);

julia> @objective(model, Min, x[1]^2 + x[1] * x[2])
x[1]² + x[1]*x[2]

julia> set_objective_coefficient(model, [x[1], x[1]], [x[1], x[2]], [2, 3])

julia> objective_function(model)
2 x[1]² + 3 x[1]*x[2]
```
"""
function set_objective_coefficient(
    model::GenericModel{T},
    variables_1::AbstractVector{<:GenericVariableRef{T}},
    variables_2::AbstractVector{<:GenericVariableRef{T}},
    coeffs::AbstractVector{<:Real},
) where {T}
    if _nlp_objective_function(model) !== nothing
        error(
            "A nonlinear objective created by the legacy `@NLobjective` is " *
            "set in the model. This does not support modification.",
        )
    end
    n1, n2, m = length(variables_1), length(variables_2), length(coeffs)
    if !(n1 == n2 == m)
        msg = "The number of variables ($n1, $n2) and coefficients ($m) must match"
        throw(DimensionMismatch(msg))
    end
    coeffs_t = convert.(T, coeffs)
    F = moi_function_type(objective_function_type(model))
    _set_objective_coefficient(model, variables_1, variables_2, coeffs_t, F)
    model.is_model_dirty = true
    return
end

function _set_objective_coefficient(
    model::GenericModel{T},
    variables_1::AbstractVector{<:V},
    variables_2::AbstractVector{<:V},
    coeffs::AbstractVector{<:T},
    ::Type{F},
) where {T,F,V<:GenericVariableRef{T}}
    new_obj = GenericQuadExpr{T,V}()
    add_to_expression!(new_obj, objective_function(model))
    for (c, x, y) in zip(coeffs, variables_1, variables_2)
        add_to_expression!(new_obj, c, x, y)
    end
    set_objective_function(model, new_obj)
    return
end

function _set_objective_coefficient(
    model::GenericModel{T},
    variables_1::AbstractVector{<:GenericVariableRef{T}},
    variables_2::AbstractVector{<:GenericVariableRef{T}},
    coeffs::AbstractVector{<:T},
    ::Type{MOI.ScalarQuadraticFunction{T}},
) where {T}
    for (i, x, y) in zip(eachindex(coeffs), variables_1, variables_2)
        if x == y
            coeffs[i] *= T(2)
        end
    end
    MOI.modify(
        backend(model),
        MOI.ObjectiveFunction{MOI.ScalarQuadraticFunction{T}}(),
        MOI.ScalarQuadraticCoefficientChange.(
            index.(variables_1),
            index.(variables_2),
            coeffs,
        ),
    )
    return
end

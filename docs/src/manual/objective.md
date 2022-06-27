```@meta
CurrentModule = JuMP
DocTestSetup = quote
    using JuMP
end
DocTestFilters = [r"≤|<=", r"≥|>=", r" == | = ", r" ∈ | in ", r"MathOptInterface|MOI"]
```

# Objectives

This page describes macros and functions related to linear and quadratic
objective functions only, unless otherwise indicated. For nonlinear objective
functions, see [Nonlinear Modeling](@ref).

## Set a linear objective

Use the [`@objective`](@ref) macro to set a linear objective function.

Use `Min` to create a minimization objective:
```jldoctest; setup = :(model=Model(); @variable(model, x))
julia> @objective(model, Min, 2x + 1)
2 x + 1
```

Use `Max` to create a maximization objective:
```jldoctest; setup = :(model=Model(); @variable(model, x))
julia> @objective(model, Max, 2x + 1)
2 x + 1
```

## Set a quadratic objective

Use the [`@objective`](@ref) macro to set a quadratic objective function.

Use `^2` to have a variable squared:
```jldoctest; setup = :(model=Model(); @variable(model, x))
julia> @objective(model, Min, x^2 + 2x + 1)
x² + 2 x + 1
```

You can also have bilinear terms between variables:
```jldoctest; setup = :(model=Model())
julia> @variable(model, x)
x

julia> @variable(model, y)
y

julia> @objective(model, Max, x * y + x + y)
x*y + x + y
```

## Query the objective function

Use [`objective_function`](@ref) to return the current objective function.
```jldoctest; setup = :(model=Model(); @variable(model, x))
julia> @objective(model, Min, 2x + 1)
2 x + 1

julia> objective_function(model)
2 x + 1
```

## Evaluate the objective function at a point

Use [`value`](@ref) to evaluate an objective function at a point specifying values for variables.

```jldoctest; setup = :(model=Model())
julia> @variable(model, x[1:2]);

julia> @objective(model, Min, 2x[1]^2 + x[1] + 0.5*x[2])
2 x[1]² + x[1] + 0.5 x[2]

julia> f = objective_function(model)
2 x[1]² + x[1] + 0.5 x[2]

julia> point = Dict(x[1] => 2.0, x[2] => 1.0);

julia> value(z -> point[z], f)
10.5
```

## Query the objective sense

Use [`objective_sense`](@ref) to return the current objective sense.
```jldoctest; setup = :(model=Model(); @variable(model, x))
julia> @objective(model, Min, 2x + 1)
2 x + 1

julia> objective_sense(model)
MIN_SENSE::OptimizationSense = 0
```

## Modify an objective

To modify an objective, call [`@objective`](@ref) with the new objective
function.
```jldoctest; setup = :(model=Model(); @variable(model, x))
julia> @objective(model, Min, 2x)
2 x

julia> @objective(model, Max, -2x)
-2 x
```

Alternatively, use [`set_objective_function`](@ref).

```jldoctest; setup = :(model=Model(); @variable(model, x))
julia> @objective(model, Min, 2x)
2 x

julia> new_objective = @expression(model, -2 * x)
-2 x

julia> set_objective_function(model, new_objective)
```

## Modify an objective coefficient

Use [`set_objective_coefficient`](@ref) to modify an objective coefficient.
```jldoctest; setup = :(model=Model(); @variable(model, x))
julia> @objective(model, Min, 2x)
2 x

julia> set_objective_coefficient(model, x, 3)

julia> objective_function(model)
3 x
```

!!! info
    There is no way to modify the coefficient of a quadratic term. Set a new
    objective instead.

## Modify the objective sense

Use [`set_objective_sense`](@ref) to modify the objective sense.
```jldoctest; setup = :(model=Model(); @variable(model, x))
julia> @objective(model, Min, 2x)
2 x

julia> objective_sense(model)
MIN_SENSE::OptimizationSense = 0

julia> set_objective_sense(model, MAX_SENSE);

julia> objective_sense(model)
MAX_SENSE::OptimizationSense = 1
```

Alternatively, call [`@objective`](@ref) and pass the existing objective
function.
```jldoctest; setup = :(model=Model(); @variable(model, x))
julia> @objective(model, Min, 2x)
2 x

julia> @objective(model, Max, objective_function(model))
2 x
```

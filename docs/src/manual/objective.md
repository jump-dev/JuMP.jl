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
```jldoctest
julia> model = Model();

julia> @variable(model, x);

julia> @objective(model, Min, 2x + 1)
2 x + 1
```

Use `Max` to create a maximization objective:
```jldoctest
julia> model = Model();

julia> @variable(model, x);

julia> @objective(model, Max, 2x + 1)
2 x + 1
```

## Set a quadratic objective

Use the [`@objective`](@ref) macro to set a quadratic objective function.

Use `^2` to have a variable squared:
```jldoctest
julia> model = Model();

julia> @variable(model, x);

julia> @objective(model, Min, x^2 + 2x + 1)
x² + 2 x + 1
```

You can also have bilinear terms between variables:
```jldoctest
julia> model = Model();

julia> @variable(model, x)
x

julia> @variable(model, y)
y

julia> @objective(model, Max, x * y + x + y)
x*y + x + y
```

## Set a nonlinear objective

Use the [`@objective`](@ref) macro to set a nonlinear objective function:

```jldoctest
julia> model = Model();

julia> @variable(model, x <= 1);

julia> @objective(model, Max, log(x))
log(x)
```

## Query the objective function

Use [`objective_function`](@ref) to return the current objective function.
```jldoctest
julia> model = Model();

julia> @variable(model, x);

julia> @objective(model, Min, 2x + 1)
2 x + 1

julia> objective_function(model)
2 x + 1
```

## Evaluate the objective function at a point

Use [`value`](@ref) to evaluate an objective function at a point specifying values for variables.

```jldoctest
julia> model = Model();

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
```jldoctest
julia> model = Model();

julia> @variable(model, x);

julia> @objective(model, Min, 2x + 1)
2 x + 1

julia> objective_sense(model)
MIN_SENSE::OptimizationSense = 0
```

## Modify an objective

To modify an objective, call [`@objective`](@ref) with the new objective
function.
```jldoctest
julia> model = Model();

julia> @variable(model, x);

julia> @objective(model, Min, 2x)
2 x

julia> @objective(model, Max, -2x)
-2 x
```

Alternatively, use [`set_objective_function`](@ref).

```jldoctest
julia> model = Model();

julia> @variable(model, x);

julia> @objective(model, Min, 2x)
2 x

julia> new_objective = @expression(model, -2 * x)
-2 x

julia> set_objective_function(model, new_objective)
```

## Modify an objective coefficient

Use [`set_objective_coefficient`](@ref) to modify an objective coefficient.
```jldoctest
julia> model = Model();

julia> @variable(model, x);

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
```jldoctest
julia> model = Model();

julia> @variable(model, x);

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
```jldoctest
julia> model = Model();

julia> @variable(model, x);

julia> @objective(model, Min, 2x)
2 x

julia> @objective(model, Max, objective_function(model))
2 x
```

## Set a vector-valued objective

Define a multi-objective optimization problem by passing a vector of objectives:

```jldoctest
julia> model = Model();

julia> @variable(model, x[1:2]);

julia> @objective(model, Min, [1 + x[1], 2 * x[2]])
2-element Vector{AffExpr}:
 x[1] + 1
 2 x[2]

julia> f = objective_function(model)
2-element Vector{AffExpr}:
 x[1] + 1
 2 x[2]
```

!!! tip
    The [Multi-objective knapsack](@ref) tutorial provides an example of
    solving a multi-objective integer program.

In most cases, multi-objective optimization solvers will return multiple
solutions, corresponding to points on the Pareto frontier. See [Multiple solutions](@ref)
for information on how to query and work with multiple solutions.

Note that you must set a single objective sense, that is, you cannot have
both minimization and maximization objectives. Work around this limitation by
choosing `Min` and negating any objectives you want to maximize:

```jldoctest
julia> model = Model();

julia> @variable(model, x[1:2]);

julia> @expression(model, obj1, 1 + x[1])
x[1] + 1

julia> @expression(model, obj2, 2 * x[1])
2 x[1]

julia> @objective(model, Min, [obj1, -obj2])
2-element Vector{AffExpr}:
 x[1] + 1
 -2 x[1]
```

Defining your objectives as expressions allows flexibility in how you can solve
variations of the same problem, with some objectives removed and constrained to
be no worse that a fixed value.

```jldoctest
julia> model = Model();

julia> @variable(model, x[1:2]);

julia> @expression(model, obj1, 1 + x[1])
x[1] + 1

julia> @expression(model, obj2, 2 * x[1])
2 x[1]

julia> @expression(model, obj3, x[1] + x[2])
x[1] + x[2]

julia> @objective(model, Min, [obj1, obj2, obj3])  # Three-objective problem
3-element Vector{AffExpr}:
 x[1] + 1
 2 x[1]
 x[1] + x[2]

julia> # optimize!(model), look at the solution, talk to stakeholders, then
       # decide you want to solve a new problem where the third objective is
       # removed and constrained to be better than 2.0.
       nothing

julia> @objective(model, Min, [obj1, obj2])   # Two-objective problem
2-element Vector{AffExpr}:
 x[1] + 1
 2 x[1]

julia> @constraint(model, obj3 <= 2.0)
x[1] + x[2] ≤ 2
```

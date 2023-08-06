```@meta
CurrentModule = JuMP
DocTestSetup = quote
    using JuMP
end
DocTestFilters = [r"≤|<=", r"≥|>=", r" == | = ", r" ∈ | in ", r"MathOptInterface|MOI"]
```

# [Nonlinear Modeling](@id new_nonlinear_interface)

!!! warning
    This page describes an experimental nonlinear interface to JuMP. The API
    described below is stable, and it will not break with future 1.X releases of
    JuMP. However, solver support may be limited, and there may be gaps in
    functionality compared with [Nonlinear Modeling](@ref). To report a bug, or
    request a missing feature, please [open an issue](https://github.com/jump-dev/JuMP.jl/issues/new/choose).

JuMP has support for general smooth nonlinear (convex and nonconvex)
optimization problems. JuMP is able to provide exact, sparse second-order
derivatives to solvers. This information can improve solver accuracy and
performance.

## Set a nonlinear objective

Use [`@objective`](@ref) to set a nonlinear objective.

```jldoctest
julia> model = Model();

julia> @variable(model, x[1:2]);

julia> @objective(model, Min, exp(x[1]) - sqrt(x[2]))
exp(x[1]) - sqrt(x[2])
```

To modify a nonlinear objective, call [`@objective`](@ref) again.

## Add a nonlinear constraint

Use [`@constraint`](@ref) to add a nonlinear constraint.

```jldoctest nonlinear_constraint
julia> model = Model();

julia> @variable(model, x[1:2]);

julia> @constraint(model, exp(x[1]) <= 1)
exp(x[1]) - 1.0 ≤ 0

julia> @constraint(model, con[i = 1:2], 2^x[i] >= i)
2-element Vector{ConstraintRef{Model, MathOptInterface.ConstraintIndex{MathOptInterface.ScalarNonlinearFunction, MathOptInterface.GreaterThan{Float64}}, ScalarShape}}:
 con[1] : (2.0 ^ x[1]) - 1.0 ≥ 0
 con[2] : (2.0 ^ x[2]) - 2.0 ≥ 0
```

Delete a nonlinear constraint using [`delete`](@ref):
```jldoctest nonlinear_constraint
julia> delete(model, con[1])
```

## Create a nonlinear expression

Use [`@expression`](@ref) to create nonlinear expression objects. The syntax
is identical to [`@expression`](@ref), except that the expression can contain
nonlinear terms.

```jldoctest nl_expression
julia> model = Model();

julia> @variable(model, x[1:2]);

julia> expr = @expression(model, exp(x[1]) + sqrt(x[2]))
exp(x[1]) + sqrt(x[2])

julia> my_anon_expr = @expression(model, [i = 1:2], sin(x[i]))
2-element Vector{NonlinearExpr}:
 sin(x[1])
 sin(x[2])

julia> @expression(model, my_expr[i = 1:2], sin(x[i]))
2-element Vector{NonlinearExpr}:
 sin(x[1])
 sin(x[2])
```

A [`NonlinearExpr`](@ref) can be used in [`@objective`](@ref),
[`@constraint`](@ref), and even nested in other [`@expression`](@ref)s.

```jldoctest nl_expression
julia> @objective(model, Min, expr^2 + 1)
((exp(x[1]) + sqrt(x[2])) ^ 2.0) + 1.0

julia> @constraint(model, [i = 1:2], my_expr[i] <= i)
2-element Vector{ConstraintRef{Model, MathOptInterface.ConstraintIndex{MathOptInterface.ScalarNonlinearFunction, MathOptInterface.LessThan{Float64}}, ScalarShape}}:
 sin(x[1]) - 1.0 ≤ 0
 sin(x[2]) - 2.0 ≤ 0

julia> @expression(model, nested[i = 1:2], sin(my_expr[i]))
2-element Vector{NonlinearExpr}:
 sin(sin(x[1]))
 sin(sin(x[2]))
```

Use [`value`](@ref) to query the value of a nonlinear expression:

```jldoctest nl_expression
julia> set_start_value(x[1], 1.0)

julia> value(start_value, nested[1])
0.7456241416655579

julia> sin(sin(1.0))
0.7456241416655579
```

## User-defined functions

In addition to a standard list of univariate and multivariate functions
recognized by the `MOI.Nonlinear` submodule, JuMP supports *user-defined*
Julia functions.

!!! warning
    User-defined functions must return a scalar output. For a work-around, see
    [User-defined functions with vector outputs](@ref).

### Register a function

Register a user-defined function using the [`@register`](@ref) macro:

```@repl
using JuMP
square(x) = x^2
f(x, y) = (x - 1)^2 + (y - 2)^2
model = Model();
@register(model, my_square, 1, square)
@register(model, my_f, 2, f)
@variable(model, x[1:2]);
@objective(model, Min, my_f(x[1], my_square(x[2])))
```

The arguments to [`@register`](@ref) are:

 1. The model in which the function is registered.
 2. A Julia symbol object which serves as the name of the user-defined function
    in JuMP expressions. This name must not be the same as that of the function.
 3. The number of scalar input arguments that the function takes.
 4. A Julia method which computes the function.

!!! warning
    User-defined functions cannot be re-registered and will not update if you
    modify the underlying Julia function. If you want to change a user-defined
    function between solves, rebuild the model or use a different name.

### Registered functions without macros

The [`@register`](@ref) macro is syntactic sugar for the
[`add_user_defined_function`](@ref) method. Thus, the non-macro version of the
preceding example is:

```@repl
using JuMP
square(x) = x^2
f(x, y) = (x - 1)^2 + (y - 2)^2
model = Model();
my_square = add_user_defined_function(model, :my_square, 1, square)
my_f = add_user_defined_function(model, :my_f, 2, f)
@variable(model, x[1:2]);
@objective(model, Min, my_f(x[1], my_square(x[2])))
```

This has two important consequences.

First, you cannot register a user-defined function with the same name as an
existing function. For example, a call to [`@register`](@ref) like:
```julia
julia> @register(model, square, 1, square)
```
will error because it is equivalent to:
```julia
julia> square = add_user_defined_function(model, :square, 1, square)
ERROR: invalid redefinition of constant square
Stacktrace:
[...]
```
and `square` already exists as a Julia function.

Second, you can construct and use [`UserDefinedFunction`](@ref)s outside the
macros.

```@repl
using JuMP
square(x) = x^2
model = Model();
@register(model, my_square, 1, square)
@variable(model, x)
typeof(my_square)
x_squared = my_square(x)
typeof(x_squared)
my_square_2 = UserDefinedFunction(:my_square)
my_square_2(x_squared)
```

### Register gradients and Hessians

By default, JuMP will use automatic differentiation to compute the gradient and
Hessian of user-defined functions. If your function is not amenable to
automatic differentiation, or you can compute analytic derivatives, you may pass
additional arguments to [`@register`](@ref) to compute the first- and
second-derivatives.

#### Univariate functions

For univariate functions, a gradient function `∇f` returns a number that
represents the first-order derivative. You may, in addition, pass a third
function which returns a number representing the second-order derivative:
```@repl
using JuMP
f(x) = x^2
∇f(x) = 2x
∇²f(x) = 2
model = Model();
@register(model, my_square, 1, f, ∇f, ∇²f)  # Providing ∇²f is optional
@variable(model, x)
@objective(model, Min, my_square(x))
```

#### Multivariate functions

For multivariate functions, the gradient function `∇f` must take an
`AbstractVector` as the first argument that is filled in-place. The Hessian
function, `∇²f`, must take an `AbstractMatrix` as the first argument, the
lower-triangular of which is filled in-place:
```@repl
using JuMP
f(x...) = (1 - x[1])^2 + 100 * (x[2] - x[1]^2)^2
function ∇f(g::AbstractVector{T}, x::T...) where {T}
    g[1] = 400 * x[1]^3 - 400 * x[1] * x[2] + 2 * x[1] - 2
    g[2] = 200 * (x[2] - x[1]^2)
    return
end
function ∇²f(H::AbstractMatrix{T}, x::T...) where {T}
    H[1, 1] = 1200 * x[1]^2 - 400 * x[2] + 2
    # H[1, 2] = -400 * x[1]  <-- Not needed. Fill the lower-triangular only.
    H[2, 1] = -400 * x[1]
    H[2, 2] = 200.0
    return
end
model = Model();
@register(model, rosenbrock, 2, f, ∇f, ∇²f)  # Providing ∇²f is optional
@variable(model, x[1:2])
@objective(model, Min, rosenbrock(x[1], x[2]))
```

You may assume the Hessian matrix `H` is initialized with zeros, and because `H`
is symmetric, you need only to fill in the non-zero of the lower-triangular
terms. The matrix type passed in as `H` depends on the automatic differentiation
system, so make sure the first argument to the Hessian function supports an
`AbstractMatrix` (it may be something other than `Matrix{Float64}`). Moreover,
you may assume only that `H` supports `size(H)` and `setindex!`. Finally, the
matrix is treated as dense, so the performance will be poor on functions with
high-dimensional input.

### User-defined functions with vector inputs

User-defined functions which take vectors as input arguments (for example,
`f(x::Vector)`) are *not* supported. Instead, use Julia's splatting syntax to
create a function with scalar arguments. For example, instead of:
```julia
f(x::Vector) = sum(x[i]^i for i in 1:length(x))
```
define:
```julia
f(x...) = sum(x[i]^i for i in 1:length(x))
```

Another approach is to define the splatted function as an anonymous function:
```@repl
using JuMP
model = Model();
@variable(model, x[1:5])
f(x::Vector) = sum(x[i]^i for i in 1:length(x))
@register(model, my_f, 5, (x...) -> f(collect(x)))
@objective(model, Min, my_f(x...))
```

### Automatic differentiation

JuMP does not support black-box optimization, so all user-defined functions must
provide derivatives in some form. Fortunately, JuMP supports automatic
differentiation of user-defined functions.

!!! info
    Automatic differentiation is *not* finite differencing. JuMP's automatically
    computed derivatives are not subject to approximation error.

JuMP uses [ForwardDiff.jl](https://github.com/JuliaDiff/ForwardDiff.jl) to
perform automatic differentiation; see the ForwardDiff.jl
[documentation](https://www.juliadiff.org/ForwardDiff.jl/v0.10.2/user/limitations.html)
for a description of how to write a function suitable for automatic
differentiation.

#### Common mistakes when writing a user-defined function

!!! warning
    Get an error like `No method matching Float64(::ForwardDiff.Dual)`? Read
    this section, and see the guidelines at [ForwardDiff.jl](https://www.juliadiff.org/ForwardDiff.jl/release-0.10/user/limitations.html).

The most common error is that your user-defined function is not generic with
respect to the number type, that is, don't assume that the input to the function
is `Float64`.
```julia
f(x::Float64) = 2 * x  # This will not work.
f(x::Real)    = 2 * x  # This is good.
f(x)          = 2 * x  # This is also good.
```

Another reason you may encounter this error is if you create arrays inside
your function which are `Float64`.
```julia
function bad_f(x...)
    y = zeros(length(x))  # This constructs an array of `Float64`!
    for i = 1:length(x)
        y[i] = x[i]^i
    end
    return sum(y)
end

function good_f(x::T...) where {T<:Real}
    y = zeros(T, length(x))  # Construct an array of type `T` instead!
    for i = 1:length(x)
        y[i] = x[i]^i
    end
    return sum(y)
end
```

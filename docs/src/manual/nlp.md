```@meta
CurrentModule = JuMP
DocTestSetup = quote
    using JuMP
end
DocTestFilters = [r"≤|<=", r"≥|>=", r" == | = ", r" ∈ | in ", r"MathOptInterface|MOI"]
```

# Nonlinear Modeling

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
-(exp(x[1]), sqrt(x[2]))
```

To modify a nonlinear objective, call [`@objective`](@ref) again.

## Add a nonlinear constraint

Use [`@constraint`](@ref) to add a nonlinear constraint.

```jldoctest nonlinear_constraint
julia> model = Model();

julia> @variable(model, x[1:2]);

julia> @constraint(model, exp(x[1]) <= 1)
-(exp(x[1]), 1.0) ≤ 0.0

julia> @constraint(model, con[i = 1:2], 2^x[i] >= i)
2-element Vector{ConstraintRef{Model, MathOptInterface.ConstraintIndex{MathOptInterface.ScalarNonlinearFunction, MathOptInterface.GreaterThan{Float64}}, ScalarShape}}:
 con[1] : -(^(2.0, x[1]), 1.0) ≥ 0.0
 con[2] : -(^(2.0, x[2]), 2.0) ≥ 0.0
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
+(exp(x[1]), sqrt(x[2]))

julia> my_anon_expr = @expression(model, [i = 1:2], sin(x[i]))
2-element Vector{NonlinearExpr}:
 sin(x[1])
 sin(x[2])

julia> @expression(model, my_expr[i = 1:2], sin(x[i]))
2-element Vector{NonlinearExpr}:
 sin(x[1])
 sin(x[2])
```

Nonlinear expression can be used in [`@objective`](@ref), [`@constraint`](@ref),
and even nested in other [`@expression`](@ref)s.

```jldoctest nl_expression
julia> @objective(model, Min, expr^2 + 1)
+(^(+(exp(x[1]), sqrt(x[2])), 2.0), 1.0)

julia> @constraint(model, [i = 1:2], my_expr[i] <= i)
2-element Vector{ConstraintRef{Model, MathOptInterface.ConstraintIndex{MathOptInterface.ScalarNonlinearFunction, MathOptInterface.LessThan{Float64}}, ScalarShape}}:
 -(sin(x[1]), 1.0) ≤ 0.0
 -(sin(x[2]), 2.0) ≤ 0.0

julia> @expression(model, nested[i = 1:2], sin(my_expr[i]))
2-element Vector{NonlinearExpr}:
 sin(sin(x[1]))
 sin(sin(x[2]))
```

Use [`value`](@ref) to query the value of a nonlinear expression:
```jldoctest nl_expression
julia> set_start_value(x[1], 1.0)

julia> 0.7456241416655579 # value(start_value, nested[1])  # TODO
0.7456241416655579

julia> sin(sin(1.0))
0.7456241416655579
```

## User-defined Functions

JuMP natively supports the set of univariate and multivariate functions
recognized by the `MOI.Nonlinear` submodule. In addition to this list of
functions, it is possible to register custom *user-defined* nonlinear functions.
User-defined functions can be used anywhere in [`@objective`](@ref),
[`@constraint`](@ref), and [`@expression`](@ref).

!!! warning
    User-defined functions must return a scalar output. For a work-around, see
    [User-defined functions with vector outputs](@ref).

### Automatic differentiation

JuMP does not support black-box optimization, so all user-defined functions must
provide derivatives in some form. Fortunately, JuMP supports **automatic
differentiation of user-defined functions**, a feature to our knowledge not
available in any comparable modeling systems.

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

### Register a function

To register a user-defined function with derivatives computed by
automatic differentiation, use the [`@register`](@ref) macro as in the following
example:

```@example
using JuMP #hide
square(x) = x^2
f(x, y) = (x - 1)^2 + (y - 2)^2
model = Model()
@register(model, my_square, 1, square)
@register(model, my_f, 2, f)
@variable(model, x[1:2] >= 0.5)
@objective(model, Min, my_f(x[1], my_square(x[2])))
```

The above code creates a JuMP model with the objective function
`(x[1] - 1)^2 + (x[2]^2 - 2)^2`. The arguments to [`@register`](@ref) are:
 1. The model for which the functions are registered.
 2. A Julia symbol object which serves as the name of the user-defined function
    in JuMP expressions.
 3. The number of input arguments that the function takes.
 4. The Julia method which computes the function

!!! warning
    User-defined functions cannot be re-registered and will not update if you
    modify the underlying Julia function. If you want to change a user-defined
    function between solves, rebuild the model or use a different name. To use
    a different name programmatically, see [Raw expression input](@ref).

### Register a function and gradient

Forward-mode automatic differentiation as implemented by ForwardDiff.jl has a
computational cost that scales linearly with the number of input dimensions. As
such, it is not the most efficient way to compute gradients of user-defined
functions if the number of input arguments is large. In this case, users may
want to provide their own routines for evaluating gradients.

#### Univariate functions

For univariate functions, the gradient function `∇f` returns a number that
represents the first-order derivative:
```@example
using JuMP #hide
f(x) = x^2
∇f(x) = 2x
model = Model()
@register(model, my_square, 1, f, ∇f)
@variable(model, x >= 0)
@objective(model, Min, my_square(x))
```

#### Multivariate functions

For multivariate functions, the gradient function `∇f` must take a gradient
vector as the first argument that is filled in-place:
```@example
using JuMP #hide
f(x, y) = (x - 1)^2 + (y - 2)^2
function ∇f(g::AbstractVector{T}, x::T, y::T) where {T}
    g[1] = 2 * (x - 1)
    g[2] = 2 * (y - 2)
    return
end
model = Model()
@register(model, my_square, 2, f, ∇f)
@variable(model, x[1:2] >= 0)
@objective(model, Min, my_square(x[1], x[2]))
```

!!! warning
    Make sure the first argument to `∇f` supports an `AbstractVector`, and do
    not assume the input is `Float64`.

### Register a function, gradient, and hessian

You can also register a function with the second-order derivative information,
which is a scalar for univariate functions, and a symmetric matrix for
multivariate functions.

#### Univariate functions

Pass a function which returns a number representing the second-order derivative:
```@example
using JuMP #hide
f(x) = x^2
∇f(x) = 2x
∇²f(x) = 2
model = Model()
@register(model, my_square, 1, f, ∇f, ∇²f)
@variable(model, x >= 0)
@objective(model, Min, my_square(x))
```

#### Multivariate functions

For multivariate functions, the hessian function `∇²f` must take an
`AbstractMatrix` as the first argument, the lower-triangular of which is filled
in-place:
```@example
using JuMP #hide
f(x...) = (1 - x[1])^2 + 100 * (x[2] - x[1]^2)^2
function ∇f(g, x...)
    g[1] = 400 * x[1]^3 - 400 * x[1] * x[2] + 2 * x[1] - 2
    g[2] = 200 * (x[2] - x[1]^2)
    return
end
function ∇²f(H, x...)
    H[1, 1] = 1200 * x[1]^2 - 400 * x[2] + 2
    # H[1, 2] = -400 * x[1]  <-- Not needed. Fill the lower-triangular only.
    H[2, 1] = -400 * x[1]
    H[2, 2] = 200.0
    return
end
model = Model()
@register(model, rosenbrock, 2, f, ∇f, ∇²f)
@variable(model, x[1:2])
@objective(model, Min, rosenbrock(x[1], x[2]))
```

!!! warning
    You may assume the Hessian matrix `H` is initialized with zeros, and because
    `H` is symmetric, you need only to fill in the non-zero of the
    lower-triangular terms. The matrix type passed in as `H` depends on the
    automatic differentiation system, so make sure the first argument to the
    Hessian function supports an `AbstractMatrix` (it may be something other
    than `Matrix{Float64}`). However, you may assume only that `H` supports
    `size(H)` and `setindex!`. Finally, the matrix is treated as dense, so the
    performance will be poor on functions with high-dimensional input.

### User-defined functions with vector inputs

User-defined functions which take vectors as input arguments (for example,
`f(x::Vector)`) are *not* supported. Instead, use Julia's splatting syntax to
create a function with scalar arguments. For example, instead of
```julia
f(x::Vector) = sum(x[i]^i for i in 1:length(x))
```
define:
```julia
f(x...) = sum(x[i]^i for i in 1:length(x))
```

This function `f` can be used in a JuMP model as follows:
```@example
using JuMP #hide
model = Model()
@variable(model, x[1:5] >= 0)
f(x...) = sum(x[i]^i for i in 1:length(x))
@register(model, my_f, 5, f)
@objective(model, Min, my_f(x...))
```

## Factors affecting solution time

The execution time when solving a nonlinear programming problem can be divided
into two parts, the time spent in the optimization algorithm (the solver) and
the time spent evaluating the nonlinear functions and corresponding derivatives.
Ipopt explicitly displays these two timings in its output, for example:

```
Total CPU secs in IPOPT (w/o function evaluations)   =      7.412
Total CPU secs in NLP function evaluations           =      2.083
```

For Ipopt in particular, one can improve the performance by installing advanced
sparse linear algebra packages, see [Installation Guide](@ref). For other
solvers, see their respective documentation for performance tips.

The function evaluation time, on the other hand, is the responsibility of the
modeling language. JuMP computes derivatives by using reverse-mode automatic
differentiation with graph coloring methods for exploiting sparsity of the
Hessian matrix. As a conservative bound, JuMP's performance here currently
may be expected to be within a factor of 5 of AMPL's. Our [paper in
SIAM Review](https://mlubin.github.io/pdf/jump-sirev.pdf) has more details.

```@meta
CurrentModule = JuMP
DocTestSetup = quote
    using JuMP
end
DocTestFilters = [r"≤|<=", r"≥|>=", r" == | = ", r" ∈ | in ", r"MathOptInterface|MOI"]
```

# [Nonlinear Modeling](@id new_nonlinear_interface)

!!! warning
    This page describes a new nonlinear interface to JuMP. It replaces the
    legacy `@NL` interface, which is documented at [Nonlinear Modeling](@ref).
    The API described below is stable, and it will not break with future 1.X
    releases of JuMP. However, solver support may be limited, and there may be
    gaps in functionality compared with the legacy interface. To report a bug,
    or request a missing feature, please [open an issue](https://github.com/jump-dev/JuMP.jl/issues/new/choose).

JuMP has support for nonlinear (convex and nonconvex) optimization problems.
JuMP is able to automatically provide exact, sparse second-order derivatives to
solvers. This information can improve solver accuracy and performance.

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

Use [`@expression`](@ref) to create nonlinear expression objects:

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

## Nonlinear expressions in detail

Nonlinear expressions in JuMP are represented by a [`NonlinearExpr`](@ref)
object.

### Constructors

Nonlinear expressions can be created using the [`NonlinearExpr`](@ref)
constructors:

```jldoctest nonlinear_expressions
julia> model = Model();

julia> @variable(model, x);

julia> expr = NonlinearExpr(:sin, Any[x])
sin(x)
```

or via operator overloading:

```jldoctest
julia> model = Model();

julia> @variable(model, x);

julia> expr = sin(x)
sin(x)
```

### Supported arguments

Nonlinear expressions can contain a mix of numbers, [`AffExpr`](@ref),
[`QuadExpr`](@ref), and other [`NonlinearExpr`](@ref):

```jldoctest
julia> model = Model();

julia> @variable(model, x);

julia> aff = x + 1;

julia> quad = x^2 + x;

julia> expr = cos(x) * sin(quad) + aff
(cos(x) * sin(x² + x)) + (x + 1)
```

### Supported operators

The list of supported operators may vary between solvers. Given an optimizer,
query the list of supported operators using
[`MOI.ListOfSupportedNonlinearOperators`](@ref):
```jldoctest; filter=[r":.+", r"[0-9]+\-element"]
julia> import Ipopt

julia> import MathOptInterface as MOI

julia> MOI.get(Ipopt.Optimizer(), MOI.ListOfSupportedNonlinearOperators())
85-element Vector{Symbol}:
 :+
 :-
 :abs
 :sqrt
 :cbrt
 :abs2
 :inv
 :log
 :log10
 :log2
 ⋮
 :min
 :max
 :&&
 :||
 :<=
 :(==)
 :>=
 :<
 :>
```

In some univariate cases, the operator is defined in [`SpecialFunctions.jl`](https://github.com/JuliaMath/SpecialFunctions.jl).
To use these functions, you must explicitly import `SpecialFunctions.jl`
```jldoctest
julia> import Ipopt

julia> op = MOI.get(Ipopt.Optimizer(), MOI.ListOfSupportedNonlinearOperators());

julia> :erfcx in op
true

julia> :dawson in op
true

julia> import SpecialFunctions

julia> model = Model();

julia> @variable(model, x)
x

julia> @expression(model, SpecialFunctions.erfcx(x))
erfcx(x)

julia> @expression(model, SpecialFunctions.dawson(x))
dawson(x)
```

### Limitations

Some nonlinear expressions cannot be created via operator overloading. For
example, to minimize the likelihood of bugs in user-code, we have not overloaded
comparisons such as `<` and `>=` between JuMP objects:

```jldoctest
julia> model = Model();

julia> @variable(model, x);

julia> x < 1
ERROR: Cannot evaluate `<` between a variable and a number.
[...]
```

Instead, wrap the expression in the [`@expression`](@ref) macro:
```jldoctest
julia> model = Model();

julia> @variable(model, x);

julia> expr = @expression(model, x < 1)
x < 1
```

For technical reasons, other operators that are not overloaded include `||`,
`&&`, and `ifelse`.

```jldoctest
julia> model = Model();

julia> @variable(model, x);

julia> expr = @expression(model, ifelse(x < -1 || x >= 1, x^2, 0.0))
ifelse((x < -1) || (x >= 1), x², 0.0)
```

As an alternative, use the `JuMP.op_` functions, which fallback to the
various comparison and logical operators:
```jldoctest
julia> model = Model();

julia> @variable(model, x);

julia> expr = op_ifelse(
           op_or(op_less_than(x, -1), op_greater_equal(x, 1)),
           x^2,
           0.0,
       )
ifelse((x < -1) || (x >= 1), x², 0.0)
```

The available functions are:

| JuMP function                     | Julia function |
| :-------------------------------- | :------------- |
| [`op_ifelse`](@ref)        | `ifelse`       |
| [`op_and`](@ref)           | `&&`           |
| [`op_or`](@ref)            | `\|\|`           |
| [`op_greater_than`](@ref)  | `>`            |
| [`op_greater_equal`](@ref) | `>=`           |
| [`op_less_than`](@ref)     | `<`            |
| [`op_less_equal`](@ref)    | `<=`           |
| [`op_equal_to`](@ref)      | `==`           |

### Fields

Each [`NonlinearExpr`](@ref) has two fields.

The `.head` field is a `Symbol` that represents the operator being called:

```jldoctest nonlinear_expressions
julia> expr.head
:sin
```

The `.args` field is a `Vector{Any}` containing the arguments to the operator:

```jldoctest nonlinear_expressions
julia> expr.args
1-element Vector{Any}:
 x
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
@register(model, op_square, 1, square)
@register(model, op_f, 2, f)
@variable(model, x[1:2]);
@objective(model, Min, op_f(x[1], op_square(x[2])))
```

The arguments to [`@register`](@ref) are:

 1. The model in which the function is registered.
 2. A Julia symbol object which serves as the name of the user-defined function
    in JuMP expressions. This name must not be the same as that of the function.
 3. The number of scalar input arguments that the function takes.
 4. A Julia method which computes the function.

!!! warning
    User-defined functions cannot be re-registered or deleted.

You can obtain a reference to the operator using the `model[:key]` syntax:

```@repl
using JuMP
square(x) = x^2
model = Model();
@register(model, op_square, 1, square)
op_square_2 = model[:op_square]
```

### Registered functions without macros

The [`@register`](@ref) macro is syntactic sugar for the
[`register_nonlinear_operator`](@ref) method. Thus, the non-macro version of the
preceding example is:

```@repl
using JuMP
square(x) = x^2
f(x, y) = (x - 1)^2 + (y - 2)^2
model = Model();
op_square = register_nonlinear_operator(model, 1, square; name = :op_square)
model[:op_square] = op_square
op_f = register_nonlinear_operator(model, 2, f; name = :op_f)
model[:op_f] = op_f
@variable(model, x[1:2]);
@objective(model, Min, op_f(x[1], op_square(x[2])))
```

### Registering with the same name as an existing function

A common error encountered is the following:
```jldoctest nonlinear_invalid_redefinition
julia> using JuMP

julia> model = Model();

julia> f(x) = x^2
f (generic function with 1 method)

julia> @register(model, f, 1, f)
ERROR: Unable to register the nonlinear operator `:f` with the same name as
an existing function.
[...]
```
This error occurs because `@register(model, f, 1, f)` is equivalent to:
```julia
julia> f = register_nonlinear_operator(model, 1, f; name = :f)
```
but `f` already exists as a Julia function.

If you evaluate the function without registering it, JuMP will trace the
function using operator overloading:
```jldoctest nonlinear_invalid_redefinition
julia> @variable(model, x);

julia> f(x)
x²
```

To force JuMP to treat `f` as a user-defined function and not trace it, register
the function using [`register_nonlinear_operator`](@ref) and define a new method
which manually creates a [`NonlinearExpr`](@ref):
```jldoctest nonlinear_invalid_redefinition
julia> _ = register_nonlinear_operator(model, 1, f; name = :f)
NonlinearOperator(:f, f)

julia> f(x::AbstractJuMPScalar) = NonlinearExpr(:f, Any[x])
f (generic function with 2 methods)

julia> @expression(model, log(f(x)))
log(f(x))
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
@register(model, op_f, 1, f, ∇f, ∇²f)  # Providing ∇²f is optional
@variable(model, x)
@objective(model, Min, op_f(x))
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
is symmetric, you need only to fill in the non-zero lower-triangular
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
@register(model, op_f, 5, (x...) -> f(collect(x)))
@objective(model, Min, op_f(x...))
```

### Automatic differentiation

JuMP does not support black-box optimization, so all user-defined functions must
provide derivatives in some form. Fortunately, JuMP supports automatic
differentiation of user-defined functions.

!!! info
    Automatic differentiation is *not* finite differencing. JuMP's automatically
    computed derivatives are not subject to approximation error.

JuMP uses [ForwardDiff.jl](https://github.com/JuliaDiff/ForwardDiff.jl) to
perform automatic differentiation of user-defined functions; see the ForwardDiff.jl
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

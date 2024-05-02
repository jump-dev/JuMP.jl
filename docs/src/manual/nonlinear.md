```@meta
CurrentModule = JuMP
DocTestSetup = quote
    using JuMP
end
DocTestFilters = [r"≤|<=", r"≥|>=", r" == | = ", r" ∈ | in ", r"MathOptInterface|MOI"]
```

# Nonlinear Modeling

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

## Add a parameter

Some solvers have explicit support for parameters, which are constants in the
model that can be efficiently updated between solves.

JuMP implements parameters by a decision variable constrained on creation to the
[`Parameter`](@ref) set.

```jldoctest
julia> model = Model();

julia> @variable(model, x);

julia> @variable(model, p[i = 1:2] in Parameter(i))
2-element Vector{VariableRef}:
 p[1]
 p[2]

julia> parameter_value(p[1])
1.0

julia> set_parameter_value(p[1], 3.5)

julia> @objective(model, Max, log(p[1] * x + p[2]))
log(p[1]*x + p[2])
```

See [Parameters](@ref variables_parameters) for more information on how to
create and manage parameters.

Parameters are most useful when solving nonlinear models in a sequence:

```@repl
using JuMP, Ipopt
model = Model(Ipopt.Optimizer);
set_silent(model)
@variable(model, x)
@variable(model, p in Parameter(1.0))
@objective(model, Min, (x - p)^2)
optimize!(model)
value(x)
set_parameter_value(p, 5.0)
optimize!(model)
value(x)
```

Using parameters can be faster than creating a new model from scratch with
updated data because JuMP is able to avoid repeating a number of steps in
processing the model before handing it off to the solver.

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

## Automatic differentiation

JuMP computes first- and second-order derivatives using sparse reverse-mode
automatic differentiation. For details, see [ReverseAD](@ref).

For a tutorial on how to construct and query the derivatives, see
[Computing Hessians](@ref)

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
           op_or(op_strictly_less_than(x, -1), op_greater_than_or_equal_to(x, 1)),
           x^2,
           0.0,
       )
ifelse((x < -1) || (x >= 1), x², 0.0)
```

The available functions are:

| JuMP function                         | Julia function |
| :------------------------------------ | :------------- |
| [`op_ifelse`](@ref)                   | `ifelse`       |
| [`op_and`](@ref)                      | `&&`           |
| [`op_or`](@ref)                       | `\|\|`         |
| [`op_greater_than_or_equal_to`](@ref) | `>=`           |
| [`op_less_than_or_equal_to`](@ref)    | `<=`           |
| [`op_equal_to`](@ref)                 | `==`           |
| [`op_strictly_greater_than`](@ref)    | `>`            |
| [`op_strictly_less_than`](@ref)       | `<`            |

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

### Forcing nonlinear expressions

The JuMP macros and operator overloading will preferentially build affine ([`GenericAffExpr`](@ref)) and quadratic ([`GenericQuadExpr`](@ref)) expressions
instead of [`GenericNonlinearExpr`](@ref). For example:
```jldoctest force_nonlinear
julia> model = Model();

julia> @variable(model, x);

julia> f = (x - 0.1)^2
x² - 0.2 x + 0.010000000000000002

julia> typeof(f)
QuadExpr (alias for GenericQuadExpr{Float64, GenericVariableRef{Float64}})
```
To override this behavior, use the [`@force_nonlinear`](@ref) macro:
```jldoctest force_nonlinear
julia> g = @force_nonlinear((x - 0.1)^2)
(x - 0.1) ^ 2

julia> typeof(g)
NonlinearExpr (alias for GenericNonlinearExpr{GenericVariableRef{Float64}})
```

!!! warning
    Use this macro only if necessary. See the docstring of [`@force_nonlinear`](@ref)
    for more details on when you should use it.

## Function tracing

Nonlinear expressions can be constructed using _function tracing_. Function
tracing is when you call a regular Julia function with JuMP variables as
arguments and the function builds a nonlinear expression via operator
overloading. For example:

```@repl
using JuMP
model = Model();
@variable(model, x[1:2]);
f(x::Vector{VariableRef}) = 2 * sin(x[1]^2) + sqrt(x[2])
y = f(x)
typeof(y)
@objective(model, Max, f(x))
```

Function tracing supports functions which return vectors or arrays of
[`NonlinearExpr`](@ref):

```@repl
using JuMP
model = Model();
@variable(model, x[1:2]);
f(x::Vector{VariableRef}) = sqrt.(x)
y = f(x)
typeof(y)
@constraint(model, f(x) .<= 2)
@objective(model, Max, sum(f(x)))
```

Because function tracing uses operator overloading, there are many functions for
which it will not work. For example:

```jldoctest
julia> using JuMP

julia> model = Model();

julia> @variable(model, x[1:2]);

julia> f(x::Vector{VariableRef}) = x[1] > 1 ? 0 : x[2]
f (generic function with 1 method)

julia> f(x)
ERROR: Cannot evaluate `>` between a variable and a number.
[...]
```

In these cases, you should define a [User-defined operator](@ref jump_user_defined_operators)
using the [`@operator`](@ref) macro.

## [User-defined operators](@id jump_user_defined_operators)

In addition to a standard list of univariate and multivariate operators
recognized by the `MOI.Nonlinear` submodule, JuMP supports user-defined
operators, which let you represent nonlinear functions that cannot (or should
not) be traced, for example, because they rely on non-Julia subroutines.

!!! warning
    User-defined operators must return a scalar output. For a work-around, see
    [User-defined operators with vector outputs](@ref).

### Add an operator

Add a user-defined operator using the [`@operator`](@ref) macro:

```@repl
using JuMP
square(x) = x^2
f(x, y) = (x - 1)^2 + (y - 2)^2
model = Model();
@operator(model, op_square, 1, square)
@operator(model, op_f, 2, f)
@variable(model, x[1:2]);
@objective(model, Min, op_f(x[1], op_square(x[2])))
```

The arguments to [`@operator`](@ref) are:

 1. The model to which the operator is added.
 2. A Julia symbol object which serves as the name of the user-defined operator
    in JuMP expressions. This name must not be the same as that of the function.
 3. The number of scalar input arguments that the function takes.
 4. A Julia method which computes the function.

!!! warning
    User-defined operators cannot be deleted.

You can obtain a reference to the operator using the `model[:key]` syntax:

```@repl
using JuMP
square(x) = x^2
model = Model();
@operator(model, op_square, 1, square)
op_square_2 = model[:op_square]
```

### Automatic differentiation

JuMP computes first- and second-order derivatives of expressions using
[ReverseAD](@ref), which implements sparse reverse-mode automatic
differentiation. However, because [ReverseAD](@ref) requires the algebraic
expression as input, JuMP cannot use [ReverseAD](@ref) to differentiate
user-defined operators.

Instead, unless [Gradients and Hessians](@ref) are explicitly provided,
user-defined operators must support automatic differentiation by
[ForwardDiff.jl](https://github.com/JuliaDiff/ForwardDiff.jl).

The use of FowardDiff.jl has two important implications:

 1. ForwardDiff.jl supports only a limited subset of Julia. If you encounter an
    error adding the operator, see [Common mistakes when writing a user-defined operator](@ref).
 2. Differentiating operators with many arguments is slow. In general, you
    should try to keep the number of arguments to less than 100, and ideally, to
    less than 10.

Because of the use of ForwardDiff, in most cases, you should prefer to use
function tracing instead of defining a user-defined operator.

### Add an operator without macros

The [`@operator`](@ref) macro is syntactic sugar for [`add_nonlinear_operator`](@ref).
Thus, the non-macro version of the preceding example is:

```@repl
using JuMP
square(x) = x^2
f(x, y) = (x - 1)^2 + (y - 2)^2
model = Model();
op_square = add_nonlinear_operator(model, 1, square; name = :op_square)
model[:op_square] = op_square
op_f = add_nonlinear_operator(model, 2, f; name = :op_f)
model[:op_f] = op_f
@variable(model, x[1:2]);
@objective(model, Min, op_f(x[1], op_square(x[2])))
```

### Operators with the same name as an existing function

A common error encountered is the following:
```jldoctest nonlinear_invalid_redefinition
julia> using JuMP

julia> model = Model();

julia> f(x) = x^2
f (generic function with 1 method)

julia> @operator(model, f, 1, f)
ERROR: Unable to add the nonlinear operator `:f` with the same name as
an existing function.
[...]
```
This error occurs because `@operator(model, f, 1, f)` is equivalent to:
```julia
julia> f = add_nonlinear_operator(model, 1, f; name = :f)
```
but `f` already exists as a Julia function.

If you evaluate the function without adding it as an operator, JuMP will trace
the function using operator overloading:
```jldoctest nonlinear_invalid_redefinition
julia> @variable(model, x);

julia> f(x)
x²
```

To force JuMP to treat `f` as a user-defined operator and not trace it, add
the operator using [`add_nonlinear_operator`](@ref) and define a new method
which manually creates a [`NonlinearExpr`](@ref):
```jldoctest nonlinear_invalid_redefinition
julia> _ = add_nonlinear_operator(model, 1, f; name = :f)
NonlinearOperator(f, :f)

julia> f(x::AbstractJuMPScalar) = NonlinearExpr(:f, Any[x])
f (generic function with 2 methods)

julia> @expression(model, log(f(x)))
log(f(x))
```

### Gradients and Hessians

By default, JuMP will use automatic differentiation to compute the gradient and
Hessian of user-defined operators. If your function is not amenable to the
default automatic differentiation, or you can compute analytic derivatives, you
may pass additional arguments to [`@operator`](@ref) to compute the first- and
second-derivatives.

!!! tip
    The tutorial [Automatic differentiation of user-defined operators](@ref)
    has examples of how to use third-party Julia packages to compute automatic
    derivatives.

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
@operator(model, op_f, 1, f, ∇f, ∇²f)  # Providing ∇²f is optional
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
@operator(model, rosenbrock, 2, f, ∇f, ∇²f)  # Providing ∇²f is optional
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

### User-defined operators with vector inputs

User-defined operators which take vectors as input arguments (for example,
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
@operator(model, op_f, 5, (x...) -> f(collect(x)))
@objective(model, Min, op_f(x...))
```

If the operator takes several vector inputs, write a function that takes the
splatted arguments and reconstructs the required vector inputs:
```@repl
using JuMP
model = Model();
@variable(model, x[1:2]);
@variable(model, y[1:2]);
@variable(model, z);
f(x::Vector, y::Vector, z) = sum((x[i] * y[i])^z for i in 1:2)
f(x, y, z)
f_splat(args...) = f(collect(args[1:2]), collect(args[3:4]), args[5])
f_splat(x..., y..., z)
@operator(model, op_f, 5, f_splat)
@objective(model, Min, op_f(x..., y..., z))
```

### Common mistakes when writing a user-defined operator

JuMP uses [ForwardDiff.jl](https://github.com/JuliaDiff/ForwardDiff.jl) to
compute the first-order derivatives of user-defined operators. ForwardDiff has a
number of limitations that you should be aware of when writing user-defined
operators.

The rest of this section provides debugging advice and explains some common
mistakes.

!!! warning
    Get an error like `No method matching Float64(::ForwardDiff.Dual)`? Read
    this section.

#### Debugging

If you add an operator that does not support ForwardDiff, a long error message
will be printed. You can review the stacktrace for more information, but it can
often be hard to understand why and where your function is failing.

It may be helpful to debug the operator outside of JuMP as follows.

If the operator is univariate, do:
```jldoctest
julia> import ForwardDiff

julia> my_operator(a) = a^2
my_operator (generic function with 1 method)

julia> ForwardDiff.derivative(my_operator, 1.0)
2.0
```

If the operator is multivariate, do:
```jldoctest
julia> import ForwardDiff

julia> my_operator(a, b) = a^2 + b^2
my_operator (generic function with 1 method)

julia> ForwardDiff.gradient(x -> my_operator(x...), [1.0, 2.0])
2-element Vector{Float64}:
 2.0
 4.0
```
Note that even though the operator takes the splatted arguments,
`ForwardDiff.gradient` requires a vector as input.

#### Operator calls something unsupported by ForwardDiff

ForwardDiff works by overloading many Julia functions for a special type
`ForwardDiff.Dual <: Real`. If your operator attempts to call a function for
which an overload has not been defined, a `MethodError` will be thrown.

For example, your operator cannot call external C functions, or be the optimal
objective value of a JuMP model.

```jldoctest
julia> import ForwardDiff

julia> my_operator_bad(x) = @ccall sqrt(x::Cdouble)::Cdouble
my_operator_bad (generic function with 1 method)

julia> ForwardDiff.derivative(my_operator_bad, 1.0)
ERROR: MethodError: no method matching Float64(::ForwardDiff.Dual{ForwardDiff.Tag{typeof(my_operator_bad), Float64}, Float64, 1})
[...]
```

Unfortunately, the list of calls supported by ForwardDiff is too large to
enumerate what is an isn't allowed, so the best advice is to try and see if it
works.

#### Operator does not accept splatted input

The operator takes `f(x::Vector)` as input, instead of the splatted `f(x...)`.

```jldoctest
julia> import ForwardDiff

julia> my_operator_bad(x::Vector) = sum(x[i]^2 for i in eachindex(x))
my_operator_bad (generic function with 1 method)

julia> my_operator_good(x...) = sum(x[i]^2 for i in eachindex(x))
my_operator_good (generic function with 1 method)

julia> ForwardDiff.gradient(x -> my_operator_bad(x...), [1.0, 2.0])
ERROR: MethodError: no method matching my_operator_bad(::ForwardDiff.Dual{ForwardDiff.Tag{var"#5#6", Float64}, Float64, 2}, ::ForwardDiff.Dual{ForwardDiff.Tag{var"#5#6", Float64}, Float64, 2})
[...]

julia> ForwardDiff.gradient(x -> my_operator_good(x...), [1.0, 2.0])
2-element Vector{Float64}:
 2.0
 4.0
```

#### Operator assumes `Float64` as input

The operator assumes `Float64` will be passed as input, but it must work for any
generic `Real` type.

```jldoctest
julia> import ForwardDiff

julia> my_operator_bad(x::Float64...) = sum(x[i]^2 for i in eachindex(x))
my_operator_bad (generic function with 1 method)

julia> my_operator_good(x::Real...) = sum(x[i]^2 for i in eachindex(x))
my_operator_good (generic function with 1 method)

julia> ForwardDiff.gradient(x -> my_operator_bad(x...), [1.0, 2.0])
ERROR: MethodError: no method matching my_operator_bad(::ForwardDiff.Dual{ForwardDiff.Tag{var"#5#6", Float64}, Float64, 2}, ::ForwardDiff.Dual{ForwardDiff.Tag{var"#5#6", Float64}, Float64, 2})
[...]

julia> ForwardDiff.gradient(x -> my_operator_good(x...), [1.0, 2.0])
2-element Vector{Float64}:
 2.0
 4.0
```

#### Operator allocates `Float64` storage

The operator allocates temporary storage using `zeros(3)` or similar. This
defaults to `Float64`, so use `zeros(T, 3)` instead.

```julia
julia> import ForwardDiff

julia> function my_operator_bad(x::Real...)
           # This line is problematic. zeros(n) is short for zeros(Float64, n)
           y = zeros(length(x))
           for i in eachindex(x)
               y[i] = x[i]^2
           end
           return sum(y)
       end
my_operator_bad (generic function with 1 method)

julia> function my_operator_good(x::T...) where {T<:Real}
           y = zeros(T, length(x))
           for i in eachindex(x)
               y[i] = x[i]^2
           end
           return sum(y)
       end
my_operator_good (generic function with 1 method)

julia> ForwardDiff.gradient(x -> my_operator_bad(x...), [1.0, 2.0])
ERROR: MethodError: no method matching Float64(::ForwardDiff.Dual{ForwardDiff.Tag{var"#1#2", Float64}, Float64, 2})
[...]

julia> ForwardDiff.gradient(x -> my_operator_good(x...), [1.0, 2.0])
2-element Vector{Float64}:
 2.0
 4.0
```

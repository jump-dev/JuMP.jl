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

There are three main changes to solve nonlinear programs in JuMP.
 * Use [`@NLobjective`](@ref) instead of [`@objective`](@ref)
 * Use [`@NLconstraint`](@ref) instead of [`@constraint`](@ref)
 * Use [`@NLexpression`](@ref) instead of [`@expression`](@ref)

!!! info
    There are some restrictions on what syntax you can use in the `@NLxxx`
    macros. Make sure to read the [Syntax notes](@ref).

## Set a nonlinear objective

Use [`@NLobjective`](@ref) to set a nonlinear objective.

```jldoctest
julia> model = Model();

julia> @variable(model, x[1:2]);

julia> @NLobjective(model, Min, exp(x[1]) - sqrt(x[2]))
```
To modify a nonlinear objective, call [`@NLobjective`](@ref) again.

## Add a nonlinear constraint

Use [`@NLconstraint`](@ref) to add a nonlinear constraint.

```jldoctest nonlinear_constraint
julia> model = Model();

julia> @variable(model, x[1:2]);

julia> @NLconstraint(model, exp(x[1]) <= 1)
exp(x[1]) - 1.0 ≤ 0

julia> @NLconstraint(model, [i = 1:2], x[i]^i >= i)
2-element Vector{NonlinearConstraintRef{ScalarShape}}:
 x[1] ^ 1.0 - 1.0 ≥ 0
 x[2] ^ 2.0 - 2.0 ≥ 0

julia> @NLconstraint(model, con[i = 1:2], prod(x[j] for j = 1:i) == i)
2-element Vector{NonlinearConstraintRef{ScalarShape}}:
 (*)(x[1]) - 1.0 = 0
 x[1] * x[2] - 2.0 = 0
```

!!! info
    You can only create nonlinear constraints with `<=`, `>=`, and `==`.
    More general `Nonlinear`-in-`Set` constraints are not supported.

Delete a nonlinear constraint using [`delete`](@ref):
```jldoctest nonlinear_constraint
julia> delete(model, con[1])
```

## Create a nonlinear expression

Use [`@NLexpression`](@ref) to create nonlinear expression objects. The syntax
is identical to [`@expression`](@ref), except that the expression can contain
nonlinear terms.

```jldoctest nl_expression
julia> model = Model();

julia> @variable(model, x[1:2]);

julia> expr = @NLexpression(model, exp(x[1]) + sqrt(x[2]))
subexpression[1]: exp(x[1]) + sqrt(x[2])

julia> my_anon_expr = @NLexpression(model, [i = 1:2], sin(x[i]))
2-element Vector{NonlinearExpression}:
 subexpression[2]: sin(x[1])
 subexpression[3]: sin(x[2])

julia> @NLexpression(model, my_expr[i = 1:2], sin(x[i]))
2-element Vector{NonlinearExpression}:
 subexpression[4]: sin(x[1])
 subexpression[5]: sin(x[2])
```

Nonlinear expression can be used in [`@NLobjective`](@ref), [`@NLconstraint`](@ref),
and even nested in other [`@NLexpression`](@ref)s.

```jldoctest nl_expression
julia> @NLobjective(model, Min, expr^2 + 1)

julia> @NLconstraint(model, [i = 1:2], my_expr[i] <= i)
2-element Vector{NonlinearConstraintRef{ScalarShape}}:
 subexpression[4] - 1.0 ≤ 0
 subexpression[5] - 2.0 ≤ 0

julia> @NLexpression(model, nested[i = 1:2], sin(my_expr[i]))
2-element Vector{NonlinearExpression}:
 subexpression[6]: sin(subexpression[4])
 subexpression[7]: sin(subexpression[5])
```

Use [`value`](@ref) to query the value of a nonlinear expression:
```jldoctest nl_expression
julia> set_start_value(x[1], 1.0)

julia> value(start_value, nested[1])
0.7456241416655579

julia> sin(sin(1.0))
0.7456241416655579
```

## Create a nonlinear parameter

For nonlinear models only, JuMP offers a syntax for explicit "parameter" objects,
which are constants in the model that can be efficiently updated between solves.

Nonlinear parameters are declared by using the [`@NLparameter`](@ref) macro
and may be indexed by arbitrary sets analogously to JuMP variables and
expressions.

The initial value of the parameter must be provided on the right-hand side of
the `==` sign.

```jldoctest nonlinear_parameters
julia> model = Model();

julia> @variable(model, x);

julia> @NLparameter(model, p[i = 1:2] == i)
2-element Vector{NonlinearParameter}:
 parameter[1] == 1.0
 parameter[2] == 2.0
```

Create anonymous parameters using the `value` keyword:
```jldoctest nonlinear_parameters
julia> anon_parameter = @NLparameter(model, value = 1)
parameter[3] == 1.0
```

!!! info
    A parameter is not an optimization variable. It must be fixed to a value with
    `==`. If you want a parameter that is `<=` or `>=`, create a variable instead
    using [`@variable`](@ref).

Use [`value`](@ref) and [`set_value`](@ref) to query or update the value of a
parameter.

```jldoctest nonlinear_parameters
julia> value.(p)
2-element Vector{Float64}:
 1.0
 2.0

julia> set_value(p[2], 3.0)
3.0

julia> value.(p)
2-element Vector{Float64}:
 1.0
 3.0
```

Nonlinear parameters can be used *within nonlinear macros* only:

```jldoctest nonlinear_parameters
julia> @objective(model, Max, p[1] * x)
ERROR: MethodError: no method matching *(::NonlinearParameter, ::VariableRef)
[...]

julia> @NLobjective(model, Max, p[1] * x)

julia> @expression(model, my_expr, p[1] * x^2)
ERROR: MethodError: no method matching *(::NonlinearParameter, ::QuadExpr)
Closest candidates are:
[...]

julia> @NLexpression(model, my_nl_expr, p[1] * x^2)
subexpression[1]: parameter[1] * x ^ 2.0
```

### When to use a parameter

Nonlinear parameters are useful when solving nonlinear models in a sequence:

```@example
using JuMP, Ipopt
model = Model(Ipopt.Optimizer)
set_silent(model)
@variable(model, z)
@NLparameter(model, x == 1.0)
@NLobjective(model, Min, (z - x)^2)
optimize!(model)
@show value(z) # Equals 1.0.

# Now, update the value of x to solve a different problem.
set_value(x, 5.0)
optimize!(model)
@show value(z) # Equals 5.0
nothing #hide
```

!!! info
    Using nonlinear parameters can be faster than creating a new model from
    scratch with updated data because JuMP is able to avoid repeating a number
    of steps in processing the model before handing it off to the solver.

## Syntax notes

The syntax accepted in nonlinear macros is more restricted than the syntax
for linear and quadratic macros. We note some important points below.

### No operator overloading

There is no operator overloading provided to build up nonlinear expressions.
For example, if `x` is a JuMP variable, the code `3x` will return an
`AffExpr` object that can be used inside of future expressions and linear
constraints. However, the code `sin(x)` is an error. All nonlinear
expressions must be inside of macros.

```jldoctest
julia> model = Model();

julia> @variable(model, x);

julia> expr = sin(x) + 1
ERROR: sin is not defined for type AbstractVariableRef. Are you trying to build a nonlinear problem? Make sure you use @NLconstraint/@NLobjective.
[...]

julia> expr = @NLexpression(model, sin(x) + 1)
subexpression[1]: sin(x) + 1.0
```

### Scalar operations only

Except for the splatting syntax discussed below, all expressions
must be simple scalar operations. You cannot use `dot`, matrix-vector products,
vector slices, etc.
```jldoctest nlp_scalar_only
julia> model = Model();

julia> @variable(model, x[1:2]);

julia> @variable(model, y);

julia> c = [1, 2];

julia> @NLobjective(model, Min, c' * x + 3y)
ERROR: Unexpected array [1 2] in nonlinear expression. Nonlinear expressions may contain only scalar expressions.
[...]
```

Translate vector operations into explicit `sum()` operations:
```jldoctest nlp_scalar_only
julia> @NLobjective(model, Min, sum(c[i] * x[i] for i = 1:2) + 3y)
```

Or use an [`@expression`](@ref):
```jldoctest nlp_scalar_only
julia> @expression(model, expr, c' * x)
x[1] + 2 x[2]

julia> @NLobjective(model, Min, expr + 3y)

```

### Splatting

The [splatting operator](https://docs.julialang.org/en/v1/manual/faq/#...-splits-one-argument-into-many-different-arguments-in-function-calls-1)
  `...` is recognized in a very restricted setting for expanding function
  arguments. The expression splatted can be *only* a symbol. More complex
  expressions are not recognized.

```jldoctest; filter=r"≤|<="
julia> model = Model();

julia> @variable(model, x[1:3]);

julia> @NLconstraint(model, *(x...) <= 1.0)
x[1] * x[2] * x[3] - 1.0 ≤ 0

julia> @NLconstraint(model, *((x / 2)...) <= 0.0)
ERROR: Unsupported use of the splatting operator. JuMP supports splatting only symbols. For example, `x...` is ok, but `(x + 1)...`, `[x; y]...` and `g(f(y)...)` are not.
```

## User-defined Functions

JuMP natively supports the set of univariate and multivariate functions recognized by the
`MOI.Nonlinear` submodule. In addition to this list of functions, it is possible
to register custom *user-defined* nonlinear functions. User-defined functions
can be used anywhere in [`@NLobjective`](@ref), [`@NLconstraint`](@ref), and
[`@NLexpression`](@ref).

JuMP will attempt to automatically register functions it detects in your
nonlinear expressions, which usually means manually registering a function is
not needed. Two exceptions are if you want to provide custom derivatives, or if
the function is not available in the scope of the nonlinear expression.

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
automatic differentiation, use the [`register`](@ref) method as in the following
example:

```@example
using JuMP #hide
square(x) = x^2
f(x, y) = (x - 1)^2 + (y - 2)^2

model = Model()

register(model, :square, 1, square; autodiff = true)
register(model, :my_f, 2, f; autodiff = true)

@variable(model, x[1:2] >= 0.5)
@NLobjective(model, Min, my_f(x[1], square(x[2])))
```

The above code creates a JuMP model with the objective function
`(x[1] - 1)^2 + (x[2]^2 - 2)^2`. The arguments to [`register`](@ref) are:
 1. The model for which the functions are registered.
 2. A Julia symbol object which serves as the name of the user-defined function
    in JuMP expressions.
 3. The number of input arguments that the function takes.
 4. The Julia method which computes the function
 5. A flag to instruct JuMP to compute exact gradients automatically.

!!! tip
    The symbol `:my_f` doesn't have to match the name of the function `f`.
    However, it's more readable if it does. Make sure you use `my_f`
    and not `f` in the macros.

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
register(model, :my_square, 1, f, ∇f; autodiff = true)
@variable(model, x >= 0)
@NLobjective(model, Min, my_square(x))
```
If `autodiff = true`, JuMP will use automatic differentiation to compute the
hessian.

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
register(model, :my_square, 2, f, ∇f)
@variable(model, x[1:2] >= 0)
@NLobjective(model, Min, my_square(x[1], x[2]))
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
register(model, :my_square, 1, f, ∇f, ∇²f)
@variable(model, x >= 0)
@NLobjective(model, Min, my_square(x))
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
register(model, :rosenbrock, 2, f, ∇f, ∇²f)
@variable(model, x[1:2])
@NLobjective(model, Min, rosenbrock(x[1], x[2]))
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
register(model, :f, 5, f; autodiff = true)
@NLobjective(model, Min, f(x...))
```

!!! tip
    Make sure to read the syntax restrictions of [Splatting](@ref).

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

## Querying derivatives from a JuMP model

For some advanced use cases, one may want to directly query the derivatives of a
JuMP model instead of handing the problem off to a solver.
Internally, JuMP implements the [`MOI.AbstractNLPEvaluator`](@ref) interface. To
obtain an NLP evaluator object from a JuMP model, use [`NLPEvaluator`](@ref).
[`index`](@ref) returns the [`MOI.VariableIndex`](@ref) corresponding to a JuMP
variable. `MOI.VariableIndex` itself is a type-safe wrapper for `Int64` (stored
in the `.value` field.)

For example:

```jldoctest derivatives
julia> raw_index(v::MOI.VariableIndex) = v.value
raw_index (generic function with 1 method)

julia> model = Model();

julia> @variable(model, x)
x

julia> @variable(model, y)
y

julia> @NLobjective(model, Min, sin(x) + sin(y))

julia> values = zeros(2)
2-element Vector{Float64}:
 0.0
 0.0

julia> x_index = raw_index(JuMP.index(x))
1

julia> y_index = raw_index(JuMP.index(y))
2

julia> values[x_index] = 2.0
2.0

julia> values[y_index] = 3.0
3.0

julia> d = NLPEvaluator(model)
Nonlinear.Evaluator with available features:
  * :Grad
  * :Jac
  * :JacVec
  * :Hess
  * :HessVec
  * :ExprGraph

julia> MOI.initialize(d, [:Grad])

julia> MOI.eval_objective(d, values)
1.0504174348855488

julia> sin(2.0) + sin(3.0)
1.0504174348855488

julia> ∇f = zeros(2)
2-element Vector{Float64}:
 0.0
 0.0

julia> MOI.eval_objective_gradient(d, ∇f, values)

julia> ∇f[x_index], ∇f[y_index]
(-0.4161468365471424, -0.9899924966004454)

julia> cos(2.0), cos(3.0)
(-0.4161468365471424, -0.9899924966004454)
```

Only nonlinear constraints (those added with [`@NLconstraint`](@ref)), and
nonlinear objectives (added with [`@NLobjective`](@ref)) exist in the scope of
the [`NLPEvaluator`](@ref).

The [`NLPEvaluator`](@ref) *does not evaluate derivatives of linear or quadratic
constraints or objectives*.

The [`index`](@ref) method applied to a nonlinear constraint reference object
returns its index as a [`MOI.Nonlinear.ConstraintIndex`](@ref). For example:

```jldoctest
julia> model = Model();

julia> @variable(model, x);

julia> @NLconstraint(model, cons1, sin(x) <= 1);

julia> @NLconstraint(model, cons2, x + 5 == 10);

julia> typeof(cons1)
NonlinearConstraintRef{ScalarShape} (alias for ConstraintRef{GenericModel{Float64}, MathOptInterface.Nonlinear.ConstraintIndex, ScalarShape})

julia> index(cons1)
MathOptInterface.Nonlinear.ConstraintIndex(1)

julia> index(cons2)
MathOptInterface.Nonlinear.ConstraintIndex(2)
```

```@meta
# TODO: Provide a link for how to access the linear and quadratic parts of the
# model.
```

Note that for one-sided nonlinear constraints, JuMP subtracts any values on the
right-hand side when computing expressions. In other words, one-sided nonlinear
constraints are always transformed to have a right-hand side of zero.

This method of querying derivatives directly from a JuMP model is convenient for
interacting with the model in a structured way, for example, for accessing derivatives
of specific variables. For example, in statistical maximum likelihood estimation
problems, one is often interested in the Hessian matrix at the optimal solution,
which can be queried using the [`NLPEvaluator`](@ref).

## Raw expression input

!!! warning
    This section requires advanced knowledge of Julia's `Expr`. You should read
    the [Expressions and evaluation](https://docs.julialang.org/en/v1/manual/metaprogramming/#Expressions-and-evaluation)
    section of the Julia documentation first.

In addition to the [`@NLexpression`](@ref), [`@NLobjective`](@ref) and
[`@NLconstraint`](@ref) macros, it is also possible to provide Julia `Expr`
objects directly by using [`add_nonlinear_expression`](@ref),
[`set_nonlinear_objective`](@ref) and [`add_nonlinear_constraint`](@ref).

This input form may be useful if the expressions are generated programmatically,
or if you experience compilation issues with the macro input (see
[Known performance issues](@ref) for more information).

### Add a nonlinear expression

Use [`add_nonlinear_expression`](@ref) to add a nonlinear expression to the model.

```jldoctest
julia> model = Model();

julia> @variable(model, x)
x

julia> expr = :($(x) + sin($(x)^2))
:(x + sin(x ^ 2))

julia> expr_ref = add_nonlinear_expression(model, expr)
subexpression[1]: x + sin(x ^ 2.0)
```
This is equivalent to
```jldoctest
julia> model = Model();

julia> @variable(model, x);

julia> expr_ref = @NLexpression(model, x + sin(x^2))
subexpression[1]: x + sin(x ^ 2.0)
```

!!! note
    You must interpolate the variables directly into the expression `expr`.

### Set the objective function

Use [`set_nonlinear_objective`](@ref) to set a nonlinear objective.

```jldoctest
julia> model = Model();

julia> @variable(model, x);

julia> expr = :($(x) + $(x)^2)
:(x + x ^ 2)

julia> set_nonlinear_objective(model, MIN_SENSE, expr)
```
This is equivalent to
```jldoctest
julia> model = Model();

julia> @variable(model, x);

julia> @NLobjective(model, Min, x + x^2)
```

!!! note
    You must use `MIN_SENSE` or `MAX_SENSE` instead of `Min` and `Max`.

### Add a constraint

Use [`add_nonlinear_constraint`](@ref) to add a nonlinear constraint.

```jldoctest
julia> model = Model();

julia> @variable(model, x);

julia> expr = :($(x) + $(x)^2)
:(x + x ^ 2)

julia> add_nonlinear_constraint(model, :($(expr) <= 1))
(x + x ^ 2.0) - 1.0 ≤ 0
```

This is equivalent to
```jldoctest
julia> model = Model();

julia> @variable(model, x);

julia> @NLconstraint(model, Min, x + x^2 <= 1)
(x + x ^ 2.0) - 1.0 ≤ 0
```

### More complicated examples

Raw expression input is most useful when the expressions are generated
programmatically, often in conjunction with user-defined functions.

As an example, we construct a model with the nonlinear constraints `f(x) <= 1`,
where `f(x) = x^2` and `f(x) = sin(x)^2`:
```jldoctest
julia> function main(functions::Vector{Function})
           model = Model()
           @variable(model, x)
           for (i, f) in enumerate(functions)
               f_sym = Symbol("f_$(i)")
               register(model, f_sym, 1, f; autodiff = true)
               add_nonlinear_constraint(model, :($(f_sym)($(x)) <= 1))
           end
           print(model)
           return
       end
main (generic function with 1 method)

julia> main([x -> x^2, x -> sin(x)^2])
Feasibility
Subject to
 f_1(x) - 1.0 ≤ 0
 f_2(x) - 1.0 ≤ 0
```

As another example, we construct a model with the constraint
`x^2 + sin(x)^2 <= 1`:
```jldoctest
julia> function main(functions::Vector{Function})
           model = Model()
           @variable(model, x)
           expr = Expr(:call, :+)
           for (i, f) in enumerate(functions)
               f_sym = Symbol("f_$(i)")
               register(model, f_sym, 1, f; autodiff = true)
               push!(expr.args, :($(f_sym)($(x))))
           end
           add_nonlinear_constraint(model, :($(expr) <= 1))
           print(model)
           return
       end
main (generic function with 1 method)

julia> main([x -> x^2, x -> sin(x)^2])
Feasibility
Subject to
 (f_1(x) + f_2(x)) - 1.0 ≤ 0
```

### Registered functions with a variable number of arguments

User defined functions require a fixed number of input arguments. However,
sometimes you will want to use a registered function like:
```jldoctest nlp_register_variable_arguments
julia> f(x...) = sum(exp(x[i]^2) for i in 1:length(x));
```
with different numbers of arguments.

The solution is to register the same function `f` for each unique number of
input arguments, making sure to use a unique name each time. For example:

```jldoctest nlp_register_variable_arguments
julia> A = [[1], [1, 2], [2, 3, 4], [1, 3, 4, 5]];

julia> model = Model();

julia> @variable(model, x[1:5]);

julia> funcs = Set{Symbol}();

julia> for a in A
           key = Symbol("f$(length(a))")
           if !(key in funcs)
               push!(funcs, key)
               register(model, key, length(a), f; autodiff = true)
           end
           add_nonlinear_constraint(model, :($key($(x[a]...)) <= 1))
       end

julia> print(model)
Feasibility
Subject to
 f1(x[1]) - 1.0 ≤ 0
 f2(x[1], x[2]) - 1.0 ≤ 0
 f3(x[2], x[3], x[4]) - 1.0 ≤ 0
 f4(x[1], x[3], x[4], x[5]) - 1.0 ≤ 0
```

## Known performance issues

The macro-based input to JuMP's nonlinear interface can cause a performance
issue if you:

 1. write a macro with a large number (hundreds) of terms
 2. call that macro from within a function instead of from the top-level in
    global scope.

The first issue does not depend on the number of resulting terms in the
mathematical expression, but rather the number of terms in the Julia `Expr`
representation of that expression. For example, the expression
`sum(x[i] for i in 1:1_000_000)` contains one million mathematical terms, but
the `Expr` representation is just a single sum.

The most common cause, other than a lot of tedious typing, is if you write a
program that automatically writes a JuMP model as a text file, which you later
execute. One example is [MINLPlib.jl](https://github.com/lanl-ansi/MINLPLib.jl)
which automatically transpiled models in the GAMS scalar format into JuMP
examples.

As a rule of thumb, if you are writing programs to automatically generate
expressions for the JuMP macros, you should target the [Raw expression input](@ref)
instead. For more information, read [MathOptInterface Issue#1997](https://github.com/jump-dev/MathOptInterface.jl/issues/1997).

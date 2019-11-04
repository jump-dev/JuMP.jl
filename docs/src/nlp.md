# Nonlinear Modeling

```@meta
CurrentModule = JuMP
DocTestSetup = quote
    using JuMP
end
```

JuMP has support for general smooth nonlinear (convex and nonconvex)
optimization problems. JuMP is able to provide exact, sparse second-order
derivatives to solvers. This information can improve solver accuracy and
performance.

Nonlinear objectives and constraints are specified by using the `@NLobjective`
and `@NLconstraint` macros. The familiar `sum()` syntax is supported within
these macros, as well as `prod()` which analogously represents the product of
the terms within. Note that the `@objective` and `@constraint` macros (and
corresponding functions) do *not* currently support nonlinear expressions.
However, a model can contain a mix of linear, quadratic, and nonlinear
contraints or objective functions. Starting points may be provided by using the
`start` keyword argument to `@variable`.

For example, we can solve the classical Rosenbrock problem (with a twist) as
follows:

```julia
using Ipopt
model = Model(Ipopt.Optimizer)
@variable(model, x, start = 0.0)
@variable(model, y, start = 0.0)

@NLobjective(model, Min, (1 - x)^2 + 100 * (y - x^2)^2)

optimize!(model)
println("x = ", value(x), " y = ", value(y))

# adding a (linear) constraint
@constraint(model, x + y == 10)
optimize!(model)
println("x = ", value(x), " y = ", value(y))
```

See the JuMP [examples directory](https://github.com/JuliaOpt/JuMP.jl/tree/bff0916a2025df64e4a0be8933b58ea7bdc5eb0b/examples)
for more examples (which include `mle.jl`, `rosenbrock.jl`, and `clnlbeam.jl`).

The [NLP solver tests](https://github.com/JuliaOpt/JuMP.jl/blob/bff0916a2025df64e4a0be8933b58ea7bdc5eb0b/test/nlp_solver.jl)
contain additional examples.

## Syntax notes

The syntax accepted in nonlinear expressions is more restricted than the syntax
for linear and quadratic expressions. We note some important points below.

- With the exception of the splatting syntax discussed below, all expressions
  must be simple scalar operations. You cannot use `dot`,
  matrix-vector products, vector slices, etc. Translate vector operations into
  explicit `sum()` operations or use the `AffExpr` plus auxiliary variable trick
  described below.
- There is no operator overloading provided to build up nonlinear expressions.
  For example, if `x` is a JuMP variable, the code `3x` will return an `AffExpr`
  object that can be used inside of future expressions and linear constraints.
  However, the code `sin(x)` is an error. All nonlinear expressions must be
  inside of macros.
- [User-defined Functions](@ref) may be used within nonlinear expressions only
  after they are registered. For example, the follow code results in an error
  because `register()` must be called first to register `my_function`.

```jldoctest
model = Model()
my_function(a, b) = exp(a) * b
@variable(model, x)
@variable(model, y)
@NLobjective(model, Min, my_function(x, y))

# output

ERROR: Unrecognized function "my_function" used in nonlinear expression.
```

- `AffExpr` and `QuadExpr` objects cannot currently be used inside nonlinear
  expressions. Instead, introduce auxiliary variables, e.g.:

```julia
    my_expr = dot(c, x) + 3y # where x and y are variables
    @variable(model, aux)
    @constraint(model, aux == my_expr)
    @NLobjective(model, Min, sin(aux))
```

- You can declare embeddable nonlinear expressions with `@NLexpression`.
  For example:

```julia
    @NLexpression(model, my_expr[i = 1:n], sin(x[i]))
    @NLconstraint(model, my_constr[i = 1:n], my_expr[i] <= 0.5)
```

- Anonymous syntax is supported in `@NLexpression` and `@NLconstraint`:

```julia
    my_expr = @NLexpression(model, [i = 1:n], sin(x[i]))
    my_constr = @NLconstraint(model, [i = 1:n], my_expr[i] <= 0.5)
```

- The [splatting operator](https://docs.julialang.org/en/v1/manual/faq/#...-splits-one-argument-into-many-different-arguments-in-function-calls-1)
  `...` is recognized in a very restricted setting for expanding function
  arguments. The expression splatted can be *only* a symbol. More complex
  expressions are not recognized.

```jldoctest; filter=r"≤|<="
julia> model = Model();

julia> @variable(model, x[1:3]);

julia> @NLconstraint(model, *(x...) <= 1.0)
x[1] * x[2] * x[3] - 1.0 ≤ 0

julia> @NLconstraint(model, *((x / 2)...) <= 0.0)
ERROR: LoadError: Unexpected expression in (*)(x / 2...). JuMP supports splatting only symbols. For example, x... is ok, but (x + 1)..., [x; y]... and g(f(y)...) are not.
```

## Nonlinear Parameters

For nonlinear models only, JuMP offers a syntax for explicit "parameter" objects
which can be used to modify a model in-place just by updating the value of the
parameter. Nonlinear parameters are declared by using the `@NLparameter` macro
and may be indexed by arbitrary sets analogously to JuMP variables and
expressions. The initial value of the parameter must be provided on the
right-hand side of the `==` sign. There is no anonymous syntax for creating
parameters.

```@docs
@NLparameter
```

You may use `value` and `set_value` to query or update the value of a parameter.
```@docs
value(::JuMP.NonlinearParameter)
set_value(::JuMP.NonlinearParameter, ::Number)
```
Nonlinear parameters can be used *within nonlinear expressions* only:

```julia
@NLparameter(model, x == 10)
@variable(model, z)
@objective(model, Max, x * z)             # Error: x is a nonlinear parameter.
@NLobjective(model, Max, x * z)           # Ok.
@expression(model, my_expr, x * z^2)      # Error: x is a nonlinear parameter.
@NLexpression(model, my_nl_expr, x * z^2) # Ok.
```

Nonlinear parameters are useful when solving nonlinear models in a sequence:

```julia
using Ipopt
model = Model(Ipopt.Optimizer)
@variable(model, z)
@NLparameter(model, x == 1.0)
@NLobjective(model, Min, (z - x)^2)
optimize!(model)
value(z) # Equals 1.0.

# Now, update the value of x to solve a different problem.
set_value(x, 5.0)
optimize!(model)
value(z) # Equals 5.0
```

Using nonlinear parameters can be faster than creating a new model from scratch
with updated data because JuMP is able to avoid repeating a number of steps in
processing the model before handing it off to the solver.

## User-defined Functions

JuMP's library of recognized univariate functions is derived from the
[Calculus.jl](https://github.com/johnmyleswhite/Calculus.jl) package. If you
encounter a standard special function not currently supported by JuMP, consider
contributing to the
[list of derivative rules](https://github.com/johnmyleswhite/Calculus.jl/blob/cb42f3699177449a42bdc3461c8aea8777aa8c39/src/differentiate.jl#L115)
there. In addition to this built-in list of functions, it is possible to
register custom (*user-defined*) nonlinear functions to use within nonlinear
expressions. JuMP does not support black-box optimization, so all user-defined
functions must provide derivatives in some form. Fortunately, JuMP supports
**automatic differentiation of user-defined functions**, a feature to our
knowledge not available in any comparable modeling systems.

!!! note
    Automatic differentiation is *not* finite differencing. JuMP's automatically
    computed derivatives are not subject to approximation error.

JuMP uses [ForwardDiff.jl](https://github.com/JuliaDiff/ForwardDiff.jl) to
perform automatic differentiation; see the ForwardDiff.jl
[documentation](http://www.juliadiff.org/ForwardDiff.jl/v0.10.2/user/limitations.html)
for a description of how to write a function suitable for automatic
differentiation.

!!! note
    If you see method errors with `ForwardDiff.Duals`, see the guidelines at
    [ForwardDiff.jl](http://www.juliadiff.org/ForwardDiff.jl/release-0.10/user/limitations.html).
    The most common error is that your user-defined function is not generic with
    respect to the number type, i.e., don't assume that the input to the function
    is `Float64`.
    ```julia
    f(x::Float64) = 2 * x  # This will not work.
    f(x::Real)    = 2 * x  # This is good.
    f(x)          = 2 * x  # This is also good.
    ```

To register a user-defined function with derivatives computed by
automatic differentiation, use the `register` method as in the following
example:

```julia
my_square(x) = x^2
my_f(x,y) = (x - 1)^2 + (y - 2)^2

model = Model()

register(model, :my_f, 2, my_f, autodiff=true)
register(model, :my_square, 1, my_square, autodiff=true)

@variable(model, x[1:2] >= 0.5)
@NLobjective(model, Min, my_f(x[1], my_square(x[2])))
```

The above code creates a JuMP model with the objective function
`(x[1] - 1)^2 + (x[2]^2 - 2)^2`. The first argument to `register` is the
model for which the functions are registered. The second argument is a Julia
symbol object which serves as the name of the user-defined function in JuMP
expressions; the JuMP name need not be the same as the name of the corresponding
Julia method. The third argument specifies how many arguments the function
takes. The fourth argument is the name of the Julia method which computes the
function, and `autodiff=true` instructs JuMP to compute exact gradients
automatically.

Forward-mode automatic differentiation as implemented by ForwardDiff.jl has a
computational cost that scales linearly with the number of input dimensions. As
such, it is not the most efficient way to compute gradients of user-defined
functions if the number of input arguments is large. In this case, users may
want to provide their own routines for evaluating gradients. The more general
syntax for `register` which accepts user-provided derivative evaluation
routines is:

```julia
JuMP.register(model::Model, s::Symbol, dimension::Integer, f::Function,
              ∇f::Function, ∇²f::Function)
```

The input differs between functions which take a single input argument and
functions which take more than one. For univariate functions, the derivative
evaluation routines should return a number which represents the first and
second-order derivatives respectively. For multivariate functions, the
derivative evaluation routines will be passed a gradient vector which they must
explicitly fill. Second-order derivatives of multivariate functions are not
currently supported; this argument should be omitted. The following example sets
up the same optimization problem as before, but now we explicitly provide
evaluation routines for the user-defined functions:

```julia
my_square(x) = x^2
my_square_prime(x) = 2x
my_square_prime_prime(x) = 2

my_f(x, y) = (x - 1)^2 + (y - 2)^2
function ∇f(g, x, y)
    g[1] = 2 * (x - 1)
    g[2] = 2 * (y - 2)
end

model = Model()

register(model, :my_f, 2, my_f, ∇f)
register(model, :my_square, 1, my_square, my_square_prime,
         my_square_prime_prime)

@variable(model, x[1:2] >= 0.5)
@NLobjective(model, Min, my_f(x[1], my_square(x[2])))
```

Once registered, user-defined functions can also be used in constraints. For
example:
```julia
@NLconstraint(model, my_square(x[1]) <= 2.0)
```

### User-defined functions with vector inputs

User-defined functions which take vectors as input arguments (e.g.
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
```julia
model = Model()
@variable(model, x[1:5] >= 0)
f(x...) = sum(x[i]^i for i in 1:length(x))
register(model, :f, 5, f; autodiff = true)
@NLobjective(model, Min, f(x...))
```

## Factors affecting solution time

The execution time when solving a nonlinear programming problem can be divided
into two parts, the time spent in the optimization algorithm (the solver) and
the time spent evaluating the nonlinear functions and corresponding derivatives.
Ipopt explicitly displays these two timings in its output, for example:

``` sourceCode
Total CPU secs in IPOPT (w/o function evaluations)   =      7.412
Total CPU secs in NLP function evaluations           =      2.083
```

For Ipopt in particular, one can improve the performance by installing advanced
sparse linear algebra packages, see [Installation Guide](@ref). For other
solvers, see their respective documentation for performance tips.

The function evaluation time, on the other hand, is the responsibility of the
modeling language. JuMP computes derivatives by using reverse-mode automatic
differentiation with graph coloring methods for exploiting sparsity of the
Hessian matrix [^1]. As a conservative bound, JuMP's performance here currently
may be expected to be within a factor of 5 of AMPL's.

## Querying derivatives from a JuMP model

For some advanced use cases, one may want to directly query the derivatives of a
JuMP model instead of handing the problem off to a solver.
Internally, JuMP implements the `AbstractNLPEvaluator` interface from
[MathOptInterface](http://www.juliaopt.org/MathOptInterface.jl/v0.9.1/apireference/#NLP-evaluator-methods-1).
To obtain an NLP evaluator object from a JuMP model, use `JuMP.NLPEvaluator`.
`JuMP.index` returns the `MOI.VariableIndex` corresponding to a JuMP variable.
`MOI.VariableIndex` itself is a type-safe wrapper for `Int64` (stored in the
`value` field.)

For example:

```jldoctest derivatives
raw_index(v::MOI.VariableIndex) = v.value
model = Model()
@variable(model, x)
@variable(model, y)
@NLobjective(model, Min, sin(x) + sin(y))
values = zeros(2)
x_index = raw_index(JuMP.index(x))
y_index = raw_index(JuMP.index(y))
values[x_index] = 2.0
values[y_index] = 3.0
d = NLPEvaluator(model)
MOI.initialize(d, [:Grad])
MOI.eval_objective(d, values) # == sin(2.0) + sin(3.0)

# output
1.0504174348855488
```

```jldoctest derivatives
∇f = zeros(2)
MOI.eval_objective_gradient(d, ∇f, values)
(∇f[x_index], ∇f[y_index]) # == (cos(2.0), cos(3.0))

# output
(-0.4161468365471424, -0.9899924966004454)
```

Only nonlinear constraints (those added with `@NLconstraint`), and nonlinear
objectives (added with `@NLobjective`) exist in the scope of the `NLPEvaluator`.
The `NLPEvaluator` *does not evaluate derivatives of linear or quadratic
constraints or objectives*. The `index` method applied to a nonlinear constraint
reference object returns its index as a `NonlinearConstraintIndex`. The `value`
field of `NonlinearConstraintIndex` stores the raw integer index. For example:

```jldoctest
julia> model = Model();

julia> @variable(model, x);

julia> @NLconstraint(model, cons1, sin(x) <= 1);

julia> @NLconstraint(model, cons2, x + 5 == 10);

julia> typeof(cons1)
ConstraintRef{Model,NonlinearConstraintIndex,ScalarShape}

julia> index(cons1)
NonlinearConstraintIndex(1)

julia> index(cons2)
NonlinearConstraintIndex(2)
```

```@meta
# TODO: Provide a link for how to access the linear and quadratic parts of the
# model.
```

Note that for one-sided nonlinear constraints, JuMP subtracts any values on the
right-hand side when computing expressions. In other words, one-sided nonlinear
constraints are always transformed to have a right-hand side of zero.

This method of querying derivatives directly from a JuMP model is convenient for
interacting with the model in a structured way, e.g., for accessing derivatives
of specific variables. For example, in statistical maximum likelihood estimation
problems, one is often interested in the Hessian matrix at the optimal solution,
which can be queried using the `NLPEvaluator`.

## Raw expression input

In addition to the `@NLobjective` and `@NLconstraint` macros, it is also
possible to provide Julia `Expr` objects directly by using
`set_NL_objective` and `add_NL_constraint`. This input form may be
useful if the expressions are generated programmatically. JuMP variables should
be spliced into the expression object. For example:

```julia
@variable(model, 1 <= x[i = 1:4] <= 5)
set_NL_objective(model, :Min, :($(x[1])*$(x[4])*($(x[1])+$(x[2])+$(x[3])) + $(x[3])))
add_NL_constraint(model, :($(x[1])*$(x[2])*$(x[3])*$(x[4]) >= 25))

# Equivalent form using traditional JuMP macros:
@NLobjective(model, Min, x[1] * x[4] * (x[1] + x[2] + x[3]) + x[3])
@NLconstraint(model, x[1] * x[2] * x[3] * x[4] >= 25)
```

See the Julia documentation for more examples and description of Julia
expressions.

## Reference

```@docs
@NLconstraint
@NLexpression
@NLobjective
```

[^1]: Dunning, Huchette, and Lubin, "JuMP: A Modeling Language for Mathematical Optimization", SIAM Review, [PDF](https://mlubin.github.io/pdf/jump-sirev.pdf).

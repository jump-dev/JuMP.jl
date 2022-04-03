# [Nonlinear](@id nonlinear_developers)

The `JuMP.Nonlinear` submodule contains data structures and functions for
working with a nonlinear program in the form of an expression tree. This page
explains the API and describes the rationale behind its design.

## Standard form

[Nonlinear programs (NLPs)](https://en.wikipedia.org/wiki/Nonlinear_programming)
are a class of optimization problems in which some of the constraints or the
objective function are nonlinear:
```math
\begin{align}
    \min_{x \in \mathbb{R}^n} & f_0(x) \\
    \;\;\text{s.t.} & l_j \le f_j(x) \le u_j & j = 1 \ldots m
\end{align}
```
There may be additional constraints, as well as things like variable bounds
and integrality restrictions, but we do not consider them here because they are
best dealt with by other components of JuMP and MathOptInterface.

## API overview

The core element of the `Nonlinear` submodule is
[`Nonlinear.NonlinearData`](@ref):
```jldoctest nonlinear_developer
julia> import JuMP: Nonlinear

julia> data = Nonlinear.NonlinearData()
NonlinearData with available features:
  * :ExprGraph
```
[`Nonlinear.NonlinearData`](@ref) is a mutable struct that stores all of the
nonlinear information added to the model.

### Decision variables

Decision variables are represented by [`MOI.VariableIndex`](@ref)s. The user is
responsible for creating these.

### [Expressions](@id Nonlinear_Expressions)

The input data-structure is a Julia `Expr`. The input expressions can
incorporate [`MOI.VariableIndex`](@ref)es, but these must be interpolated into
the expression with `$`:
```jldoctest nonlinear_developer
julia> import JuMP: MOI

julia> x = MOI.VariableIndex(1)
MathOptInterface.VariableIndex(1)

julia> input = :(1 + sin($x)^2)
:(1 + sin(MathOptInterface.VariableIndex(1)) ^ 2)
```
There are a number of restrictions on the input `Expr`:
 * It cannot contain macros
 * It cannot contain broadcasting
 * It cannot contain splatting (except in limited situations)
 * It cannot contain linear algebra, such as matrix-vector products
 * It cannot contain generator expressions, including `sum(i for i in S)`

Given an input expression, add an expression using
[`Nonlinear.add_expression`](@ref):
```jldoctest nonlinear_developer
julia> expr = Nonlinear.add_expression(data, input)
JuMP.Nonlinear.ExpressionIndex(1)
```
The return value, `expr`, is a [`Nonlinear.ExpressionIndex`](@ref) that can
then be interpolated into other input expressions.

### [Parameters](@id Nonlinear_Parameters)

In addition to constant literals like `1` or `1.23`, you can create parameters.
Parameter are constants that you can change before passing the expression to the
solver. Create a parameter using [`Nonlinear.add_parameter`](@ref), which
accepts a default value:
```jldoctest nonlinear_developer
julia> p = Nonlinear.add_parameter(data, 1.23)
JuMP.Nonlinear.ParameterIndex(1)
```
The return value, `p`, is a [`Nonlinear.ParameterIndex`](@ref) that can then be
interpolated into other input expressions.

Update a parameter using [`Nonlinear.set_parameter`](@ref):
```jldoctest nonlinear_developer
julia> Nonlinear.set_parameter(data, p, 4.56)
```

### [Objectives](@id Nonlinear_Objectives)

Set a nonlinear objective using [`Nonlinear.set_objective`](@ref):
```jldoctest nonlinear_developer
julia> Nonlinear.set_objective(data, :($p + $expr + $x))
```

### [Constraints](@id Nonlinear_Constraints)

Add a constraint using [`Nonlinear.add_constraint`](@ref):
```jldoctest nonlinear_developer
julia> c = Nonlinear.add_constraint(data, :(1 + sqrt($x) <= 2.0))
JuMP.Nonlinear.ConstraintIndex(1)
```
The return value, `c`, is a [`Nonlinear.ConstraintIndex`](@ref) that is a unique
identifier for the constraint. Interval constraints are also supported:
```jldoctest nonlinear_developer
julia> c2 = Nonlinear.add_constraint(data, :(-1.0 <= 1 + sqrt($x) <= 2.0))
JuMP.Nonlinear.ConstraintIndex(2)
```

Delete a constraint using [`Nonlinear.delete`](@ref):
```jldoctest nonlinear_developer
julia> Nonlinear.delete(data, c2)
```

## User-defined operators

By default, `Nonlinear` supports a wide range of univariate and multivariate
operators. However, you can also define your own operators by _registering_
them.

### Univariate operators

Register a univariate user-defined operator using
[`Nonlinear.register_operator`](@ref):
```jldoctest nonlinear_developer
julia> f(x) = 1 + sin(x)^2
f (generic function with 1 method)

julia> Nonlinear.register_operator(data, :my_f, 1, f)
```
Now, you can use `:my_f` in expressions:
```jldoctest nonlinear_developer
julia> new_expr = Nonlinear.add_expression(data, :(my_f($x + 1)))
JuMP.Nonlinear.ExpressionIndex(2)
```
By default, `Nonlinear` will compute first- and second-derivatives of the
registered operator using `ForwardDiff.jl`. Override this by passing functions
which compute the respective derivative:
```jldoctest nonlinear_developer
julia> f′(x) = 2 * sin(x) * cos(x)
f′ (generic function with 1 method)

julia> Nonlinear.register_operator(data, :my_f2, 1, f, f′)
```
or
```jldoctest nonlinear_developer
julia> f′′(x) = 2 * (cos(x)^2 - sin(x)^2)
f′′ (generic function with 1 method)

julia> Nonlinear.register_operator(data, :my_f3, 1, f, f′, f′′)
```

### Multivariate operators

Register a multivariate user-defined operator using
[`Nonlinear.register_operator`](@ref):
```jldoctest nonlinear_developer
julia> g(x...) = x[1]^2 + x[1] * x[2] + x[2]^2
g (generic function with 1 method)

julia> Nonlinear.register_operator(data, :my_g, 2, g)
```
Now, you can use `:my_f` in expressions:
```jldoctest nonlinear_developer
julia> new_expr = Nonlinear.add_expression(data, :(my_g($x + 1, $x)))
JuMP.Nonlinear.ExpressionIndex(3)
```
By default, `Nonlinear` will compute the gradient of the registered
operator using `ForwardDiff.jl`. (Hessian information is not supported.)
Over-ride this by passing a function to compute the gradient:
```jldoctest nonlinear_developer
julia> function ∇g(ret, x...)
           ret[1] = 2 * x[1] + x[2]
           ret[2] = x[1] + 2 * x[2]
           return
       end
∇g (generic function with 1 method)

julia> Nonlinear.register_operator(data, :my_g2, 2, g, ∇g)
```

## [MathOptInterface](@id Nonlinear_MOI_interface)

`Nonlinear` implements the MathOptInterface API to allow solvers to query the
function and derivative information of our nonlinear model `data`. However,
before we can call [`MOI.initialize`](@ref), we need to set an
[`Nonlinear.AbstractAutomaticDifferentiation`](@ref).

There are two to choose from within JuMP, although other packages may add more
options by sub-typing [`Nonlinear.AbstractAutomaticDifferentiation`](@ref):
 * [`Nonlinear.Default`](@ref)

If we set [`Nonlinear.Default`](@ref), then we get access to `:ExprGraph`:
```jldoctest nonlinear_developer
julia> Nonlinear.set_differentiation_backend(data, Nonlinear.Default(), [x])

julia> data
NonlinearData with available features:
  * :ExprGraph
```
!!! note
    [`Nonlinear.set_differentiation_backend`](@ref) requires an ordered list of
    the variables that are included in the model. This order corresponds to the
    the order of the primal decision vector `x` which is passed to the various
    functions in MOI's nonlinear API.

The `:ExprGraph` feature means we can call [`MOI.objective_expr`](@ref) and
[`MOI.constraint_expr`](@ref) to retrieve the expression graph of the problem.
However, we cannot call gradient terms such as
[`MOI.eval_objective_gradient`](@ref) because [`Nonlinear.Default`](@ref) does
not know how to differentiate a nonlinear expression.

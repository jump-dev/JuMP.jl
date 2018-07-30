Interacting with solvers
========================

A JuMP model keeps a [MathOptInterface (MOI)](https://github.com/JuliaOpt/MathOptInterface.jl)
backend internally that stores the optimization problem and acts as the
optimization solver (the backend can also not support optimization, e.g. it can
simply store the model in a file). JuMP can be viewed as a lightweight
user-friendly layer on top of the MOI backend:

* JuMP does not maintain any copy of the model outside this MOI backend.
* JuMP variable (resp. constraint) references are simple structures containing
  both a reference to the JuMP model and the MOI index of the variable (resp.
  constraint).
* JuMP gives the constraints to the MOI backend in the form provided by the user
  without doing any automatic reformulation.
* variables additions, constraints additions/modifications and objective
  modifications are directly applied to the MOI backend thus expecting the
  backend to support such modifications.

While this allows JuMP API to to be a thin wrapper on top of the solver API,
as mentioned in the last point above, this seems rather demanding on the
solver. Indeed, while some solvers support incremental building of the model and
modifications before and after solve, other solvers only support the model being
copied at once before solve. Moreover it seems to require all solvers to
implement all possible reformulations independently which seems both very
ambitious and might generate a lot of duplicated code.

These apparent limitations are in fact addressed at the MOI level in a manner
that is completely transparent to JuMP. While the MOI API may seem very
demanding, it allows MOI models to be a succession of lightweight MOI layers
that fill the gap between JuMP requirements and the solver capabilities.

JuMP models can be created in three different modes: Automatic, Manual and
Direct.

## Automatic and Manual modes

In Automatic and Manual modes, two MOI layers are automatically applied to the
optimizer:

* `CachingOptimizer`: it maintain a cache of the model so that when the
  the optimizer does not support an incremental change to the model, the
  optimizer's internal model can be discarded and restored from the cache just
  before optimization. The `CachingOptimizer` has two different modes: Automatic
  and Manual corresponding to the two JuMP modes with the same names.
* `LazyBridgeOptimizer`: when a constraint added is not supported by the
  optimizer, it tries transform the constraint into an equivalent form,
  possibly adding new variables and constraints that are supported by the
  optimizer. The applied transformations are selected among known recipes
  which are called bridges. A few default bridges are defined in MOI but new
  ones can be defined and added to the `LazyBridgeOptimizer` used by JuMP.

See the [MOI documentation](http://www.juliaopt.org/MathOptInterface.jl/stable/)
for more details on these two MOI layers.

To create a fresh new JuMP model (or a fresh new copy of a JuMP model), JuMP
needs to create a new empty optimizer instance. New optimizer instances can
be obtained using a factory that can be created using the
[`with_optimizer`](@ref) function:
```@docs
with_optimizer
```

The factory can be set to the JuMP model using the [`JuMP.setoptimizer`](@ref)
function:
```@docs
JuMP.setoptimizer
```

New JuMP models are created using the [`Model`](@ref) constructor:
```@docs
Model()
Model(::JuMP.OptimizerFactory)
```

## Direct mode

JuMP models can be created in Direct mode using the [`JuMP.direct_model`](@ref)
function.
```@docs
JuMP.direct_model
```

TODO: How to set parameters (solver
specific and generic). Status codes. Accessing the result.
How to accurately measure the solve time.

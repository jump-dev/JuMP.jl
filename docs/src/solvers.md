Interacting with solvers
========================

A JuMP model keeps a [MathOptInterface (MOI)](https://github.com/jump-dev/MathOptInterface.jl)
*backend* of type `MOI.ModelLike` that stores the optimization
problem and acts as the optimization solver. We call it an MOI *backend* and not
optimizer as it can also be a wrapper around an optimization file format such as
MPS that writes the JuMP model in a file. From JuMP, the MOI
backend can be accessed using the [`backend`](@ref) function. JuMP can be
viewed as a lightweight, user-friendly layer on top of the MOI backend, in the
sense that:

* JuMP does not maintain any copy of the model outside this MOI backend.
* JuMP variable (resp. constraint) references are simple structures containing
  both a reference to the JuMP model and the MOI index of the variable (resp.
  constraint).
* JuMP gives the constraints to the MOI backend in the form provided by the user
  without doing any automatic reformulation.
* variables additions, constraints additions/modifications and objective
  modifications are directly applied to the MOI backend thus expecting the
  backend to support such modifications.

While this allows JuMP to be a thin wrapper on top of the solver API, as
mentioned in the last point above, this seems rather demanding on the solver.
Indeed, while some solvers support incremental building of the model and
modifications before and after solve, other solvers only support the model being
copied at once before solve. Moreover, it seems to require all solvers to
implement all possible reformulations independently which seems both very
ambitious and might generate a lot of duplicated code.

These apparent limitations are addressed at level of MOI in a manner
that is completely transparent to JuMP. While the MOI API may seem very
demanding, it allows MOI models to be a succession of lightweight MOI layers
that fill the gap between JuMP requirements and the solver capabilities. The
remainder of this section describes how JuMP interacts with the MOI backend.

JuMP models can be created in three different modes: `AUTOMATIC`, `MANUAL` and
`DIRECT`.

## Automatic and Manual modes

In `AUTOMATIC` and `MANUAL` modes, two MOI layers are automatically applied to
the optimizer:

* `CachingOptimizer`: maintains a cache of the model so that when the optimizer
  does not support an incremental change to the model, the optimizer's internal
  model can be discarded and restored from the cache just before optimization.
  The `CachingOptimizer` has two different modes: `AUTOMATIC` and `MANUAL`
  corresponding to the two JuMP modes with the same names.
* `LazyBridgeOptimizer` (this can be disabled using the `bridge_constraints`
  keyword argument to [`Model`](@ref) constructor): when a constraint added is
  not supported by the optimizer, it attempts to transform the constraint into
  an equivalent form, possibly adding new variables and constraints that are
  supported by the optimizer. The applied transformations are selected among
  known recipes which are called bridges. A few default bridges are defined in
  MOI but new ones can be defined and added to the `LazyBridgeOptimizer` used by
  JuMP.

See the [MOI documentation](http://jump.dev/MathOptInterface.jl/v0.9.1/)
for more details on these two MOI layers.

To attach an optimizer to a JuMP model, JuMP needs to be able to create a new
empty optimizer instance. For this reason, we provide JuMP with a function
that creates a new optimizer (i.e., an optimizer factory), instead of a concrete
optimizer object.

The factory can be provided either at model construction time by calling
[`set_optimizer`](@ref). An optimizer must be set before a call to
[`optimize!`](@ref). The optimizer can be grouped with attributes to be
set before optimization with [`optimizer_with_attributes`](@ref).
```@docs
set_optimizer
optimizer_with_attributes
NoOptimizer
JuMP.optimize!
```

New JuMP models are created using the [`Model`](@ref) constructor:
```@docs
Model()
Model(::Any)
```

A JuMP model may be reused by emptying it first:
```@docs
Base.empty!(::Model)
```

```@meta
# TODO: how to control the caching optimizer states
```

## Direct mode

JuMP models can be created in `DIRECT` mode using the
[`JuMP.direct_model`](@ref) function.
```@docs
JuMP.direct_model
```

```@docs
JuMP.backend
```

## Solver attributes

Some solver attributes can be queried and set through JuMP models.

```@docs
solver_name

bridge_constraints

get_optimizer_attribute
set_optimizer_attribute
set_optimizer_attributes
set_silent
unset_silent
set_time_limit_sec
unset_time_limit_sec
time_limit_sec
```

## File formats

JuMP can write models to a variety of file-formats using [`write_to_file`](@ref)
and [`Base.write`](@ref).

The file formats are defined [within MathOptInterface in the FileFormats enumeration](https://jump.dev/MathOptInterface.jl/v0.9/apireference/#MathOptInterface.FileFormats.FileFormat).

```@docs
write_to_file
Base.write(::IO, ::Model; ::MOI.FileFormats.FileFormat)
```

JuMP models can be created from file formats using [`read_from_file`](@ref) and
[`Base.read`](@ref).
```@docs
read_from_file
Base.read(::IO, ::Type{Model}; ::MOI.FileFormats.FileFormat)
```

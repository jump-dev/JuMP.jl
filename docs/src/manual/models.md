```@meta
CurrentModule = JuMP
DocTestSetup = quote
    using JuMP, GLPK
end
DocTestFilters = [r"≤|<=", r"≥|>=", r" == | = ", r" ∈ | in ", r"MathOptInterface|MOI"]
```

# Models

A JuMP model keeps a [MathOptInterface (MOI)](https://github.com/jump-dev/MathOptInterface.jl)
*backend* of type `MOI.ModelLike` that stores the optimization
problem and acts as the optimization solver. We call it an MOI *backend* and not
optimizer as it can also be a wrapper around an optimization file format such as
MPS that writes the JuMP model in a file. From JuMP, the MOI
backend can be accessed using the [`backend`](@ref) function.

JuMP can be viewed as a lightweight, user-friendly layer on top of the MOI
backend, in the sense that:

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

## Create a model

Create a model by passing an optimizer to [`Model`](@ref):
```jldoctest
julia> model = Model(GLPK.Optimizer)
A JuMP Model
Feasibility problem with:
Variables: 0
Model mode: AUTOMATIC
CachingOptimizer state: EMPTY_OPTIMIZER
Solver name: GLPK
```
or by calling [`set_optimizer`](@ref) on an empty [`Model`](@ref):
```jldoctest
julia> model = Model()
A JuMP Model
Feasibility problem with:
Variables: 0
Model mode: AUTOMATIC
CachingOptimizer state: NO_OPTIMIZER
Solver name: No optimizer attached.

julia> set_optimizer(model, GLPK.Optimizer)
```

!!! info
    JuMP uses "optimizer" as a synonym for "solver." Our convention is to use
    "solver" to refer to the underlying software, and use "optimizer" to refer
    to the Julia object that wraps the solver. For example, `GLPK` is a solver,
    and `GLPK.Optimizer` is an optimizer.

Use [`optimizer_with_attributes`](@ref) to create an optimizer with some
attributes initialized:
```jldoctest
julia> model = Model(optimizer_with_attributes(GLPK.Optimizer, "msg_lev" => 0))
A JuMP Model
Feasibility problem with:
Variables: 0
Model mode: AUTOMATIC
CachingOptimizer state: EMPTY_OPTIMIZER
Solver name: GLPK
```

Alternatively, you can create a function which takes no arguments and returns
an initialized `Optimizer` object:
```jldoctest
julia> function my_optimizer()
           model = GLPK.Optimizer()
           MOI.set(model, MOI.RawParameter("msg_lev"), 0)
           return model
       end
my_optimizer (generic function with 1 method)

julia> model = Model(my_optimizer)
A JuMP Model
Feasibility problem with:
Variables: 0
Model mode: AUTOMATIC
CachingOptimizer state: EMPTY_OPTIMIZER
Solver name: GLPK
```

## Direct mode

JuMP models can be created in `MOI.DIRECT` mode using [`direct_model`](@ref):
```jldoctest direct_mode
julia> model = direct_model(GLPK.Optimizer())
A JuMP Model
Feasibility problem with:
Variables: 0
Model mode: DIRECT
Solver name: GLPK
```

The benefit of using [`direct_model`](@ref) is that there are no extra layers
between `model` and the provided optimizer:
```jldoctest direct_mode
julia> backend(model)
A GLPK model

julia> typeof(backend(model))
GLPK.Optimizer
```

A downside of direct mode is that there is no bridging layer. Therefore, only
constraints which are natively supported by the solver are supported. For
example, `GLPK.jl` does not implement constraints of the form `l <= a' x <= u`.
```julia direct_mode
julia> @variable(model, x)
x

julia> @constraint(model, 1 <= x <= 2)
ERROR: Constraints of type MathOptInterface.ScalarAffineFunction{Float64}-in-MathOptInterface.Interval{Float64} are not supported by the solver.
Stacktrace:
 [1] error(::String) at ./error.jl:33
 [2] moi_add_constraint(::GLPK.Optimizer, ::MathOptInterface.ScalarAffineFunction{Float64}, ::MathOptInterface.Interval{Float64}) at /Users/oscar/.julia/dev/JuMP/src/constraints.jl:459
 [3] add_constraint(::Model, ::ScalarConstraint{GenericAffExpr{Float64,VariableRef},MathOptInterface.Interval{Float64}}, ::String) at /Users/oscar/.julia/dev/JuMP/src/constraints.jl:473
 [4] top-level scope at /Users/oscar/.julia/dev/JuMP/src/macros.jl:448
 [5] top-level scope at REPL[43]:1
```

!!! info
    If the model was created as `Model(GLPK.Optimizer)`, the constraint
    `1 <= x <= 2` would be reformulated into two constraints: `x >= 1` and
    `x <= 2`.

## Turn off output

Use [`set_silent`](@ref) and [`unset_silent`](@ref) to disable or enable
printing output from the solver.
```jldoctest
julia> model = Model(GLPK.Optimizer);

julia> set_silent(model)
true

julia> unset_silent(model)
false
```

## Set a time limit

Use [`set_time_limit_sec`](@ref), [`unset_time_limit_sec`](@ref), and
[`time_limit_sec`](@ref) to manage time limits.
```jldoctest
julia> model = Model(GLPK.Optimizer);

julia> set_time_limit_sec(model, 60.0)
60.0

julia> time_limit_sec(model)
60.0

julia> unset_time_limit_sec(model)

julia> time_limit_sec(model)
2.147483647e6
```

## Write a model to file

JuMP can write models to a variety of file-formats using [`write_to_file`](@ref)
and [`Base.write`](@ref).

```jldoctest file_formats; setup=:(model = Model(); io = IOBuffer())
julia> write_to_file(model, "model.mps")

julia> write(io, model; format = MOI.FileFormats.FORMAT_MPS)
```

!!! info
    The supported file formats are defined by the [MOI.FileFormats.FileFormat](https://jump.dev/MathOptInterface.jl/v0.9/apireference/#MathOptInterface.FileFormats.FileFormat)
    enum.
    ```jldoctest
    julia> MOI.FileFormats.FileFormat
    Enum MathOptInterface.FileFormats.FileFormat:
    FORMAT_AUTOMATIC = 0
    FORMAT_CBF = 1
    FORMAT_LP = 2
    FORMAT_MOF = 3
    FORMAT_MPS = 4
    FORMAT_SDPA = 5
    ```

## Read a model from file

JuMP models can be created from file formats using [`read_from_file`](@ref) and
[`Base.read`](@ref).

```jldoctest file_formats
julia> model = read_from_file("model.mps")
A JuMP Model
Minimization problem with:
Variables: 0
Objective function type: GenericAffExpr{Float64,VariableRef}
Model mode: AUTOMATIC
CachingOptimizer state: NO_OPTIMIZER
Solver name: No optimizer attached.

julia> seekstart(io);

julia> model2 = read(io, Model; format = MOI.FileFormats.FORMAT_MPS)
A JuMP Model
Minimization problem with:
Variables: 0
Objective function type: GenericAffExpr{Float64,VariableRef}
Model mode: AUTOMATIC
CachingOptimizer state: NO_OPTIMIZER
Solver name: No optimizer attached.
```

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

See the [MOI documentation](https://jump.dev/MathOptInterface.jl/v0.9.1/)
for more details on these two MOI layers.

```@meta
# TODO: how to control the caching optimizer states
```

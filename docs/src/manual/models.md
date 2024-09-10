```@meta
CurrentModule = JuMP
DocTestSetup = quote
    using JuMP, HiGHS, SCS
end
DocTestFilters = [r"≤|<=", r"≥|>=", r" == | = ", r" ∈ | in ", r"MathOptInterface|MOI"]
```

# [Models](@id jump_models)

JuMP models are the fundamental building block that we use to construct
optimization problems. They hold things like the variables and constraints, as
well as which solver to use and even solution information.

!!! info
    JuMP uses "optimizer" as a synonym for "solver." Our convention is to use
    "solver" to refer to the underlying software, and use "optimizer" to refer
    to the Julia object that wraps the solver. For example, `HiGHS` is a solver,
    and `HiGHS.Optimizer` is an optimizer.

!!! tip
    See [Supported solvers](@ref) for a list of available solvers.

## Create a model

Create a model by passing an optimizer to [`Model`](@ref):
```jldoctest
julia> model = Model(HiGHS.Optimizer)
A JuMP Model
├ solver: HiGHS
├ objective_sense: FEASIBILITY_SENSE
├ num_variables: 0
├ num_constraints: 0
└ Names registered in the model: none
```

If you don't know which optimizer you will be using at creation time, create a
model without an optimizer, and then call [`set_optimizer`](@ref) at any time
prior to [`optimize!`](@ref):
```jldoctest
julia> model = Model()
A JuMP Model
├ solver: none
├ objective_sense: FEASIBILITY_SENSE
├ num_variables: 0
├ num_constraints: 0
└ Names registered in the model: none

julia> set_optimizer(model, HiGHS.Optimizer)
```

!!! tip
    Don't know what the fields `Model mode` and `CachingOptimizer state` mean?
    Read the [Backends](@ref) section.

### What is the difference?

For most models, there is no difference between passing the optimizer to
[`Model`](@ref), and calling [`set_optimizer`](@ref).

However, if an optimizer does not support a constraint in the model, the timing
of when an error will be thrown can differ:

 * If you pass an optimizer, an error will be thrown when you try to add the
   constraint.
 * If you call [`set_optimizer`](@ref), an error will be thrown when you try to
   solve the model via [`optimize!`](@ref).

Therefore, most users should pass an optimizer to [`Model`](@ref) because it
provides the earliest warning that your solver is not suitable for the model you
are trying to build. However, if you are modifying a problem by adding and
deleting different constraint types, you may need to use
[`set_optimizer`](@ref). See [Switching optimizer for the relaxed problem](@ref)
for an example of when this is useful.

### Reducing time-to-first-solve latency

By default, JuMP uses [bridges](@ref LazyBridgeOptimizer) to reformulate the
model you are building into an equivalent model supported by the solver.

However, if your model is already supported by the solver, bridges add latency
(read [The "time-to-first-solve" issue](@ref)). This is particularly noticeable
for small models.

To reduce the "time-to-first-solve,s" try passing `add_bridges = false`.
```jldoctest
julia> model = Model(HiGHS.Optimizer; add_bridges = false);
```
or
```jldoctest
julia> model = Model();

julia> set_optimizer(model, HiGHS.Optimizer; add_bridges = false)
```

However, be wary. If your model and solver combination needs bridges, an error
will be thrown:
```jldoctest
julia> model = Model(SCS.Optimizer; add_bridges = false);


julia> @variable(model, x)
x

julia> @constraint(model, 2x <= 1)
ERROR: Constraints of type MathOptInterface.ScalarAffineFunction{Float64}-in-MathOptInterface.LessThan{Float64} are not supported by the solver.

If you expected the solver to support your problem, you may have an error in your formulation. Otherwise, consider using a different solver.

The list of available solvers, along with the problem types they support, is available at https://jump.dev/JuMP.jl/stable/installation/#Supported-solvers.
[...]
```

### Solvers which expect environments

Some solvers accept (or require) positional arguments such as a license
environment or a path to a binary executable. For these solvers, you can pass
a function to [`Model`](@ref) which takes zero arguments and returns an instance
of the optimizer.

A common use-case for this is passing an environment or sub-solver to the
optimizer:
```jldoctest
julia> import HiGHS

julia> import MultiObjectiveAlgorithms as MOA

julia> model = Model(() -> MOA.Optimizer(HiGHS.Optimizer))
A JuMP Model
├ solver: MOA[algorithm=MultiObjectiveAlgorithms.Lexicographic, optimizer=HiGHS]
├ objective_sense: FEASIBILITY_SENSE
├ num_variables: 0
├ num_constraints: 0
└ Names registered in the model: none
```

## [Solver options](@id solver_options)

JuMP uses "attribute" as a synonym for "option." Use
[`optimizer_with_attributes`](@ref) to create an optimizer with some attributes
initialized:
```jldoctest
julia> model = Model(
           optimizer_with_attributes(HiGHS.Optimizer, "output_flag" => false),
       )
A JuMP Model
├ solver: HiGHS
├ objective_sense: FEASIBILITY_SENSE
├ num_variables: 0
├ num_constraints: 0
└ Names registered in the model: none
```

Alternatively, use [`set_attribute`](@ref) to set an attribute after
the model has been created:
```jldoctest
julia> model = Model(HiGHS.Optimizer);

julia> set_attribute(model, "output_flag", false)

julia> get_attribute(model, "output_flag")
false
```

You can also modify attributes within an [`optimizer_with_attributes`](@ref)
object:
```jldoctest
julia> solver = optimizer_with_attributes(HiGHS.Optimizer, "output_flag" => true);

julia> get_attribute(solver, "output_flag")
true

julia> set_attribute(solver, "output_flag", false)

julia> get_attribute(solver, "output_flag")
false

julia> model = Model(solver);
```

## Changing the number types

By default, the coefficients of affine and quadratic expressions are numbers
of type either `Float64` or `Complex{Float64}` (see [Complex number support](@ref)).

The type `Float64` can be changed using the [`GenericModel`](@ref) constructor:

```jldoctest
julia> model = GenericModel{Rational{BigInt}}();

julia> @variable(model, x)
x

julia> @expression(model, expr, 1 // 3 * x)
1//3 x

julia> typeof(expr)
GenericAffExpr{Rational{BigInt}, GenericVariableRef{Rational{BigInt}}}
```

Using a `value_type` other than `Float64` is an advanced operation and should be
used only if the underlying solver actually solves the problem using the
provided value type.

!!! warning
    [Nonlinear Modeling](@ref) is currently restricted to the `Float64` number
    type.

## Print the model

By default, `show(model)` will print a summary of the problem:
```jldoctest model_print
julia> model = Model(); @variable(model, x >= 0); @objective(model, Max, x);

julia> model
A JuMP Model
├ solver: none
├ objective_sense: MAX_SENSE
│ └ objective_function_type: VariableRef
├ num_variables: 1
├ num_constraints: 1
│ └ VariableRef in MOI.GreaterThan{Float64}: 1
└ Names registered in the model
  └ :x
```

Use `print` to print the formulation of the model (in IJulia, this will render
as LaTeX.
```jldoctest model_print
julia> print(model)
Max x
Subject to
 x ≥ 0
```

!!! warning
    This format is specific to JuMP and may change in any future release. It is
    not intended to be an instance format. To write the model to a file, use
    [`write_to_file`](@ref) instead.


Use [`latex_formulation`](@ref) to display the model in LaTeX form.
```jldoctest model_print
julia> latex_formulation(model)
$$ \begin{aligned}
\max\quad & x\\
\text{Subject to} \quad & x \geq 0\\
\end{aligned} $$
```

In IJulia (and Documenter), ending a cell in with [`latex_formulation`](@ref)
will render the model in LaTeX:

```@example
using JuMP                # hide
model = Model()           # hide
@variable(model, x >= 0)  # hide
@objective(model, Max, x) # hide
latex_formulation(model)
```

## Turn off output

Use [`set_silent`](@ref) and [`unset_silent`](@ref) to disable or enable
printing output from the solver.
```jldoctest
julia> model = Model(HiGHS.Optimizer);

julia> set_silent(model)

julia> unset_silent(model)
```

!!! tip
    Most solvers will also have a [solver-specific option](@ref solver_options)
    to provide finer-grained control over the output. Consult their README's for
    details.

## Set a time limit

Use [`set_time_limit_sec`](@ref), [`unset_time_limit_sec`](@ref), and
[`time_limit_sec`](@ref) to manage time limits.
```jldoctest time_limit
julia> model = Model(HiGHS.Optimizer);

julia> set_time_limit_sec(model, 60.0)


julia> time_limit_sec(model)
60.0

julia> unset_time_limit_sec(model)

julia> limit = time_limit_sec(model)

julia> limit === nothing
true
```

If your time limit is encoded as a `Dates.Period` object, use the following code
to convert it to `Float64` for [`set_time_limit_sec`](@ref):
```jldoctest time_limit
julia> import Dates

julia> seconds(x::Dates.Period) = 1e-3 * Dates.value(round(x, Dates.Millisecond))
seconds (generic function with 1 method)

julia> set_time_limit_sec(model, seconds(Dates.Hour(1)))

julia> time_limit_sec(model)
3600.0
```

!!! info
    Some solvers do not support time limits. In these cases, an error will be
    thrown.

## Write a model to file

JuMP can write models to a variety of file-formats using [`write_to_file`](@ref)
and [`Base.write`](@ref).

For most common file formats, the file type will be detected from the extension.

For example, here is how to write an [MPS file](https://en.wikipedia.org/wiki/MPS_(format)):
```jldoctest file_formats
julia> model = Model();

julia> write_to_file(model, "model.mps")
```

Other supported file formats include:

 * `.cbf` for the [Conic Benchmark Format](https://docs.mosek.com/latest/capi/cbf-format.html)
 * `.lp` for the [LP file format](https://docs.mosek.com/latest/capi/lp-format.html)
 * `.mof.json` for the [MathOptFormat](https://jump.dev/MathOptFormat/)
 * `.nl` for [AMPL's NL file format](https://en.wikipedia.org/wiki/Nl_(format))
 * `.rew` for the [REW file format](https://www.gurobi.com/documentation/9.5/refman/rew_format.html)
 * `.sdpa` and ".dat-s" for the [SDPA file format](http://plato.asu.edu/ftp/sdpa_format.txt)

To write to a specific `io::IO`, use [`Base.write`](@ref). Specify the file type
by passing a [`MOI.FileFormats.FileFormat`](@ref) enum.
```jldoctest file_formats
julia> model = Model();

julia> io = IOBuffer();

julia> write(io, model; format = MOI.FileFormats.FORMAT_MPS)
```

## Read a model from file

JuMP models can be created from file formats using [`read_from_file`](@ref) and
`Base.read`.

```jldoctest file_formats
julia> model = read_from_file("model.mps")
A JuMP Model
├ solver: none
├ objective_sense: MIN_SENSE
│ └ objective_function_type: AffExpr
├ num_variables: 0
├ num_constraints: 0
└ Names registered in the model: none

julia> seekstart(io);

julia> model2 = read(io, Model; format = MOI.FileFormats.FORMAT_MPS)
A JuMP Model
├ solver: none
├ objective_sense: MIN_SENSE
│ └ objective_function_type: AffExpr
├ num_variables: 0
├ num_constraints: 0
└ Names registered in the model: none
```

!!! note
    Because file formats do not serialize the containers of JuMP variables and
    constraints, the names in the model will _not_ be registered. Therefore, you
    cannot access named variables and constraints via `model[:x]`. Instead, use
    [`variable_by_name`](@ref) or [`constraint_by_name`](@ref) to access
    specific variables or constraints.

## Relax integrality

Use [`relax_integrality`](@ref) to remove any integrality constraints from the
model, such as integer and binary restrictions on variables.
[`relax_integrality`](@ref) returns a function that can be later called with
zero arguments to re-add the removed constraints:
```jldoctest
julia> model = Model();

julia> @variable(model, x, Int)
x

julia> num_constraints(model, VariableRef, MOI.Integer)
1

julia> undo = relax_integrality(model);

julia> num_constraints(model, VariableRef, MOI.Integer)
0

julia> undo()

julia> num_constraints(model, VariableRef, MOI.Integer)
1
```

### Switching optimizer for the relaxed problem

A common reason for relaxing integrality is to compute dual variables of the
relaxed problem. However, some mixed-integer linear solvers (for example, Cbc) do not
return dual solutions, even if the problem does not have integrality
restrictions.

Therefore, after [`relax_integrality`](@ref) you should call
[`set_optimizer`](@ref) with a solver that does support dual solutions, such as
Clp.

For example, instead of:
```julia
using JuMP, Cbc
model = Model(Cbc.Optimizer)
@variable(model, x, Int)
undo = relax_integrality(model)
optimize!(model)
reduced_cost(x)  # Errors
```
do:
```julia
using JuMP, Cbc, Clp
model = Model(Cbc.Optimizer)
@variable(model, x, Int)
undo = relax_integrality(model)
set_optimizer(model, Clp.Optimizer)
optimize!(model)
reduced_cost(x)  # Works
```

## Get the matrix representation

Use [`lp_matrix_data`](@ref) to return a data structure that represents the
matrix form of a linear program.

```jldoctest
julia> begin
           model = Model()
           @variable(model, x >= 1, Bin)
           @variable(model, 2 <= y)
           @variable(model, 3 <= z <= 4, Int)
           @constraint(model, x == 5)
           @constraint(model, 2x + 3y <= 6)
           @constraint(model, -4y >= 5z + 7)
           @constraint(model, -1 <= x + y <= 2)
           @objective(model, Max, 1 + 2x)
       end;

julia> data = lp_matrix_data(model);

julia> data.A
4×3 SparseArrays.SparseMatrixCSC{Float64, Int64} with 7 stored entries:
 1.0    ⋅     ⋅
  ⋅   -4.0  -5.0
 2.0   3.0    ⋅
 1.0   1.0    ⋅

julia> data.b_lower
4-element Vector{Float64}:
   5.0
   7.0
 -Inf
  -1.0

julia> data.b_upper
4-element Vector{Float64}:
  5.0
 Inf
  6.0
  2.0

julia> data.x_lower
3-element Vector{Float64}:
 1.0
 2.0
 3.0

julia> data.x_upper
3-element Vector{Float64}:
 Inf
 Inf
  4.0

julia> data.c
3-element Vector{Float64}:
 2.0
 0.0
 0.0

julia> data.c_offset
1.0

julia> data.sense
MAX_SENSE::OptimizationSense = 1

julia> data.integers
1-element Vector{Int64}:
 3

julia> data.binaries
1-element Vector{Int64}:
 1
```

!!! warning
    [`lp_matrix_data`](@ref) is intentionally limited in the types of problems
    that it supports and the structure of the matrices it outputs. It is mainly
    intended as a pedagogical and debugging tool. It should not be used to
    interface solvers, see [Implementing a solver interface](@ref) instead.

## Backends

!!! info
    This section discusses advanced features of JuMP. For new users, you may
    want to skip this section. You don't need to know how JuMP manages problems
    behind the scenes to create and solve JuMP models.

A JuMP [`Model`](@ref) is a thin layer around a *backend* of type
[`MOI.ModelLike`](@ref) that stores the optimization problem and acts as the
optimization solver.

However, if you construct a model like `Model(HiGHS.Optimizer)`, the backend is
not a `HiGHS.Optimizer`, but a more complicated object.

From JuMP, the MOI backend can be accessed using the [`backend`](@ref) function.
Let's see what the [`backend`](@ref) of a JuMP [`Model`](@ref) is:
```jldoctest models_backends
julia> model = Model(HiGHS.Optimizer);

julia> b = backend(model)
MOIU.CachingOptimizer
├ state: EMPTY_OPTIMIZER
├ mode: AUTOMATIC
├ model_cache: MOIU.UniversalFallback{MOIU.Model{Float64}}
│ ├ ObjectiveSense: FEASIBILITY_SENSE
│ ├ ObjectiveFunctionType: MOI.ScalarAffineFunction{Float64}
│ ├ NumberOfVariables: 0
│ └ NumberOfConstraints: 0
└ optimizer: MOIB.LazyBridgeOptimizer{HiGHS.Optimizer}
  ├ Variable bridges: none
  ├ Constraint bridges: none
  ├ Objective bridges: none
  └ model: A HiGHS model with 0 columns and 0 rows.
```

Uh oh. Even though we passed a `HiGHS.Optimizer`, the backend is a much more
complicated object.

### CachingOptimizer

A `MOIU.CachingOptimizer` is a layer that abstracts the difference between
solvers that support incremental modification (for example, they support adding
variables one-by-one), and solvers that require the entire problem in a single
API call (for example, they only accept the `A`, `b` and `c` matrices of a linear
program).

It has two parts:

 1. A cache, where the model can be built and modified incrementally
    ```jldoctest models_backends
    julia> b.model_cache
    MOIU.UniversalFallback{MOIU.Model{Float64}}
    ├ ObjectiveSense: FEASIBILITY_SENSE
    ├ ObjectiveFunctionType: MOI.ScalarAffineFunction{Float64}
    ├ NumberOfVariables: 0
    └ NumberOfConstraints: 0
    ```
 2. An optimizer, which is used to solve the problem
    ```jldoctest models_backends
    julia> b.optimizer
    MOIB.LazyBridgeOptimizer{HiGHS.Optimizer}
    ├ Variable bridges: none
    ├ Constraint bridges: none
    ├ Objective bridges: none
    └ model: A HiGHS model with 0 columns and 0 rows.
    ```

!!! info
    The [LazyBridgeOptimizer](@ref) section explains what a
    `LazyBridgeOptimizer` is.

The `CachingOptimizer` has logic to decide when to copy the problem from the
cache to the optimizer, and when it can efficiently update the optimizer
in-place.

A `CachingOptimizer` may be in one of three possible states:

* `NO_OPTIMIZER`: The CachingOptimizer does not have any optimizer.
* `EMPTY_OPTIMIZER`: The CachingOptimizer has an empty optimizer, and it is not
  synchronized with the cached model.
* `ATTACHED_OPTIMIZER`: The CachingOptimizer has an optimizer, and it is
  synchronized with the cached model.

A `CachingOptimizer` has two modes of operation:

* `AUTOMATIC`: The `CachingOptimizer` changes its state when necessary. For
  example, [`optimize!`](@ref) will automatically call `attach_optimizer` (an
  optimizer must have been previously set). Attempting to add a constraint or
  perform a modification not supported by the optimizer results in a drop to
  `EMPTY_OPTIMIZER` mode.
* `MANUAL`: The user must change the state of the `CachingOptimizer` using
  [`MOIU.reset_optimizer(::JuMP.Model)`](@ref),
  [`MOIU.drop_optimizer(::JuMP.Model)`](@ref), and
  [`MOIU.attach_optimizer(::JuMP.Model)`](@ref). Attempting to perform
  an operation in the incorrect state results in an error.

By default [`Model`](@ref) will create a `CachingOptimizer` in `AUTOMATIC` mode.

### LazyBridgeOptimizer

The second layer that JuMP applies automatically is a [`MOI.Bridges.LazyBridgeOptimizer`](@ref).
A [`MOI.Bridges.LazyBridgeOptimizer`](@ref) is an MOI layer that attempts to
transform the problem from the formulation provided by the user into an
equivalent problem supported by the solver. This may involve adding new
variables and constraints to the optimizer. The transformations are selected
from a set of known recipes called _bridges_.

A common example of a bridge is one that splits an interval constraint like
`@constraint(model, 1 <= x + y <= 2)` into two constraints,
`@constraint(model, x + y >= 1)` and `@constraint(model, x + y <= 2)`.

Use the `add_bridges = false` keyword to remove the bridging layer:
```jldoctest
julia> model = Model(HiGHS.Optimizer; add_bridges = false)
A JuMP Model
├ solver: HiGHS
├ objective_sense: FEASIBILITY_SENSE
├ num_variables: 0
├ num_constraints: 0
└ Names registered in the model: none

julia> backend(model)
MOIU.CachingOptimizer
├ state: EMPTY_OPTIMIZER
├ mode: AUTOMATIC
├ model_cache: MOIU.UniversalFallback{MOIU.Model{Float64}}
│ ├ ObjectiveSense: FEASIBILITY_SENSE
│ ├ ObjectiveFunctionType: MOI.ScalarAffineFunction{Float64}
│ ├ NumberOfVariables: 0
│ └ NumberOfConstraints: 0
└ optimizer: A HiGHS model with 0 columns and 0 rows.
```

Bridges can be added and removed from a [`MOI.Bridges.LazyBridgeOptimizer`](@ref)
using [`add_bridge`](@ref) and [`remove_bridge`](@ref). Use
[`print_active_bridges`](@ref) to see which bridges are used to reformulate the
model. Read the [Example: ellipsoid approximation](@ref) tutorial for more details.

### Unsafe backend

In some advanced use-cases, it is necessary to work with the inner optimization
model directly. To access this model, use [`unsafe_backend`](@ref):
```jldoctest models_backends
julia> backend(model)
MOIU.CachingOptimizer
├ state: EMPTY_OPTIMIZER
├ mode: AUTOMATIC
├ model_cache: MOIU.UniversalFallback{MOIU.Model{Float64}}
│ ├ ObjectiveSense: FEASIBILITY_SENSE
│ ├ ObjectiveFunctionType: MOI.ScalarAffineFunction{Float64}
│ ├ NumberOfVariables: 0
│ └ NumberOfConstraints: 0
└ optimizer: MOIB.LazyBridgeOptimizer{HiGHS.Optimizer}
  ├ Variable bridges: none
  ├ Constraint bridges: none
  ├ Objective bridges: none
  └ model: A HiGHS model with 0 columns and 0 rows.

julia> unsafe_backend(model)
A HiGHS model with 0 columns and 0 rows.
```

!!! warning
    [`backend`](@ref) and [`unsafe_backend`](@ref) are advanced routines. Read
    their docstrings to understand the caveats of their usage, and only call
    them if you wish to access low-level solver-specific functions.

## Direct mode

Using a `CachingOptimizer` results in an additional copy of the model being
stored by JuMP in the `.model_cache` field. To avoid this overhead, create a
JuMP model using [`direct_model`](@ref):
```jldoctest direct_mode
julia> model = direct_model(HiGHS.Optimizer())
A JuMP Model
├ mode: DIRECT
├ solver: HiGHS
├ objective_sense: FEASIBILITY_SENSE
├ num_variables: 0
├ num_constraints: 0
└ Names registered in the model: none
```

!!! warning
    Solvers that do not support incremental modification do not support
    `direct_model`. An error will be thrown, telling you to use a
    `CachingOptimizer` instead.

The benefit of using [`direct_model`](@ref) is that there are no extra layers
(for example, `Cachingoptimizer` or `LazyBridgeOptimizer`) between `model` and the
provided optimizer:
```jldoctest direct_mode
julia> backend(model)
A HiGHS model with 0 columns and 0 rows.
```

A downside of direct mode is that there is no bridging layer. Therefore, only
constraints which are natively supported by the solver are supported. For
example, `HiGHS.jl` does not implement quadratic constraints:
```jldoctest direct_mode
julia> model = direct_model(HiGHS.Optimizer());

julia> set_silent(model)

julia> @variable(model, x[1:2]);

julia> @constraint(model, x[1]^2 + x[2]^2 <= 2)
ERROR: Constraints of type MathOptInterface.ScalarQuadraticFunction{Float64}-in-MathOptInterface.LessThan{Float64} are not supported by the solver.

If you expected the solver to support your problem, you may have an error in your formulation. Otherwise, consider using a different solver.

The list of available solvers, along with the problem types they support, is available at https://jump.dev/JuMP.jl/stable/installation/#Supported-solvers.
Stacktrace:
```

!!! warning
    Another downside of direct mode is that the behavior of querying solution
    information after modifying the problem is solver-specific. This can lead to
    errors, or the solver silently returning an incorrect value. See
    [OptimizeNotCalled errors](@ref) for more information.

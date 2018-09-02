Quick Start Guide
=================

This quick start guide will introduce the main concepts of JuMP. If you are
familiar with another modeling language embedded in a high-level language such
as PuLP (Python) or a solver-specific interface you will find most of this
familiar. If you are coming from an AMPL or similar background, you may find
some of the concepts novel but the general appearance will still be familiar.

The example in this guide is deliberately kept simple. There are more complex
examples in the [`JuMP/examples/` folder](https://github.com/JuliaOpt/JuMP.jl/tree/master/examples).

Once JuMP is installed, to use JuMP in your programs, you just need to say:
```jldoctest quickstart_example
julia> using JuMP
```

You also need to include a Julia package which provides an appropriate solver.
One such solver is `GLPK.Optimizer`, which is provided by the
[GLPK.jl package](https://github.com/JuliaOpt/GLPK.jl).
```julia
julia> using GLPK
```
See [Installation Guide](@ref) for a list of other solvers you can use.

Models are created with the `Model()` function. The `with_optimizer` syntax is
used to specify the optimizer to be used:
```julia
julia> model = Model(with_optimizer(GLPK.Optimizer))
A JuMP Model
```

```@meta
DocTestSetup = quote
    # Using a caching optimizer removes the need to # load a solver such as GLPK
    # for building the documentation.
    const MOI = JuMP.MathOptInterface
    model = Model(with_optimizer(MOI.Utilities.MockOptimizer,
                                 JuMP.JuMPMOIModel{Float64}(),
                                 evalobjective=false))
end
```
!!! note
    Your model doesn't have to be called `model` - it's just a name.

There are a few options for defining a variable, depending on whether you want
to have lower bounds, upper bounds, both bounds, or even no bounds. The
following commands will create two variables, `x` and `y`, with both lower and
upper bounds. Note the first argument is our model variable ``model``. These
variables are associated with this model and cannot be used in another model.
```jldoctest quickstart_example
julia> @variable(model, 0 <= x <= 2)
x

julia> @variable(model, 0 <= y <= 30)
y
```
See the [Variables](@ref) section for more information on creating variables.

```@meta
DocTestSetup = nothing
```

Next we'll set our objective. Note again the `model`, so we know which model's
objective we are setting! The objective sense, `Max` or `Min`, should be
provided as the second argument. Note also that we don't have a multiplication
`*` symbol between `5` and our variable `x` - Julia is smart enough to not need
it! Feel free to stick with `*` if it makes you feel more comfortable, as we
have done with `3 * y`:
```jldoctest quickstart_example
julia> @objective(model, Max, 5x + 3 * y)
```

Adding constraints is a lot like setting the objective. Here we create a
less-than-or-equal-to constraint using `<=`, but we can also create equality
constraints using `==` and greater-than-or-equal-to constraints with `>=`:
```jldoctest quickstart_example; filter=r"≤|<="
julia> con = @constraint(model, 1x + 5y <= 3)
x + 5 y <= 3.0
```
Note that we bind the constraint to the Julia variable `con` for later analysis.

Models are solved with the `JuMP.optimize!` function:
```jldoctest quickstart_example
julia> JuMP.optimize!(model)
```

```@meta
DocTestSetup = quote
    # Now we load in the solution. Using a caching optimizer removes the need to
    # load a solver such as GLPK for building the documentation.
    mock = JuMP.caching_optimizer(model).optimizer
    MOI.set!(mock, MOI.TerminationStatus(), MOI.Success)
    MOI.set!(mock, MOI.PrimalStatus(), MOI.FeasiblePoint)
    MOI.set!(mock, MOI.DualStatus(), MOI.FeasiblePoint)
    MOI.set!(mock, MOI.ResultCount(), 1)
    MOI.set!(mock, MOI.ObjectiveValue(), 10.6)
    MOI.set!(mock, MOI.VariablePrimal(), JuMP.optimizer_index(x), 2.0)
    MOI.set!(mock, MOI.VariablePrimal(), JuMP.optimizer_index(y), 0.2)
    MOI.set!(mock, MOI.ConstraintDual(), JuMP.optimizer_index(con), -0.6)
    MOI.set!(mock, MOI.ConstraintDual(), JuMP.optimizer_index(JuMP.UpperBoundRef(x)), -4.4)
    MOI.set!(mock, MOI.ConstraintDual(), JuMP.optimizer_index(JuMP.LowerBoundRef(y)), 0.0)
end
```

After the call to `JuMP.optimize!` has finished, we need to understand why the
optimizer stopped. This can be for a number of reasons. First, the solver might
have found the optimal solution, or proved that the problem is infeasible.
However, it might also have run into numerical difficulties, or terminated due
to a setting such as a time limit. We can ask the solver why it stopped using
the `JuMP.termination_status` function:
```jldoctest quickstart_example
julia> JuMP.termination_status(model)
Success::TerminationStatusCode = 0
```
In this case, `GLPK` returned `Success`. This does not mean that it has found
the optimal solution. Instead, it indicates that GLPK has finished running and
did not encounter any errors or user-provided termination limits.

```@meta
DocTestSetup = nothing
```

To understand the reason for termination in more detail, we need to query
`JuMP.primalstatus`:
```jldoctest quickstart_example
julia> JuMP.primal_status(model)
FeasiblePoint::ResultStatusCode = 0
```
This indicates that GLPK has found a `FeasiblePoint` to the primal problem.
Coupled with the `Success` from `JuMP.termination_status`, we can infer that GLPK
has indeed found the optimal solution. We can also query `JuMP.dual_status`:
```jldoctest quickstart_example
julia> JuMP.dual_status(model)
FeasiblePoint::ResultStatusCode = 0
```
Like the `primal_status`, GLPK indicates that it has found a `FeasiblePoint` to
the dual problem.

Finally, we can query the result of the optimization. First, we can query the
objective value:
```jldoctest quickstart_example
julia> JuMP.objective_value(model)
10.6
```
We can also query the primal result values of the `x` and `y` variables:
```jldoctest quickstart_example
julia> JuMP.result_value(x)
2.0

julia> JuMP.result_value(y)
0.2
```

We can also query the value of the dual variable associated with the constraint
`con` (which we bound to a Julia variable when defining the constraint):
```jldoctest quickstart_example
julia> JuMP.result_dual(con)
-0.6
```

To query the dual variables associated with the variable bounds, things are a
little trickier as we first need to obtain a reference to the constraint:
```jldoctest quickstart_example; filter=r"≤|<="
julia> x_upper = JuMP.UpperBoundRef(x)
x <= 2.0

julia> JuMP.result_dual(x_upper)
-4.4
```
A similar process can be followed to obtain the dual of the lower bound
constraint on `y`:
```jldoctest quickstart_example; filter=r"≥|>="
julia> y_lower = JuMP.LowerBoundRef(y)
y >= 0.0

julia> JuMP.result_dual(y_lower)
0.0
```

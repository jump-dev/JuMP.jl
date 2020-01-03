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

Models are created with the [`Model`](@ref) function. The optimizer can be set
either in `Model()` or by calling [`set_optimizer`](@ref):
```julia
julia> model = Model(GLPK.Optimizer)
A JuMP Model
Feasibility problem with:
Variables: 0
Model mode: AUTOMATIC
CachingOptimizer state: NO_OPTIMIZER
Solver name: GLPK
```
equivalently,
```julia
julia> model = Model();
julia> set_optimizer(model, GLPK.Optimizer);
julia> model
A JuMP Model
Feasibility problem with:
Variables: 0
Model mode: AUTOMATIC
CachingOptimizer state: NO_OPTIMIZER
Solver name: GLPK
```

!!! note
    The term "solver" is used as a synonym for "optimizer". The convention in
    code, however, is to always use "optimizer", e.g., `GLPK.Optimizer`.

```@meta
DocTestSetup = quote
    # Using a mock optimizer removes the need to load a solver such as GLPK for
    # building the documentation.
    const MOI = JuMP.MathOptInterface
    model = Model(() -> MOI.Utilities.MockOptimizer(
                            MOIU.Model{Float64}(),
                            eval_objective_value = false,
                            eval_variable_constraint_dual = false))
end
```
!!! note
    Your model doesn't have to be called `model` - it's just a name.

The following commands will create two variables, `x` and `y`, with both lower
and upper bounds. Note the first argument is our model `model`. These variables
(`x` and `y`) are associated with this model and cannot be used in another
model.
```jldoctest quickstart_example
julia> @variable(model, 0 <= x <= 2)
x

julia> @variable(model, 0 <= y <= 30)
y
```
See the [Variables](@ref) section for more information on creating variables,
including the syntax for specifying different combinations of bounds, i.e.,
only lower bounds, only upper bounds, or no bounds.

```@meta
DocTestSetup = nothing
```

Next we'll set our objective. Note again the `model`, so we know which model's
objective we are setting! The objective sense, `Max` or `Min`, should be
provided as the second argument. Note also that we don't have a multiplication
`*` symbol between `5` and our variable `x` - Julia is smart enough to not need
it! Feel free to use `*` if it makes you feel more comfortable, as we have done
with `3 * y`. (We have been intentionally inconsistent here to demonstrate
different syntax; however, it is good practice to pick one way or the other
consistently in your code.)
```jldoctest quickstart_example
julia> @objective(model, Max, 5x + 3 * y)
5 x + 3 y
```

Adding constraints is a lot like setting the objective. Here we create a
less-than-or-equal-to constraint using `<=`, but we can also create equality
constraints using `==` and greater-than-or-equal-to constraints with `>=`:
```jldoctest quickstart_example; filter=r"≤|<="
julia> @constraint(model, con, 1x + 5y <= 3)
con : x + 5 y <= 3.0
```
Note that in a similar manner to the `@variable` macro, we have named the
constraint `con`. This will bind the constraint to the Julia variable `con` for
later analysis.

Models are solved with the `JuMP.optimize!` function:
```jldoctest quickstart_example
julia> optimize!(model)
```

```@meta
DocTestSetup = quote
    # Now we load in the solution. Using a caching optimizer removes the need to
    # load a solver such as GLPK for building the documentation.
    mock = JuMP.backend(model).optimizer.model
    MOI.set(mock, MOI.TerminationStatus(), MOI.OPTIMAL)
    MOI.set(mock, MOI.PrimalStatus(), MOI.FEASIBLE_POINT)
    MOI.set(mock, MOI.DualStatus(), MOI.FEASIBLE_POINT)
    MOI.set(mock, MOI.ResultCount(), 1)
    MOI.set(mock, MOI.ObjectiveValue(), 10.6)
    MOI.set(mock, MOI.VariablePrimal(), JuMP.optimizer_index(x), 2.0)
    MOI.set(mock, MOI.VariablePrimal(), JuMP.optimizer_index(y), 0.2)
    MOI.set(mock, MOI.ConstraintDual(), JuMP.optimizer_index(con), -0.6)
    MOI.set(mock, MOI.ConstraintDual(), JuMP.optimizer_index(UpperBoundRef(x)), -4.4)
    MOI.set(mock, MOI.ConstraintDual(), JuMP.optimizer_index(LowerBoundRef(y)), 0.0)
end
```

After the call to `JuMP.optimize!` has finished, we need to query what happened.
The solve could terminate for a number of reasons. First, the solver might
have found the optimal solution or proved that the problem is infeasible.
However, it might also have run into numerical difficulties or terminated due
to a setting such as a time limit. We can ask the solver why it stopped using
the `JuMP.termination_status` function:
```jldoctest quickstart_example
julia> termination_status(model)
OPTIMAL::TerminationStatusCode = 1
```
In this case, `GLPK` returned `OPTIMAL`, this mean that it has found the optimal
solution.

```@meta
DocTestSetup = nothing
```

As the solver found an optimal solution, we expect the solution returned to be
a primal-dual pair of feasible solutions with zero duality gap.
We can verify the primal and dual status as follows to confirm this:
```jldoctest quickstart_example
julia> primal_status(model)
FEASIBLE_POINT::ResultStatusCode = 1

julia> dual_status(model)
FEASIBLE_POINT::ResultStatusCode = 1
```
Note that the primal and dual status only inform that the primal and dual
solutions are feasible and it is only because we verified that the termination
status is `OPTIMAL` that we can conclude that they form an optimal solution.

Finally, we can query the result of the optimization. First, we can query the
objective value:
```jldoctest quickstart_example
julia> objective_value(model)
10.6
```
We can also query the primal result values of the `x` and `y` variables:
```jldoctest quickstart_example
julia> value(x)
2.0

julia> value(y)
0.2
```

We can also query the value of the dual variable associated with the constraint
`con` (which we bound to a Julia variable when defining the constraint):
```jldoctest quickstart_example
julia> dual(con)
-0.6
```

!!! info
    See the [duality section](@ref constraint_duality) for a discussion
    of the convention that JuMP uses for signs of duals.

To query the dual variables associated with the variable bounds, things are a
little trickier as we first need to obtain a reference to the constraint:
```jldoctest quickstart_example; filter=r"≤|<="
julia> x_upper = UpperBoundRef(x)
x <= 2.0

julia> dual(x_upper)
-4.4
```
A similar process can be followed to obtain the dual of the lower bound
constraint on `y`:
```jldoctest quickstart_example; filter=r"≥|>="
julia> y_lower = LowerBoundRef(y)
y >= 0.0

julia> dual(y_lower)
0.0
```

```@meta
CurrentModule = JuMP
DocTestSetup = quote
    using JuMP
    const MOI = JuMP.MathOptInterface
end
```

# Querying Solutions

So far we have seen all the elements and constructs related to writing a JuMP
optimization model. In this section we reach the point of what to do with a
solved problem. Suppose your model is named `model`. Right after the call to
`JuMP.optimize!(model)` (which might take a while) it's possible to ask JuMP
questions about the finished optimization step. Typical questions include:
 - Why has the optimization process stopped?
 - Do I have a solution to my problem?
 - Is it optimal?
 - Do I have a dual solution?

JuMP follows closely the concepts defined in [MathOptInterface (MOI)](https://github.com/JuliaOpt/MathOptInterface.jl)
to answer user questions about a finished call to `JuMP.optimize!(model)`. There
are three main steps in querying a solution:

First, we can query the [`termination_status`](@ref) which will tell us why the
optimization stopped. This could be due to a number of reasons. For example, the
solver found an optimal solution, the problem was proven to be infeasible, or a
user-provided limit such as a time limit was encountered. For more information,
see the [Termination statuses](@ref) section below.

Second, we can query the [`primal_status`](@ref) and [`dual_status`](@ref),
which will tell us what kind of result do we have for our primal and dual
solution. This might be an optimal primal-dual pair, a primal solution without a
corresponding dual solution, or a certificate of primal or dual infeasibility.
For more information, see the [Solution statuses](@ref) section below.

Third, we can query [`JuMP.value`](@ref) and [`JuMP.dual`](@ref) to obtain the
primal and dual values of the optimization variables and constraints (if there
are values to be queried).

## Termination statuses

The reason why the optimization of `model` was finished is given by
```julia
JuMP.termination_status(model)
```

This function will return a `MOI.TerminationStatusCode` `enum`.

```@docs
MOI.TerminationStatusCode
```

## Solution statuses

These statuses indicate what kind of result is available to be queried
with [`JuMP.value`](@ref) and [`JuMP.dual`](@ref). Its possible that no result
is available to be queried.

We can obtain these statuses by calling [`JuMP.primal_status`](@ref) for the
primal status, and [`JuMP.dual_status`](@ref) for the dual status. Both will
return a `MOI.ResultStatusCode` `enum`.

```@docs
MOI.ResultStatusCode
```

Common status situations are described in the [`MathOptInterface` docs](http://www.juliaopt.org/MathOptInterface.jl/v0.8/apimanual/#Common-status-situations-1).

## Obtaining solutions

Provided the primal status is not (`MOI.NO_SOLUTION`), the primal solution can
be obtained by calling [`JuMP.value`](@ref). For the dual solution, the function
is [`JuMP.dual`](@ref). One fast way to check if the status is not
`MOI.NO_SOLUTION` is via [`JuMP.has_values`](@ref) for the primal status and
[`JuMP.has_duals`](@ref) for the dual solution.

It is important to note that if `has_values` returns false, calls to
[`JuMP.value`](@ref) and [`JuMP.dual`](@ref) might throw an error or return
arbitrary values.

The container type (e.g., scalar, vector, or matrix) of the returned solution
(primal or dual) depends on the type of the variable or constraint. See
[`AbstractShape`](@ref) and [`dual_shape`](@ref) for details.

To call [`JuMP.value`](@ref) or [`JuMP.dual`](@ref) on container of
[`VariableRef`](@ref) or [`ConstraintRef`](@ref), use the broadcast syntax,
e.g., `JuMP.value.(x)`.

The objective value of a solved problem can be obtained via
[`JuMP.objective_value`](@ref). The best known bound on the optimal objective
value can be obtained via [`JuMP.objective_bound`](@ref).

A recommended workflow for solving a model and querying the solution is the
following:
```julia
using JuMP, MathOptInterface
const MOI = MathOptInterface
model = Model()
@variable(model, x[1:10] >= 0)
# ... other constraints ...
JuMP.optimize!(model)

if JuMP.termination_status(model) == MOI.OPTIMAL
    optimal_solution = JuMP.value.(x)
    optimal_objective = JuMP.objective_value(model)
elseif JuMP.termination_status(model) == MOI.TIME_LIMIT && JuMP.has_values(model)
    suboptimal_solution = JuMP.value.(x)
    suboptimal_objective = JuMP.objective_value(model)
else
    error("The model was not solved correctly.")
end
```

## Reference

```@docs
JuMP.termination_status
JuMP.primal_status
JuMP.has_values
JuMP.value
JuMP.dual_status
JuMP.has_duals
JuMP.dual
MOI.optimize!
```

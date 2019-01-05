```@meta
CurrentModule = JuMP
DocTestSetup = quote
    using JuMP
end
```

Querying Solutions
==================

So far we have seen all the elements and constructs related to writing a JuMP
optimization model. In this section we reach the point of what to do
with a solved problem. Suppose your model is named `model`, right after the
call to `JuMP.optimize!(model)` (which might take a while) it's possible to
ask JuMP questions about the finished optimization step. Typical questions
include: Why has the optimization process stopped? Do I have a solution to my
problem? Is it optimal? Do I have a dual solution?

JuMP follows closely the concepts defined in [MathOptInterface (MOI)](https://github.com/JuliaOpt/MathOptInterface.jl)
to answer user questions about a finished call to `JuMP.optimize!(model)`.
There are 4 main steps in querying a solution:

* Firstly, we obtain the [`termination_status`](@ref) which will tell the user the reason why
optimization stopped. An optimal solution could have been found, or the problem
was proven to be infeasible, or the solver could not converge or JuMP got some
solver specific error etc.
* Secondly, we should figure out the [`primal_status`](@ref) and [`dual_status`](@ref), which will tell us what
kind of result do we have for our primal and dual solution. Maybe we have an
optimal primal-dual pair, or maybe just primal variable, or nothing at all,
or some infesability cetificate or even more.
* The third step is more straight forward, we should figure out if there really
is something available to be queryed in the primal solution with [`has_values`](@ref)
and in the dual soltion with [`has_duals`](@ref).
* Finally, we can [`JuMP.value`](@ref) and [`JuMP.dual`](@ref) to obtain primal and dual values of
the optimization variables and constraints.

## Termination Statuses

The possible reason why the optimization of `model` was finished is given by the call
```julia
JuMP.termination_status(model)
```

This function will return an instance of the `MOI.TerminationStatus` `enum`. All the possible outcome can be seen by calling help:
```julia
?MOI.TerminationStatus
```

```@docs
JuMP.termination_status
MOI.TerminationStatus
```

## Primal and Dual Statuses

We can obtain these statuses by callin the following functions on our model `model`:
```julia
JuMP.primal_status(model)
JuMP.dual_status(model)
```

Both will return an instance of the `MOI.ResultStatusCode` `enum`. All the possible outcome can be seen by calling help:
```julia
?MOI.ResultStatusCode
```

```@docs
JuMP.primal_status
JuMP.dual_status
MOI.ResultStatusCode
```

Common status situations are described in [`MathOptInterface`'s docs](http://www.juliaopt.org/MathOptInterface.jl/stable/apimanual/#Common-status-situations-1).

## Obtaining solutions

Solution are typically queried with the functions [`JuMP.value`](@ref) and [`JuMP.dual`](@ref).
These functions should be applied to references: [`VariableRef`](@ref) or a [`ConstraintRef`](@ref).
Depending on the type of variable and constraint associated with the reference the
return of the function is typically a scalar or an array but it can be an arbitrary object (see [`AbstractShape`](@ref) and [`dual_shape`](@ref)).

!!! warn
    If we desire to obtain the [`JuMP.value`](@ref)/[`JuMP.dual`](@ref) from a container of [`VariableRef`](@ref)/[`ConstraintRef`](@ref)  we should use
    the broadcast syntax.

```@docs
JuMP.has_values
JuMP.value
JuMP.has_duals
JuMP.dual
JuMP.shadow_price
```
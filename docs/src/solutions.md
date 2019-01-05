```@meta
CurrentModule = JuMP
DocTestSetup = quote
    using JuMP
    const MOI = JuMP.MathOptInterface
end
```

Querying Solutions
==================

So far we have seen all the elements and constructs related to writing a JuMP
optimization model. In this section we reach the point of what to do
with a solved problem. Suppose your model is named `model`. Right after the
call to `JuMP.optimize!(model)` (which might take a while) it's possible to
ask JuMP questions about the finished optimization step. Typical questions
include: Why has the optimization process stopped? Do I have a solution to my
problem? Is it optimal? Do I have a dual solution?

JuMP follows closely the concepts defined in [MathOptInterface (MOI)](https://github.com/JuliaOpt/MathOptInterface.jl)
to answer user questions about a finished call to `JuMP.optimize!(model)`.
There are 3 main steps in querying a solution:

* Firstly, we obtain the [`termination_status`](@ref) which will tell the user the reason why
optimization stopped. An optimal solution could have been found, or the problem
was proven to be infeasible, or the solver could not converge or JuMP got some
solver specific error etc.
* Secondly, we should figure out the [`primal_status`](@ref) and [`dual_status`](@ref), which will tell us what
kind of result do we have for our primal and dual solution. Maybe we have an
optimal primal-dual pair, or maybe just primal variable, or nothing at all,
or some infesability cetificate or even more.
* Finally, we can [`JuMP.value`](@ref) and [`JuMP.dual`](@ref) to obtain primal and dual values of
the optimization variables and constraints (if there are values to be queryed).

## Termination Statuses

The possible reason why the optimization of `model` was finished is given by the call
```julia
JuMP.termination_status(model)
```

This function will return an instance of the `MOI.TerminationStatus` `enum`. All the possible outcomes can be seen by calling help:
```julia
?MOI.TerminationStatus
```

```@docs
JuMP.termination_status
MOI.TerminationStatus
```

TODO test these strings

## Primal and Dual Statuses

These statuses indicate what kind of result is available to be queried
with [`JuMP.value`](@ref) and [`JuMP.dual`](@ref). Its possible that no result
is available to be queried.
We can obtain these statuses by calling the following functions on our model `model`:
```julia
JuMP.primal_status(model)
JuMP.dual_status(model)
```

Both will return an instance of the `MOI.ResultStatusCode` `enum`. All the possible outcomes can be seen by calling help:
```julia
?MOI.ResultStatusCode
```

```@docs
JuMP.primal_status
JuMP.dual_status
MOI.ResultStatusCode
```

TODO test these strings

Common status situations are described in [`MathOptInterface` docs](http://www.juliaopt.org/MathOptInterface.jl/v0.8/apimanual/#Common-status-situations-1).

## Obtaining solutions

Solution are queried with the functions [`JuMP.value`](@ref) and [`JuMP.dual`](@ref).
It is important to note that there might not be a solution available, consequently
calls to [`JuMP.value`](@ref) and [`JuMP.dual`](@ref) might return errors or not meaningful
results. One fast way to check availability of results is to use the functions:
[`JuMP.has_values`](@ref) and [`JuMP.has_duals`](@ref).
[`JuMP.value`](@ref) and [`JuMP.dual`](@ref) should be applied to references: [`VariableRef`](@ref)
or a [`ConstraintRef`](@ref).
Depending on the type of variable and constraint associated with the reference the
return of the function is typically a scalar or an array but it can be an arbitrary
object (see [`AbstractShape`](@ref) and [`dual_shape`](@ref)).

!!! info
    To obtain the [`JuMP.value`](@ref)/[`JuMP.dual`](@ref) from a container of [`VariableRef`](@ref)/[`ConstraintRef`](@ref)  we should use
    the broadcast syntax.

TODO: add example of broadcasting.

```@docs
JuMP.has_values
JuMP.value
JuMP.has_duals
JuMP.dual
JuMP.shadow_price
```
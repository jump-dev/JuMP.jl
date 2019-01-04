Querying Solutions
==================

So far we have seen all the elements and constructs related to writing a JuMP
optimization model, in this section we reach the point of what to do
with a solved problem. Suppose your model is named `model`, right after the
call to `JuMP.optimize!(model)` (which might take a while) its possible to
ask JuMP questions about the finished optimization step. Typical questions
include: Why has the optimization process stopped? Do I have a solution to my
problem? Is it optimal? Do I have a dual solution?

JuMP follows closely the concepts defined in [MathOptInterface (MOI)](https://github.com/JuliaOpt/MathOptInterface.jl)
to answer user questions about a finished call to `JuMP.optimize!(model)`.
There are 4 main steps in querying a solution:

* Firstly, we obtain the `termination_status` which will tell the user the reason why
optimization stopped. An optimal solution could have been found, or the problem
was proven to be infeasible, or the solver could not converge or JuMP got some
solver specific error etc.
* Secondly, we should figure out the `primal_status` and `dual_status`, which will tell us what
kind of result do we have for our primal and dual solution. Maybe we have an
optimal primal-dual pair, or maybe just primal variable, or nothing at all,
or some infesability cetificate or even more.
* The third step is more straight forward, we should figure out if there really
is something available to be queryed in the primal solution with `has_values`
and in the dual soltion with `has_duals`.
* Finally, we can `JuMP.value` and `JuMP.dual` to obtain primal and duas values of
the optimization variables and constraints.

## Termination Statuses

The possible reason why the optimization of `model` was finished is is given by the call

```julia
    JuMP.termination_status(model)
```

This function will return an instance of the `MOI.TerminationStatus` `enum` with one of
the following values:

TODO
Do we add an exautive table here?
We can check the validity of the statuses with doctests
but how do we make sure the list is complete?

Or we put docstrings here?

## Primal and Dual Statuses

We can obtain these statuses by callin the following functions on our model `model`:

```julia
    JuMP.primal_status(model)
    JuMP.dual_status(model)
```

Both will return an instance of the `MOI.ResultStatusCode` `enum` with one of
the following values:

TODO
same as above

TODO
do we add the tables from MOI? (they are very useful)

## Obtaining solutions

Solution are typically queried with the functions `JuMP.value` and `JuMP.dual`.
These functions should be applied to references: `VariableRef` or a `ConstraintRef`.
Depending on the type of variable and constraint associated with the reference the
return of the function, it can be a scalar or an array.

If we desire to obtain the `JuMP.value` of a vector of `VariableRef` we should use
the broadcast syntax.
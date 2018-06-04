How do I ...? (FAQ)
===================

**Q: I'm using a solver that supports warm-starts for integer problems. How do I communicate the initial solution to the solver?**

TODO: use `start=` in `@variable` or the `VariablePrimalStart()` attribute.

**Q: How can I suppress output in a solver-independent way?**

TODO: Update answer for JuMP 0.19.

A: JuMP does not currently support generic parameters independent of the chosen
solver object. To suppress output with Gurobi, for example, one would say

```julia
m = Model(solver=GurobiSolver(OutputFlag=0))
```

When a solver is not specified, i.e., the model is created with ``m = Model()``,
 there's no option to suppress output. A workaround is to redirect STDOUT before
 and after the call to ``solve(m)``:

```julia
TT = STDOUT # save original STDOUT stream
redirect_stdout()
solve(m)
redirect_stdout(TT) # restore STDOUT
```

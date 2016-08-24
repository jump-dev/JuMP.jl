JuMP
====
#### Frequently Asked Questions

**Q: I'm using a solver that supports warm-starts for integer problems. How do I communicate the initial solution to the solver?**

A: You can use the ``setvalue`` function on a variable to set its initial value, e.g.

```julia
@variable(m, x[1:3], Int)
setvalue(x[1], 10.0)
setvalue(x[3], 2.0)
```

If you don't set the value of a variable, the behavior is solver-dependent. Some solvers attempt to "complete" solutions that aren't completely specified.

**Q: How can I suppress output in a solver-independent way?**

A: JuMP does not currently support generic parameters independent of the chosen solver object. To suppress output with Gurobi, for example, one would say

```julia
m = Model(solver=GurobiSolver(OutputFlag=0))
```

When a solver is not specified, i.e., the model is created with ``m = Model()``, there's no option to suppress output. A workaround is to redirect STDOUT before and after the call to ``solve(m)``:

```julia
TT = STDOUT # save original STDOUT stream
redirect_stdout()
solve(m)
redirect_stdout(TT) # restore STDOUT
```

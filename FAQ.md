JuMP
====
#### Frequently Asked Questions

**Q: I'm using a solver that supports warm-starts for integer problems. How do I communicate the initial solution to the solver?**

A: You can use the ``setValue!`` function on a variable to set its initial value, e.g.

```julia
@defVar(m, x[1:3], Int)
setValue!(x[1], 10.0)
setValue!(x[3], 2.0)
```

If you don't set the value of a variable, the behaviour is solver-dependent. Some solvers attempt to "complete" solutions that aren't completely specified.



# Development roadmap

This page is not JuMP documentation *per se* but are notes for the JuMP
community. The JuMP developers have compiled this roadmap document to
share their plans and goals. Contributions to roadmap issues are especially
invited.

## JuMP 1.0

JuMP 1.0 will be ready to release roughly when all of these tasks are completed.
Some but not all of these tasks are summarized in the
[JuMP 1.0 milestone](https://github.com/jump-dev/JuMP.jl/milestone/12).

- Create a website for JuMP (**Done**: [jump.dev](https://jump.dev))
- Deprecate the JuliaOpt organization and move repositories to the
  [JuMP-dev](https://github.com/JuMP-dev) organization (**Done**)
- Address major regressions from JuMP 0.18
  - Performance ([#1403](https://github.com/jump-dev/JuMP.jl/issues/1403),
                 [#1654](https://github.com/jump-dev/JuMP.jl/issues/1654),
                 [#1607](https://github.com/jump-dev/JuMP.jl/issues/1607))
  - Callbacks (**Done**: see `examples/callbacks.jl`)
  - Column generation syntax (**Done**: see `examples/cutting_stock_column_generation.jl`)
  - Support for second-order cones in Gurobi, CPLEX, and Xpress (**Done**)
- Fix issues that we promised MOI would fix
  - Checking feasibility of solutions (**Done**: [#2466](https://github.com/jump-dev/JuMP.jl/pull/2466))
  - Accessing IIS (**Done**: see [Conflicts](@ref))
  - Accessing multiple results from solvers (**Done**: [Gurobi#392](https://github.com/jump-dev/Gurobi.jl/pull/392))
  - Dual warm-starts (**Done**: [#2214](https://github.com/jump-dev/JuMP.jl/pull/2214))
- Address "easy" usability issues
  - Line numbers in error messages (**Done**: [#2276](https://github.com/jump-dev/JuMP.jl/pull/2276))
  - LP sensitivity summary (**Done**: see [Sensitivity analysis for LP](@ref))
  - Inferred element types for collections in macros (**Done**: [#2070](https://github.com/jump-dev/JuMP.jl/pull/2070))
  - Expose solver-independent options from JuMP (**Done**: see [`set_silent`](@ref) etc.)
- Improve the documentation ([#1062](https://github.com/jump-dev/JuMP.jl/issues/1062))
  - Separate how-to, concept explanation, and technical reference following the
    [Divio recommendations](https://www.divio.com/blog/documentation/) (**Done**)
  - Fully integrate [JuMPTutorials](https://github.com/jump-dev/JuMPTutorials.jl)
    with JuMP's documentation (**Done**)
- Developer experience
  - Get JuMP's unit tests running faster. See [#1745](https://github.com/jump-dev/JuMP.jl/pull/1745). (**Done**)
- All solvers should complete the transition to MOI (**Done**)
- Provide packages for installing Bonmin and Couenne (**Done**)
- [MathOptFormat](https://github.com/odow/MathOptFormat.jl) 1.0 (**Done**)

## MOI 1.0

```@meta
# TODO: List MOI 1.0 items here.
```

## Beyond JuMP 1.0

```@meta
# TODO: Copy over list of items not tied to JuMP 1.0. These should have more
# elaborate explanations so that potential contributors know what we mean,
# i.e., a few sentences each or a link to a document/issue.
```

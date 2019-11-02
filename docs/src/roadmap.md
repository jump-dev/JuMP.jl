# Development roadmap

This page is not JuMP documentation *per se* but are notes for the JuMP
community. The JuMP developers have compiled this roadmap document to
share their plans and goals. Contributions to roadmap issues are especially
invited.

## JuMP 1.0

JuMP 1.0 will be ready to release roughly when all of these tasks are completed.
Some but not all of these tasks are summarized in the
[JuMP 1.0 milestone](https://github.com/JuliaOpt/JuMP.jl/milestone/12).

- Create a website for JuMP.
- Deprecate the JuliaOpt organization and move repositories to the
  [JuMP-dev](https://github.com/JuMP-dev) organization.
- Address major regressions from JuMP 0.18.
  - Performance ([#1403](https://github.com/JuliaOpt/JuMP.jl/issues/1403),
                 [#1654](https://github.com/JuliaOpt/JuMP.jl/issues/1654),
                 [#1607](https://github.com/JuliaOpt/JuMP.jl/issues/1607))
  - Callbacks
  - Column generation syntax (**Done**: see `examples/cutting_stock_column_generation.jl`)
  - Support for second-order cones in Gurobi, CPLEX, and Xpress.
- Fix issues that we promised MOI would fix.
  - Checking feasibility of solutions ([#693](https://github.com/JuliaOpt/JuMP.jl/issues/693))
  - Accessing IIS ([#1053](https://github.com/JuliaOpt/JuMP.jl/issues/1035))
  - Accessing multiple results from solvers
  - Dual warm-starts ([#2094](https://github.com/JuliaOpt/JuMP.jl/issues/2094))
- Address “easy” usability issues
  - Line numbers in error messages ([#1174](https://github.com/JuliaOpt/JuMP.jl/issues/1174))
  - LP sensitivity summary (**Done**: see [Sensitivity analysis for LP](@ref))
  - Inferred element types for collections in macros (**Done**: [#2070](https://github.com/JuliaOpt/JuMP.jl/pull/2070))
  - Expose solver-independent options from JuMP (**Done**: see [`set_silent`](@ref) etc.)
- Improve the documentation ([#1062](https://github.com/JuliaOpt/JuMP.jl/issues/1062))
  - Separate how-to, concept explanation, and technical reference following the
    [Divio recommendations](https://www.divio.com/blog/documentation/)
  - Fully integrate [JuMPTutorials](https://github.com/JuliaOpt/JuMPTutorials.jl)
    with JuMP's documentation.
- Developer experience
  - Get JuMP’s unit tests running in less than two minutes. See [#1745](https://github.com/JuliaOpt/JuMP.jl/pull/1745).
- All solvers should complete the transition to MOI.
- Provide packages for installing Bonmin and Couenne.
- [MathOptFormat](https://github.com/odow/MathOptFormat.jl) 1.0

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

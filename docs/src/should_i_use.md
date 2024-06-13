# Should you use JuMP?

JuMP is an [algebraic modeling language](@ref algebraic-modeling-language) for
mathematical optimization written in the [Julia language](https://julialang.org).

This page explains when you should consider using JuMP, and importantly, when
you should _not_ use JuMP.

## When should you use JuMP?

You should use JuMP if you have a constrained optimization problem that is
formulated using the language of mathematical programming, that is, the problem
has:

 * a set of real- or complex-valued decision variables
 * a scalar- or vector-valued real objective function
 * a set of constraints.

Key reasons to use JuMP include:

 - User friendliness
   - JuMP has syntax that mimics natural mathematical expressions. (See the
     section on [algebraic modeling languages](@ref algebraic-modeling-language).)
 - Solver independence
   - JuMP uses a generic solver-independent interface provided by the
     [MathOptInterface](https://github.com/jump-dev/MathOptInterface.jl)
     package, making it easy to change between a number of open-source and
     commercial optimization software packages ("solvers"). The
     [Supported solvers](@ref) section contains a table of the currently
     supported solvers.
 - Ease of embedding
   - JuMP itself is written purely in Julia. Solvers are the only binary
     dependencies.
   - JuMP provides automatic installation of most solvers.
   - Because it is embedded in a general-purpose programming language, JuMP
     makes it easy to solve optimization problems as part of a larger workflow,
     for example, inside a simulation, behind a web server, or as a subproblem
     in a decomposition algorithm. As a trade-off, JuMP's syntax is constrained
     by the syntax and functionality available in Julia.
   - JuMP is [MPL](https://www.mozilla.org/MPL/2.0/) licensed, meaning that it
     can be embedded in commercial software that complies with the terms of the
     license.
 - Speed
   - Benchmarking has shown that JuMP can create problems at similar speeds to
     special-purpose modeling languages such as [AMPL](https://ampl.com/).
   - JuMP communicates with most solvers in memory, avoiding the need to write
     intermediary files.
 - Access to advanced algorithmic techniques
   - JuMP supports efficient _in-memory_ re-solves of models.
   - JuMP provides access to solver-independent and solver-dependent
     [Callbacks](@ref callbacks_manual).

## When should you not use JuMP?

JuMP supports a broad range of optimization classes. However, there are still
some that it doesn't support, or that are better supported by other software
packages.

### You want to optimize a complicated Julia function

Packages in Julia compose well. It's common for people to pick two unrelated
packages and use them in conjunction to create novel behavior. JuMP isn't one of
those packages.

If you want to optimize an ordinary differential equation from
[DifferentialEquations.jl](https://github.com/SciML/DifferentialEquations.jl)
or tune a neural network from [Flux.jl](https://github.com/FluxML/Flux.jl),
consider using other packages such as:

 * [Optim.jl](https://github.com/JuliaNLSolvers/Optim.jl)
 * [Optimization.jl](https://github.com/SciML/Optimization.jl)
 * [NLPModels.jl](https://github.com/JuliaSmoothOptimizers/NLPModels.jl)
 * [Nonconvex.jl](https://github.com/JuliaNonconvex/Nonconvex.jl)

### Black-box, derivative free, or unconstrained optimization

JuMP does support nonlinear programs with constraints and objectives containing
user-defined operators. However, the functions must be automatically
differentiable, or need to provide explicit derivatives. (See
[User-defined operators](@ref jump_user_defined_operators) for more information.)

If your function is a black-box that is non-differentiable (for example, it is
the output of a simulation written in C++), JuMP is not the right tool for the
job. This also applies if you want to use a derivative free method.

Even if your problem is differentiable, if it is unconstrained there is limited
benefit (and downsides in the form of more overhead) to using JuMP over tools
which are only concerned with function minimization.

Alternatives to consider are:

 * [Optim.jl](https://github.com/JuliaNLSolvers/Optim.jl)
 * [Optimization.jl](https://github.com/SciML/Optimization.jl)
 * [NLopt.jl](https://github.com/JuliaOpt/NLopt.jl)

### Disciplined convex programming

JuMP does not support [disciplined convex programming (DCP)](https://dcp.stanford.edu).

Alternatives to consider are:

 * [Convex.jl](https://github.com/jump-dev/Convex.jl)
 * [CVXPY [Python]](https://github.com/cvxpy/cvxpy)
 * [YALMIP [MATLAB]](https://yalmip.github.io)

!!! note
    `Convex.jl` is also built on MathOptInterface, and shares the same set of
    underlying solvers. However, you input problems differently, and Convex.jl
    checks that the problem is DCP.

### Stochastic programming

JuMP requires deterministic input data.

If you have stochastic input data, consider using a JuMP extension such as:

 * [InfiniteOpt.jl](https://github.com/infiniteopt/InfiniteOpt.jl)
 * [StochasticPrograms.jl](https://github.com/martinbiel/StochasticPrograms.jl)
 * [SDDP.jl](https://github.com/odow/SDDP.jl)

### Polyhedral computations

JuMP does not provide tools for working with the polyhedron formed by the set
of linear constraints.

Alternatives to consider are:

 * [Polyhedra.jl](https://github.com/JuliaPolyhedra/Polyhedra.jl) (See the [documentation](https://juliapolyhedra.github.io/Polyhedra.jl/v0.7.6/optimization/#Creating-a-polyhedron-from-the-feasible-set-of-a-JuMP-model)
   to create a polyhedron from a JuMP model.)

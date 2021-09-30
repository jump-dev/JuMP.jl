# Should I use JuMP?

JuMP is an [algebraic modeling language](@ref algebraic-modeling-language) for
mathematical optimization written in the [Julia language](https://julialang.org).

## When should I use JuMP?

You should use JuMP if you have a constrained optimization problem for which you
can formulate a set of decision variables, a scalar objective function, and a
set of constraints.

Key reasons to use JuMP include:

 - User friendliness
   - Syntax that mimics natural mathematical expressions. (See the section on
     [algebraic modeling languages](@ref algebraic-modeling-language).)
 - Speed
   - Benchmarking has shown that JuMP can create problems at similar speeds to
     special-purpose modeling languages such as [AMPL](https://ampl.com/).
   - JuMP communicates with most solvers in memory, avoiding the need to write
     intermediary files.
 - Solver independence
   - JuMP uses a generic solver-independent interface provided by the
     [MathOptInterface](https://github.com/jump-dev/MathOptInterface.jl)
     package, making it easy to change between a number of open-source and
     commercial optimization software packages ("solvers"). The
     [Supported solvers](@ref) section contains a table of the currently
     supported solvers.
 - Access to advanced algorithmic techniques
   - Efficient _in-memory_ LP re-solves which previously required using
     solver-specific and/or low-level C++ libraries.
   - Access to solver-independent and solver-dependent [Callbacks](@ref callbacks_manual).
 - Ease of embedding
   - JuMP itself is written purely in Julia. Solvers are the only binary
     dependencies.
   - Automated install of many solver dependencies.
     - JuMP provides automatic installation of many open-source solvers. This is
       different to modeling languages in Python which require you to download
       and install a solver yourself.
   - Being embedded in a general-purpose programming language makes it easy to
     solve optimization problems as part of a larger workflow (e.g., inside a
     simulation, behind a web server, or as a subproblem in a decomposition
     algorithm).
     - As a trade-off, JuMP's syntax is constrained by the syntax available in
       Julia.
   - JuMP is [MPL](https://www.mozilla.org/MPL/2.0/) licensed, meaning that it
     can be embedded in commercial software that complies with the terms of the
     license.

## When should I not use JuMP?

JuMP supports a broad range of optimization classes. However, there are still
some that it doesn't support, or that are better supported by other software
packages.

### Black-box, derivative free, or unconstrained optimization

JuMP does support nonlinear programs with constraints and objectives containing
user-defined functions. However, the functions must be automatically
differentiable, or need to provide explicit derivatives. (See
[User-defined Functions](@ref) for more information.)

If your function is a black-box that is non-differentiable (e.g., the output of
a simulation written in C++), JuMP is not the right tool for the job. This also
applies if you want to use a derivative free method.

Even if your problem is differentiable, if it is unconstrained there is limited
benefit (and downsides in the form of more overhead) to using JuMP over tools
which are only concerned with function minimization.

Alternatives to consider are:
 * [Optim.jl](https://github.com/JuliaNLSolvers/Optim.jl)
 * [NLopt.jl](https://github.com/JuliaOpt/NLopt.jl)

### Multiobjective programs

If your problem has more than one objective, JuMP is not the right tool for the
job. However, [we're working on fixing this!](https://github.com/jump-dev/JuMP.jl/issues/2099).

Alternatives to consider are:
 * [vOptGeneric.jl](https://github.com/vOptSolver/vOptGeneric.jl)

### Disciplined convex programming

JuMP does not support [disciplined convex programming (DCP)](https://dcp.stanford.edu).

Alternatives to consider are:
 * [Convex.jl](https://github.com/jump-dev/Convex.jl)

!!! note
    `Convex.jl` is also built on MathOptInterface, and shares the same set of
    underlying solvers. However, you input problems differently, and Convex.jl
    checks that the problem is DCP.

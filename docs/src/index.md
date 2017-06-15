JuMP --- Julia for Mathematical Optimization
============================================

[JuMP](https://github.com/JuliaOpt/JuMP.jl) is a domain-specific modeling language for [mathematical optimization](http://en.wikipedia.org/wiki/Mathematical_optimization) embedded in [Julia](http://julialang.org/). It currently supports a number of open-source and commercial solvers (see below) for a variety of problem classes, including **linear programming**, **mixed-integer programming**, **second-order conic programming**, **semidefinite programming**, and **nonlinear programming**. JuMP's features include:

-   User friendliness
    -   Syntax that mimics natural mathematical expressions.
    -   Complete documentation.
-   Speed
    -   Benchmarking has shown that JuMP can create problems at similar speeds to special-purpose modeling languages such as [AMPL](http://www.ampl.com/).
    -   JuMP communicates with solvers in memory, avoiding the need to write intermediary files.
-   Solver independence
    -   JuMP uses a generic solver-independent interface provided by the [MathProgBase](https://github.com/mlubin/MathProgBase.jl) package, making it easy to change between a number of open-source and commercial optimization software packages ("solvers").
    -   Currently supported solvers include [Artelys Knitro](http://artelys.com/en/optimization-tools/knitro), [Bonmin](https://projects.coin-or.org/Bonmin), [Cbc](https://projects.coin-or.org/Cbc), [Clp](https://projects.coin-or.org/Clp), [Couenne](https://projects.coin-or.org/Couenne), [CPLEX](http://www-01.ibm.com/software/commerce/optimization/cplex-optimizer/), [ECOS](https://github.com/ifa-ethz/ecos), [FICO Xpress](http://www.fico.com/en/products/fico-xpress-optimization-suite), [GLPK](http://www.gnu.org/software/glpk/), [Gurobi](http://www.gurobi.com), [Ipopt](https://projects.coin-or.org/Ipopt), [MOSEK](http://www.mosek.com/), [NLopt](http://ab-initio.mit.edu/wiki/index.php/NLopt), and [SCS](https://github.com/cvxgrp/scs).
-   Access to advanced algorithmic techniques
    -   Including efficient LP re-solves &lt;probmod&gt; and callbacks for mixed-integer programming &lt;callbacks&gt; which previously required using solver-specific and/or low-level C++ libraries.
-   Ease of embedding
    -   JuMP itself is written purely in Julia. Solvers are the only binary dependencies.
    -   Being embedded in a general-purpose programming language makes it easy to solve optimization problems as part of a larger workflow (e.g., inside a simulation, behind a web server, or as a subproblem in a decomposition algorithm).
        -   As a trade-off, JuMP's syntax is constrained by the syntax available in Julia.
    -   JuMP is [MPL](https://www.mozilla.org/MPL/2.0/) licensed, meaning that it can be embedded in commercial software that complies with the terms of the license.

While neither Julia nor JuMP have reached version 1.0 yet, the releases are stable enough for everyday use and are being used in a number of research projects and neat applications by a growing community of users who are early adopters. JuMP remains under active development, and we welcome your feedback, suggestions, and bug reports.

Installing JuMP
---------------

If you are familiar with Julia you can get started quickly by using the package manager to install JuMP:

    julia> Pkg.add("JuMP")

And a solver, e.g.:

    julia> Pkg.add("Clp")  # Will install Cbc as well

Then read the quick-start and/or see a simple-example. The subsequent sections detail the complete functionality of JuMP.

Contents
--------

```@contents
Pages = ["installation.md",
    "quickstart.md",
    "refmodel.md",
    "refvariable.md",
    "refexpr.md",
    "probmod.md",
    "callbacks.md",
    "nlp.md"]
Depth = 2
```

### Citing JuMP

If you find JuMP useful in your work, we kindly request that you cite the following [paper](http://dx.doi.org/10.1137/15M1020575):

``` sourceCode
@article{DunningHuchetteLubin2017,
author = {Iain Dunning and Joey Huchette and Miles Lubin},
title = {JuMP: A Modeling Language for Mathematical Optimization},
journal = {SIAM Review},
volume = {59},
number = {2},
pages = {295-320},
year = {2017},
doi = {10.1137/15M1020575},
}
```

A preprint of this paper is freely available on [arXiv](http://arxiv.org/abs/1508.01982).

For an earlier work where we presented a prototype implementation of JuMP, see [here](http://dx.doi.org/10.1287/ijoc.2014.0623):

``` sourceCode
@article{LubinDunningIJOC,
author = {Miles Lubin and Iain Dunning},
title = {Computing in Operations Research Using Julia},
journal = {INFORMS Journal on Computing},
volume = {27},
number = {2},
pages = {238-248},
year = {2015},
doi = {10.1287/ijoc.2014.0623},
}
```

A preprint of this paper is also [freely available](http://arxiv.org/abs/1312.1431).

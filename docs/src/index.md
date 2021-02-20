![JuMP logo](assets/jump-logo-with-text.svg)
===

[![Powered by NumFOCUS](https://img.shields.io/badge/powered%20by-NumFOCUS-orange.svg?style=flat&colorA=E1523D&colorB=007D8A)](https://numfocus.org)

```@meta
# These comments do not display in the HTML output.
# See https://github.com/JuliaDocs/Documenter.jl/issues/674.
```

!!! warning
    Between versions 0.18 and 0.19, JuMP underwent a major transition in its
    underlying solver abstraction API, from
    [MathProgBase](https://github.com/JuliaOpt/MathProgBase.jl) to
    [MathOptInterface](https://github.com/jump-dev/MathOptInterface.jl).
    See [NEWS.md](https://github.com/jump-dev/JuMP.jl/blob/master/NEWS.md) for
    a comprehensive list of changes between the two versions, many of which
    are breaking. This documentation is for JuMP/MathOptInterface.
    For the documentation of JuMP 0.18, see
    [here](https://jump.dev/JuMP.jl/0.18/).

## What is JuMP?

[JuMP](https://github.com/jump-dev/JuMP.jl) is a domain-specific modeling
language for [mathematical optimization](https://en.wikipedia.org/wiki/Mathematical_optimization)
embedded in [Julia](https://julialang.org/). It currently supports a number of
open-source and commercial solvers for a variety of problem classes, including
linear, mixed-integer, second-order conic, semidefinite, and nonlinear
programming.

## Resources for getting started

* Checkout the [Installation Guide](@ref).
* Read the introductory tutorials [Getting started with Julia](@ref) and
  [Getting started with JuMP](@ref).
* Browse some of our examples, including classics such as [The diet problem](@ref),
  or the [Maximum likelihood estimation](@ref) problem using nonlinear
  programming.
* Work through more in-depth tutorials at [JuMPTutorials.jl](https://github.com/jump-dev/JuMPTutorials.jl)

!!! tip
    Need help? Join the [community forum](https://discourse.julialang.org/c/domain/opt/13)
    to search for questions to commonly asked questions.

    Before asking a question, make sure to read the post [make it easier to help you](https://discourse.julialang.org/t/psa-make-it-easier-to-help-you/14757),
    which contains a number of tips on how to ask a good question.

## Key features

JuMP's features include:

-   User friendliness
    -   Syntax that mimics natural mathematical expressions.
    -   Complete documentation (WIP!)
-   Speed
    -   Benchmarking has shown that JuMP can create problems at similar speeds
        to special-purpose modeling languages such as
        [AMPL](https://ampl.com/).
    -   JuMP communicates with most solvers in memory, avoiding the need to
        write intermediary files.
-   Solver independence
    -   JuMP uses a generic solver-independent interface provided by the
        [MathOptInterface](https://github.com/jump-dev/MathOptInterface.jl)
        package, making it easy to change between a number of open-source and
        commercial optimization software packages ("solvers"). The [Supported solvers](@ref)
        section contains a table of the currently supported solvers.
-   Access to advanced algorithmic techniques
    -   Including efficient LP re-solves which previously required using
        solver-specific and/or low-level C++ libraries.
-   Ease of embedding
    -   JuMP itself is written purely in Julia. Solvers are the only binary
        dependencies.
    -   Being embedded in a general-purpose programming language makes it easy
        to solve optimization problems as part of a larger workflow (e.g.,
        inside a simulation, behind a web server, or as a subproblem in a
        decomposition algorithm).
        -   As a trade-off, JuMP's syntax is constrained by the syntax available
            in Julia.
    -   JuMP is [MPL](https://www.mozilla.org/MPL/2.0/) licensed, meaning that
        it can be embedded in commercial software that complies with the terms
        of the license.

### Citing JuMP

If you find JuMP useful in your work, we kindly request that you cite the
following paper ([pdf](https://mlubin.github.io/pdf/jump-sirev.pdf)):

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

For an earlier work where we presented a prototype implementation of JuMP, see
[here](https://dx.doi.org/10.1287/ijoc.2014.0623):

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

A preprint of this paper is [freely available](https://arxiv.org/abs/1312.1431).

---

![NumFOCUS logo](assets/numfocus-logo.png)

JuMP is a Sponsored Project of NumFOCUS, a 501(c)(3) nonprofit charity in the
United States. NumFOCUS provides JuMP with fiscal, legal, and administrative
support to help ensure the health and sustainability of the project. Visit
[numfocus.org](https://numfocus.org) for more information.

You can support JuMP by [donating](https://numfocus.salsalabs.org/donate-to-jump/index.html).

Donations to JuMP are managed by NumFOCUS. For donors in the United States,
your gift is tax-deductible to the extent provided by law. As with any donation,
you should consult with your tax adviser about your particular tax situation.

JuMP's largest expense is the annual JuMP-dev workshop. Donations will help us
provide travel support for JuMP-dev attendees and take advantage of other
opportunities that arise to support JuMP development.

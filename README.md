![JuMP logo](https://www.juliaopt.org/images/jump-logo-with-text.svg "JuMP logo")
---

[![Powered by NumFOCUS](https://img.shields.io/badge/powered%20by-NumFOCUS-orange.svg?style=flat&colorA=E1523D&colorB=007D8A)](https://numfocus.org)

JuMP is a domain-specific modeling language for **[mathematical optimization]**
embedded in **[Julia]**. It currently supports a number of open-source and
commercial solvers ([Artelys Knitro], [BARON], [Bonmin], [Cbc], [CDCS], [CDD],
[Clp], [COSMO], [Couenne], [CPLEX], [CSDP], [ECOS], [FICO Xpress], [GLPK],
[Gurobi], [Ipopt], [Juniper], [MOSEK], [NLopt], [OSQP], [ProxSDP], [SCIP],
[SCS], [SDPA], [SDPT3], [SeDuMi], [Tulip]) for a variety of problem classes, including
**[linear programming]**, **[(mixed) integer programming]**,
**[second-order conic programming]**, **[semidefinite programming]**, and **[nonlinear programming]**.

[mathematical optimization]: https://en.wikipedia.org/wiki/Mathematical_optimization
[Julia]: https://julialang.org/
[Artelys Knitro]: https://artelys.com/en/optimization-tools/knitro
[BARON]: http://archimedes.cheme.cmu.edu/?q=baron
[Bonmin]: https://projects.coin-or.org/Bonmin
[Cbc]: https://github.com/coin-or/Cbc
[CDCS]: https://github.com/oxfordcontrol/CDCS
[CDD]: https://github.com/cddlib/cddlib
[Clp]: https://github.com/coin-or/Clp
[COSMO]: https://github.com/oxfordcontrol/COSMO.jl
[Couenne]: https://projects.coin-or.org/Couenne
[CPLEX]: https://www.ibm.com/analytics/cplex-optimizer
[CSDP]: https://projects.coin-or.org/Csdp/
[ECOS]: https://github.com/ifa-ethz/ecos
[FICO Xpress]: https://www.fico.com/en/products/fico-xpress-optimization-suite
[GLPK]: http://www.gnu.org/software/glpk/
[Gurobi]: https://www.gurobi.com/
[Ipopt]: https://github.com/coin-or/Ipopt
[Juniper]: https://github.com/lanl-ansi/Juniper.jl
[MOSEK]: https://mosek.com/
[NLopt]: https://nlopt.readthedocs.io/en/latest/
[OSQP]: https://osqp.org/
[ProxSDP]: https://github.com/mariohsouto/ProxSDP.jl
[SCIP]: https://scip.zib.de/
[SCS]: https://github.com/cvxgrp/scs
[SDPA]: http://sdpa.sourceforge.net/
[SDPT3]: https://blog.nus.edu.sg/mattohkc/softwares/sdpt3/
[SeDuMi]: http://sedumi.ie.lehigh.edu/
[Tulip]: https://github.com/ds4dm/Tulip.jl
[linear programming]: https://en.wikipedia.org/wiki/Linear_programming
[(mixed) integer programming]: https://en.wikipedia.org/wiki/Integer_programming
[second-order conic programming]: https://en.wikipedia.org/wiki/Second-order_cone_programming
[semidefinite programming]: https://en.wikipedia.org/wiki/Semidefinite_programming
[nonlinear programming]: https://en.wikipedia.org/wiki/Nonlinear_programming

JuMP makes it easy to specify and **solve optimization problems without expert knowledge**, yet at the same time allows experts to implement advanced algorithmic techniques such as exploiting efficient hot-starts in linear programming or using callbacks to interact with branch-and-bound solvers. JuMP is also **fast** - benchmarking has shown that it can create problems at similar speeds to special-purpose commercial tools such as AMPL while maintaining the expressiveness of a generic high-level programming language. JuMP can be easily embedded in complex work flows including simulations and web servers.

Our documentation includes an installation guide, quick-start guide, and reference manual. The **[JuMPTutorials.jl]** repository contains a small but growing collection of contributed examples. Submissions are welcome!

[JuMPTutorials.jl]: https://github.com/jump-dev/JuMPTutorials.jl

**See [NEWS](https://github.com/jump-dev/JuMP.jl/tree/master/NEWS.md) for
a list of the significant breaking changes in the JuMP 0.19 release.**

**Latest Release**: 0.21.5 (`release-0.21` branch)
  * [Documentation](https://jump.dev/JuMP.jl/v0.21.5/)
  * Testing status:
    * Github Actions: [![Build Status](https://github.com/jump-dev/JuMP.jl/workflows/CI/badge.svg?branch=release-0.21)](https://github.com/jump-dev/JuMP.jl/actions?query=workflow%3ACI)


**Development version** (`master` branch):
  * [Documentation](https://jump.dev/JuMP.jl/dev/)
  * Testing status:
    * Github Actions: [![Build Status](https://github.com/jump-dev/JuMP.jl/workflows/CI/badge.svg?branch=master)](https://github.com/jump-dev/JuMP.jl/actions?query=workflow%3ACI)
    * Test coverage: [![codecov](https://codecov.io/gh/jump-dev/JuMP.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/jump-dev/JuMP.jl)
  * Changes: see [NEWS](https://github.com/jump-dev/JuMP.jl/tree/master/NEWS.md)
  * [Developer chatroom](https://gitter.im/JuliaOpt/JuMP-dev)


## Installation

JuMP can be installed through the Julia package manager:

```julia
julia> Pkg.add("JuMP")
```

For full installation instructions, including how to install solvers, see the documentation linked above.


## Supported problem classes

Mathematical optimization encompasses a large variety of problem classes.
We list below what is currently supported. See the documentation for more information.

**Objective types**

* Linear
* Convex Quadratic
* Nonlinear (convex and nonconvex)

**Constraint types**

* Linear
* Convex Quadratic
* Second-order Conic
* Semidefinite
* Nonlinear (convex and nonconvex)

**Variable types**

* Continuous
* Integer-valued
* Semicontinuous
* Semi-integer


## Bug reports and support

Please report any issues via the Github **[issue tracker]**. All types of issues are welcome and encouraged; this includes bug reports, documentation typos, feature requests, etc. The **[Optimization (Mathematical)]** category on Discourse is appropriate for general discussion, including "how do I do this?" questions.

[issue tracker]: https://github.com/jump-dev/JuMP.jl/issues
[Optimization (Mathematical)]: https://discourse.julialang.org/c/domain/opt


## Citing JuMP

If you find JuMP useful in your work, we kindly request that you cite the following paper ([pdf](https://mlubin.github.io/pdf/jump-sirev.pdf)):

```bibtex
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

For an earlier work where we presented a prototype implementation of JuMP, see [here](https://dx.doi.org/10.1287/ijoc.2014.0623):

```bibtex
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

![NumFOCUS logo](https://jump.dev/JuMP.jl/dev/assets/numfocus-logo.png)

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

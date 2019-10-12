![JuMP logo](https://www.juliaopt.org/images/jump-logo-with-text.svg "JuMP logo")
---

[![Powered by NumFOCUS](https://img.shields.io/badge/powered%20by-NumFOCUS-orange.svg?style=flat&colorA=E1523D&colorB=007D8A)](http://numfocus.org)

JuMP is a domain-specific modeling language for **[mathematical optimization]**
embedded in **[Julia]**. It currently supports a number of open-source and
commercial solvers ([Artelys Knitro], [BARON], [Bonmin], [Cbc], [CDCS], [CDD],
[Clp], [COSMO], [Couenne], [CPLEX], [CSDP], [ECOS], [FICO Xpress], [GLPK],
[Gurobi], [Ipopt], [Juniper], [MOSEK], [NLopt], [OSQP], [ProxSDP], [SCIP],
[SCS], [SDPA], [SDPT3], [SeDuMi], [Tulip]) for a variety of problem classes, including
**[linear programming]**, **[(mixed) integer programming]**,
**[second-order conic programming]**, **[semidefinite programming]**, and **[nonlinear programming]**.

[mathematical optimization]: http://en.wikipedia.org/wiki/Mathematical_optimization
[Julia]: http://julialang.org/
[Artelys Knitro]: http://artelys.com/en/optimization-tools/knitro
[BARON]: http://archimedes.cheme.cmu.edu/?q=baron
[Bonmin]: https://projects.coin-or.org/Bonmin
[Cbc]: https://projects.coin-or.org/Cbc
[CDCS]: https://github.com/oxfordcontrol/CDCS
[CDD]: https://github.com/cddlib/cddlib
[Clp]: https://projects.coin-or.org/Clp
[COSMO]: https://github.com/oxfordcontrol/COSMO.jl
[Couenne]: https://projects.coin-or.org/Couenne
[CPLEX]: http://www-01.ibm.com/software/commerce/optimization/cplex-optimizer/
[CSDP]: https://projects.coin-or.org/Csdp/
[ECOS]: https://github.com/ifa-ethz/ecos
[FICO Xpress]: http://www.fico.com/en/products/fico-xpress-optimization-suite
[GLPK]: http://www.gnu.org/software/glpk/
[Gurobi]: http://www.gurobi.com/
[Ipopt]: https://projects.coin-or.org/Ipopt
[Juniper]: https://github.com/lanl-ansi/Juniper.jl
[MOSEK]: http://mosek.com/
[NLopt]: http://ab-initio.mit.edu/wiki/index.php/NLopt
[OSQP]: https://osqp.org/
[ProxSDP]: https://github.com/mariohsouto/ProxSDP.jl
[SCIP]: https://scip.zib.de/
[SCS]: https://github.com/cvxgrp/scs
[SDPA]: http://sdpa.sourceforge.net/
[SDPT3]: https://blog.nus.edu.sg/mattohkc/softwares/sdpt3/
[SeDuMi]: http://sedumi.ie.lehigh.edu/
[Tulip]: https://github.com/ds4dm/Tulip.jl
[linear programming]: http://en.wikipedia.org/wiki/Linear_programming
[(mixed) integer programming]: http://en.wikipedia.org/wiki/Integer_programming
[second-order conic programming]: http://en.wikipedia.org/wiki/Second-order_cone_programming
[semidefinite programming]: https://en.wikipedia.org/wiki/Semidefinite_programming
[nonlinear programming]: http://en.wikipedia.org/wiki/Nonlinear_programming

JuMP makes it easy to specify and **solve optimization problems without expert knowledge**, yet at the same time allows experts to implement advanced algorithmic techniques such as exploiting efficient hot-starts in linear programming or using callbacks to interact with branch-and-bound solvers. JuMP is also **fast** - benchmarking has shown that it can create problems at similar speeds to special-purpose commercial tools such as AMPL while maintaining the expressiveness of a generic high-level programming language. JuMP can be easily embedded in complex work flows including simulations and web servers.

Our documentation includes an installation guide, quick-start guide, and reference manual. The **[juliaopt-notebooks]** repository contains a small but growing collection of contributed examples. Submissions are welcome!

[juliaopt-notebooks]: https://github.com/JuliaOpt/juliaopt-notebooks

**See [NEWS](https://github.com/JuliaOpt/JuMP.jl/tree/master/NEWS.md) for
a list of the significant breaking changes in the JuMP 0.19 release.**

**Latest Release**: 0.20.0 (`release-0.20` branch)
  * [Documentation](http://www.juliaopt.org/JuMP.jl/v0.20.0/)
  * [Examples](https://github.com/JuliaOpt/JuMP.jl/tree/release-0.20/examples)
  * Testing status:
    * TravisCI: [![Build Status](https://travis-ci.org/JuliaOpt/JuMP.jl.svg?branch=release-0.20)](https://travis-ci.org/JuliaOpt/JuMP.jl)


**Development version** (`master` branch):
  * [Documentation](http://www.juliaopt.org/JuMP.jl/dev/)
  * [Examples](https://github.com/JuliaOpt/JuMP.jl/tree/master/examples)
  * Testing status:
    * TravisCI: [![Build Status](https://travis-ci.org/JuliaOpt/JuMP.jl.svg?branch=master)](https://travis-ci.org/JuliaOpt/JuMP.jl)
    * Test coverage:
      [![Coverage Status](https://coveralls.io/repos/JuliaOpt/JuMP.jl/badge.svg?branch=master)](https://coveralls.io/r/JuliaOpt/JuMP.jl?branch=master)
      [![codecov](https://codecov.io/gh/JuliaOpt/JuMP.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/JuliaOpt/JuMP.jl)
  * Changes: see [NEWS](https://github.com/JuliaOpt/JuMP.jl/tree/master/NEWS.md)
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

[issue tracker]: https://github.com/JuliaOpt/JuMP.jl/issues
[Optimization (Mathematical)]: https://discourse.julialang.org/c/domain/opt


## Citing JuMP

If you find JuMP useful in your work, we kindly request that you cite the following paper ([pdf](https://mlubin.github.io/pdf/jump-sirev.pdf)):

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

For an earlier work where we presented a prototype implementation of JuMP, see [here](http://dx.doi.org/10.1287/ijoc.2014.0623):

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

A preprint of this paper is [freely available](http://arxiv.org/abs/1312.1431).

---

![NumFOCUS logo](http://www.juliaopt.org/JuMP.jl/dev/assets/numfocus-logo.png)

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

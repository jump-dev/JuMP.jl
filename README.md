JuMP
====
#### Julia for Mathematical Programming

JuMP is a domain-specific modeling language for **[mathematical programming]**
embedded in **[Julia]**. It currently supports a number of open-source and
commercial solvers ([CPLEX], [COIN Clp], [COIN Cbc], [ECOS], [GLPK],
[Gurobi], [Ipopt], [KNITRO], [MOSEK], and [NLopt]) for a variety of problem classes, including
**[linear programming]**, **[(mixed) integer programming]**,
**[second-order conic programming]**, and **[nonlinear programming]**.

[mathematical programming]: http://en.wikipedia.org/wiki/Mathematical_optimization
[Julia]: http://julialang.org/
[COIN Clp]: https://projects.coin-or.org/Clp
[COIN Cbc]: https://projects.coin-or.org/Cbc
[ECOS]: https://github.com/ifa-ethz/ecos
[GLPK]: http://www.gnu.org/software/glpk/
[Gurobi]: http://www.gurobi.com/
[MOSEK]: http://mosek.com/
[CPLEX]: http://www-01.ibm.com/software/commerce/optimization/cplex-optimizer/
[Ipopt]: https://projects.coin-or.org/Ipopt
[KNITRO]: http://www.ziena.com/knitro.htm
[NLopt]: http://ab-initio.mit.edu/wiki/index.php/NLopt
[linear programming]: http://en.wikipedia.org/wiki/Linear_programming
[(mixed) integer programming]: http://en.wikipedia.org/wiki/Integer_programming
[second-order conic programming]: http://en.wikipedia.org/wiki/Second-order_cone_programming
[nonlinear programming]: http://en.wikipedia.org/wiki/Nonlinear_programming

JuMP makes it easy to specify and **solve optimization problems without expert knowledge**, yet at the same time allows experts to implement advanced algorithmic techniques such as exploiting efficient hot-starts in linear programming or using callbacks to interact with branch-and-bound solvers. JuMP is also **fast** - benchmarking has shown that it can create problems at similar speeds to special-purpose commercial tools such as AMPL while maintaining the expressiveness of a generic high-level programming language. JuMP can be easily embedded in complex work flows including simulations and web servers.

Our documentation includes an installation guide, quick-start guide, and reference manual.

**Latest Release**: 0.7.3 (via ``Pkg.add``)
  * [documentation](https://jump.readthedocs.org/en/release-0.7)
  * [examples](https://github.com/JuliaOpt/JuMP.jl/tree/release-0.7/examples)
  * Testing status: [![Build Status](https://travis-ci.org/JuliaOpt/JuMP.jl.png?branch=release-0.7)](https://travis-ci.org/JuliaOpt/JuMP.jl) [![JuMP](http://pkg.julialang.org/badges/JuMP_release.svg)](http://pkg.julialang.org/?pkg=JuMP&ver=release)


**Development version**:
  * [documentation](https://jump.readthedocs.org/en/latest)
  * [examples](https://github.com/JuliaOpt/JuMP.jl/tree/master/examples)
  * Testing status: [![Build Status](https://travis-ci.org/JuliaOpt/JuMP.jl.png?branch=master)](https://travis-ci.org/JuliaOpt/JuMP.jl) [![Coverage Status](https://coveralls.io/repos/JuliaOpt/JuMP.jl/badge.png)](https://coveralls.io/r/JuliaOpt/JuMP.jl)
  * Changes: see [NEWS](https://github.com/JuliaOpt/JuMP.jl/tree/master/NEWS.md)

## Installation

JuMP can be installed through the Julia package manager (version 0.3 required)

```julia
julia> Pkg.add("JuMP")
```

For full installation instructions, including how to install solvers, see the documentation linked above.



## Supported problem classes

Mathematical programming encompasses a large variety of problem classes.
We list below what is currently supported. See the documentation for more information.

**Objective types**

* Linear
* Convex Quadratic
* Nonlinear (convex and nonconvex)

**Constraint types**

* Linear
* Convex Quadratic
* Second-order Conic
* Nonlinear (convex and nonconvex)

**Variable types**

* Continuous
* Integer-valued
* Semicontinuous
* Semi-integer

## Bug reports and support

Please report any issues via the Github **[issue tracker]**. All types of issues are welcome and encouraged; this includes bug reports, documentation typos, feature requests, etc. The **[julia-opt]** mailing list is appropriate for general discussion, including "how do I do this?" questions.


[issue tracker]: https://github.com/JuliaOpt/JuMP.jl/issues
[julia-opt]: https://groups.google.com/forum/#!forum/julia-opt

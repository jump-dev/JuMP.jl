JuMP
====
#### Julia for Mathematical Programming

JuMP is a domain-specific modeling language for **[mathematical programming]**
embedded in **[Julia]**. It currently supports a number of open-source and
commercial solvers ([COIN Clp], [COIN Cbc], [GNU GLPK], [Gurobi], [MOSEK], [CPLEX]) via a 
[generic solver-independent interface](https://github.com/JuliaOpt/MathProgBase.jl). 

One the best features of JuMP is its **speed** - benchmarking has shown that it
can create problems at similar speeds to special-purpose modeling languages
such as AMPL while maintaining the expressiveness of a generic high-level 
programming language. JuMP communicates with solvers in-memory, 
avoiding the need to write intermediary files and enabling access to **advanced
features** such as efficient LP re-solves and callbacks for mixed-integer programming.

JuMP has recently enabled support for nonlinear programming for functions that can be expressed in closed algebraic form. JuMP computes exact sparse second-order derivatives needed by efficient interior-point solvers. Currently, [Ipopt] is the only supported solver, with more coming soon. 

Our documentation includes an installation guide, quick-start guide, and reference manual. 

**Latest Release**: 0.5.1 (via ``Pkg.add``)
  * [documentation](https://jump.readthedocs.org/en/release-0.5)
  * [examples](https://github.com/JuliaOpt/JuMP.jl/tree/release-0.5/examples)
  * Testing status: [![Build Status](https://travis-ci.org/JuliaOpt/JuMP.jl.png?branch=release-0.5)](https://travis-ci.org/JuliaOpt/JuMP.jl) [![Build status](https://ci.appveyor.com/api/projects/status/val81xkp6y6uiw8g/branch/release-0.5)](https://ci.appveyor.com/project/mlubin/jump-jl) [![PackageEvaluator.jl](http://iainnz.github.io/packages.julialang.org/badges/JuMP_0.2.svg)](http://iainnz.github.io/packages.julialang.org/?pkg=JuMP&ver=0.2)


**Development version**: 
  * [documentation](https://jump.readthedocs.org/en/latest)
  * [examples](https://github.com/JuliaOpt/JuMP.jl/tree/master/examples) 
  * Testing status: [![Build Status](https://travis-ci.org/JuliaOpt/JuMP.jl.png?branch=master)](https://travis-ci.org/JuliaOpt/JuMP.jl) [![Build status](https://ci.appveyor.com/api/projects/status/val81xkp6y6uiw8g/branch/master)](https://ci.appveyor.com/project/mlubin/jump-jl) [![PackageEvaluator.jl](http://iainnz.github.io/packages.julialang.org/badges/JuMP_0.3.svg)](http://iainnz.github.io/packages.julialang.org/?pkg=JuMP&ver=0.3)
  * Changes: see [NEWS](https://github.com/JuliaOpt/JuMP.jl/tree/master/NEWS.md)

*JuMP was formerly known as MathProg.jl*

## Installation

JuMP can be installed through the Julia package manager (version 0.2 required)

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

## Bug reports and support

Please report any issues via the Github **[issue tracker]**. All types of issues are welcome and encouraged; this includes bug reports, documentation typos, feature requests, etc. The **[julia-opt]** mailing list is appropriate for general discussion, including "how do I do this?" questions.


[issue tracker]: https://github.com/JuliaOpt/JuMP.jl/issues
[mathematical programming]: http://en.wikipedia.org/wiki/Mathematical_optimization
[Julia]: http://julialang.org/
[COIN Clp]: https://github.com/mlubin/Clp.jl
[COIN Cbc]: https://github.com/mlubin/Cbc.jl
[GNU GLPK]: http://www.gnu.org/software/glpk/
[Gurobi]: http://www.gurobi.com/
[MOSEK]: http://mosek.com/
[CPLEX]: http://www-01.ibm.com/software/commerce/optimization/cplex-optimizer/
[Ipopt]: https://projects.coin-or.org/Ipopt
[julia-opt]: https://groups.google.com/forum/#!forum/julia-opt

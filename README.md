JuMP
====
#### Julia for Mathematical Programming

JuMP is a domain-specific modeling language for **[mathematical programming]**
embedded in **[Julia]**. It currently supports a number of open-source and
commercial solvers ([COIN Clp], [COIN Cbc], [GNU GLPK], [Gurobi]) via a 
[generic solver-independent interface](https://github.com/JuliaOpt/MathProgBase.jl). 
One the best features of JuMP is its speed - benchmarking has shown that it
can create problems at similar speeds to special-purpose modeling languages
such as AMPL while maintaing the expressiveness of a generic high-level 
programming language.

Our **documentation**, including installation guide, quick-start guide, and
reference material is available [here](https://jump.readthedocs.org/en/latest/jump.html#installation-guide).

We have a comprehensive test suite: the current status of the development version of
JuMP is [![Build Status](https://travis-ci.org/JuliaOpt/JuMP.jl.png?branch=master)](https://travis-ci.org/JuliaOpt/JuMP.jl)

*JuMP was formerly known as MathProg.jl*

## Installation

JuMP can be installed through the Julia package manager (version 0.2 required)

```julia
julia> Pkg.add("JuMP")
```

For full installation instructions, including how to install solvers, see the [installation guide](https://jump.readthedocs.org/en/latest/jump.html#installation-guide)

[mathematical programming]: http://en.wikipedia.org/wiki/Mathematical_optimization
[Julia]: http://julialang.org/
[COIN Clp]: https://github.com/mlubin/Clp.jl
[COIN Cbc]: https://github.com/mlubin/Cbc.jl
[GNU GLPK]: http://www.gnu.org/software/glpk/
[Gurobi]: http://www.gurobi.com/

## Supported problem classes

Mathematical programming encompasses a large variety of problem classes. 
We list below what is currently supported. See the documentation for more information. 

**Objective types**

* Linear
* Convex Quadratic (Gurobi only)

**Constraint types**

* Linear
* Convex Quadratic (Gurobi only)
* Second-order Conic (Gurobi only)

**Variable types**

* Continous
* Integer-valued

## Bug reports and support

Please report any issues via the Github **[issue tracker]**. All types of issues are welcome and encouraged; this includes bug reports, documentation typos, "how do I do this?" questions, feature requests, etc.


[issue tracker]: https://github.com/JuliaOpt/JuMP.jl/issues

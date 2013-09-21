JuMP
====
#### Julia for Mathematical Programming

JuMP is a domain-specific modeling language
for **[mathematical programming]** 
embedded in **[Julia]**. 
It it currently supports a number of open-source and commerical 
solvers (COIN Clp, COIN Cbc, GNU GLPK, and Gurobi) via a generic solver-independent
interface provided by the **[MathProgBase]** package. 
One the best features of JuMP is its speed -
benchmarking has shown that it can create problems at similar
speeds to special-purpose modeling languages such as AMPL
while maintaing the expressiveness of a generic high-level programming language.

JuMP was formerly known as MathProg.jl

[![Build Status](https://travis-ci.org/IainNZ/JuMP.jl.png)](https://travis-ci.org/IainNZ/JuMP.jl)

# Installation

JuMP can be installed through the Julia package manager (version 0.2 required)

```
    julia> Pkg.add("JuMP")
```

This will install JuMP's default solvers: **[Clp]** and **[Cbc]**. All Julia-supported platforms (including Linux, OS X, and Windows) are supported by JuMP. JuMP itself is written in pure Julia; however, external solvers introduce binary dependencies. See the README for **[Cbc]** for platform-specifc installation details.

[MathProgBase]: https://github.com/mlubin/MathProgBase.jl
[Clp]: https://github.com/mlubin/Clp.jl
[Cbc]: https://github.com/mlubin/Cbc.jl
[matematical programming]: http://en.wikipedia.org/wiki/Mathematical_optimization
[Julia]: http://julialang.org/

# Supported problem classes

Mathematical programming encompasses a large variety of problem classes. 
We list below what is currently supported. See the documentation for more information. 

## Objective types 

* Linear
* Convex Quadratic (Gurobi only)

## Constraint types

* Linear
* Convex Quadratic (Gurobi only)
* Second-order Conic (Gurobi only)

## Variable types

* Continous
* Integer-valued

# Documentation

JuMP documentation, including a **[quick-start guide]** is available at https://jump.readthedocs.org/en/latest/jump.html or in the doc/ directory. 

[quick-start guide]: https://jump.readthedocs.org/en/latest/jump.html#quick-start-guide

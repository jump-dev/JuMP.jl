Installation Guide
==================

DRAFT: Install Julia. Install JuMP.

Getting Solvers
---------------

TODO: Not yet updated for MOI.

Solver support in Julia is currently provided by writing a solver-specific
package that provides a very thin wrapper around the solver's C interface and
providing a standard interface that JuMP can call. If you are interested in
providing an interface to your solver, please get in touch. The table below
lists the currently supported solvers and their capabilities.


| Solver                                                                           | Julia Package                                                                   | `solver=`                 | License     | LP    | SOCP  | MILP  | NLP   | MINLP | SDP   |
| -------------------------------------------------------------------------------- | ------------------------------------------------------------------------------- | ------------------------- | ----------- | ----- | ----- | ----- | ----- | ----- | ----- |
| [Artelys Knitro](http://artelys.com/en/optimization-tools/knitro)                | [KNITRO.jl](https://github.com/JuliaOpt/KNITRO.jl)                              | `KnitroSolver()`          |  Comm.      |       |       |       | X     | X     |       |
| [BARON](http://archimedes.cheme.cmu.edu/?q=baron>)                               | [BARON.jl](https://github.com/joehuchette/BARON.jl)                             | `BaronSolver()`           |  Comm.      |       |       |       | X     | X     |       |
| [Bonmin](https://projects.coin-or.org/Bonmin)                                    | [AmplNLWriter.jl](https://github.com/JuliaOpt/AmplNLWriter.jl)                | `AmplNLWriter(CoinOptServices.bonmin)` *      |  EPL        | X     |       | X     | X     | X     |       |
| ''                                                                               | [CoinOptServices.jl](https://github.com/JuliaOpt/CoinOptServices.jl)            | `OsilBonminSolver()`      |  ''         |       |       |       |       |       |       |
| [Cbc](https://projects.coin-or.org/Cbc)                                          | [Cbc.jl](https://github.com/JuliaOpt/Cbc.jl)                                    | `CbcSolver()`             |  EPL        |       |       | X     |       |       |       |
| [Clp](https://projects.coin-or.org/Clp)                                          | [Clp.jl](https://github.com/JuliaOpt/Clp.jl)                                    | `ClpSolver()`             |  EPL        | X     |       |       |       |       |       |
| [Couenne](https://projects.coin-or.org/Couenne)                                  | [AmplNLWriter.jl](https://github.com/JuliaOpt/AmplNLWriter.jl)                | `AmplNLWriter(CoinOptServices.couenne)` *     |  EPL        | X     |       | X     | X     | X     |       |
| ''                                                                               | [CoinOptServices.jl](https://github.com/JuliaOpt/CoinOptServices.jl)            | `OsilCouenneSolver()`     |  ''         |       |       |       |       |       |       |
| [CPLEX](http://www-01.ibm.com/software/commerce/optimization/cplex-optimizer/)   | [CPLEX.jl](https://github.com/JuliaOpt/CPLEX.jl)                                | `CplexSolver()`           |  Comm.      | X     | X     | X     |       |       |       |
| [ECOS](https://github.com/ifa-ethz/ecos)                                         | [ECOS.jl](https://github.com/JuliaOpt/ECOS.jl)                                  | `ECOSSolver()`            |  GPL        | X     | X     |       |       |       |       |
| [FICO Xpress](http://www.fico.com/en/products/fico-xpress-optimization-suite)    | [Xpress.jl](https://github.com/JuliaOpt/Xpress.jl)                              | `XpressSolver()`          |  Comm.      | X     | X     | X     |       |       |       |
| [GLPK](http://www.gnu.org/software/glpk/)                                        | [GLPKMath...](https://github.com/JuliaOpt/GLPKMathProgInterface.jl)             | `GLPKSolver[LP\|MIP]()`   |  GPL        | X     |       | X     |       |       |       |
| [Gurobi](http://gurobi.com)                                                      | [Gurobi.jl](https://github.com/JuliaOpt/Gurobi.jl)                              | `GurobiSolver()`          |  Comm.      | X     | X     | X     |       |       |       |
| [Ipopt](https://projects.coin-or.org/Ipopt)                                      | [Ipopt.jl](https://github.com/JuliaOpt/Ipopt.jl)                                | `IpoptSolver()`           |  EPL        | X     |       |       | X     |       |       |
| [MOSEK](http://www.mosek.com/)                                                   | [Mosek.jl](https://github.com/JuliaOpt/Mosek.jl)                                | `MosekSolver()`           |  Comm.      | X     | X     | X     | X     |       | X     |
| [NLopt](http://ab-initio.mit.edu/wiki/index.php/NLopt)                           | [NLopt.jl](https://github.com/JuliaOpt/NLopt.jl)                                | `NLoptSolver()`           |  LGPL       |       |       |       | X     |       |       |
| [SCS](https://github.com/cvxgrp/scs>)                                            | [SCS.jl](https://github.com/JuliaOpt/SCS.jl)                                    | `SCSSolver()`             |  MIT        | X     | X     |       |       |       | X     |


Where:

-   LP = Linear programming
-   SOCP = Second-order conic programming (including problems with convex quadratic constraints and/or objective)
-   MILP = Mixed-integer linear programming
-   NLP = Nonlinear programming
-   MINLP = Mixed-integer nonlinear programming
-   SDP = Semidefinite programming

\* requires CoinOptServices installed, see below.

To install Gurobi, for example, and use it with a JuMP model `m`, run:

```julia
Pkg.add("Gurobi")
using JuMP
using Gurobi

m = Model(solver=GurobiSolver())
```

Setting solver options is discussed in the Model &lt;ref-model&gt; section.

Solver-specific notes follow below.

### Artelys Knitro

Requires a license. The KNITRO.jl interface currently supports only nonlinear problems.

### BARON

Requires a license. A trial version is available for small problem instances.

### COIN-OR Clp and Cbc

Binaries for Clp and Cbc are provided on OS X and Windows (32- and 64-bit) by default. On Linux, they will be compiled from source (be sure to have a C++ compiler installed). Cbc supports "SOS" constraints but does *not* support MIP callbacks.

### CPLEX

Requires a working installation of CPLEX with a license (free for faculty members and graduate teaching assistants). The interface requires using CPLEX as a shared library, which is unsupported by the CPLEX developers. Special installation steps are required on OS X. CPLEX supports MIP callbacks and "SOS" constraints.

### ECOS

ECOS can be used by JuMP to solve LPs and SOCPs. ECOS does not support general quadratic objectives or constraints, only second-order conic constraints specified by using `norm` or the quadratic form `x'x <= y^2`.

### FICO Xpress

Requires a working installation of Xpress with an active license (it is possible to get license for academic use, see [FICO Academic Partner Program](http://subscribe.fico.com/Academic-Partner-Program)). Supports SOCP and "SOS" constraints. The interface is experimental, but it does pass all JuMP and MathProgBase tests. Callbacks are not yet supported.

!!! warning
    If you are using 64-bit Xpress, you must use 64-bit Julia (and similarly with 32-bit Xpress).

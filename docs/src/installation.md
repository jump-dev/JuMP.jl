Installation Guide
==================

JuMP is a package for [Julia](https://julialang.org). To use JuMP, first
[download and install](https://julialang.org/downloads/) Julia or open up
a remote notebook at [JuliaBox](https://www.juliabox.com/) or similar services.

This version of JuMP is compatible with Julia 0.7 and later.
The following instructions are based on Julia 1.0.

From Julia, JuMP is installed by using the built-in package manager:
```julia
import Pkg
Pkg.add("JuMP")
```

!!! note
    The installation instructions above assume that JuMP 0.19 has already been
    released. Until that time, see the JuMP
    [README](https://github.com/JuliaOpt/JuMP.jl/blob/master/README.md) for
    instructions on installing a development release that's compatible with
    Julia 1.0.


Getting Solvers
---------------

JuMP depends on solvers to solve optimization problems. Most solvers are not
written in Julia, and some require commercial licenses to use, so installation
is often more complex. We list below the currently available solvers.

!!! note
    This list is open for new contributions. See also
    [Interacting with solvers](@ref) and the
    [MathOptInterface docs](http://www.juliaopt.org/MathOptInterface.jl/v0.6.1/)
    for more details on how JuMP interacts with solvers. Please get in touch
    with any questions about connecting new solvers with JuMP.


| Solver                                                                         | Julia Package                                                                    | License | Supports                    |
| ------------------------------------------------------------------------------ | -------------------------------------------------------------------------------- | ------- | --------------------------- |
| [Cbc](https://projects.coin-or.org/Cbc)                                        | [Cbc.jl](https://github.com/JuliaOpt/Cbc.jl)                                     | EPL     | MILP                        |
| [Clp](https://projects.coin-or.org/Clp)                                        | [Clp.jl](https://github.com/JuliaOpt/Clp.jl)                                     | EPL     | LP                          |
| [CPLEX](http://www-01.ibm.com/software/commerce/optimization/cplex-optimizer/) | [CPLEX.jl](https://github.com/JuliaOpt/CPLEX.jl)                                 | Comm.   | LP, MILP, SOCP, MISOCP      |
| [CSDP](https://projects.coin-or.org/Csdp/)                                     | [CSDP.jl](https://github.com/JuliaOpt/CSDP.jl)                                   | EPL     | LP, SDP                     |
| [ECOS](https://github.com/ifa-ethz/ecos)                                       | [ECOS.jl](https://github.com/JuliaOpt/ECOS.jl)                                   | GPL     | LP, SOCP                    |
| [FICO Xpress](http://www.fico.com/en/products/fico-xpress-optimization-suite)  | [Xpress.jl](https://github.com/JuliaOpt/Xpress.jl)                               | Comm.   | LP, MILP, SOCP, MISOCP      |
| [GLPK](http://www.gnu.org/software/glpk/)                                      | [GLPK.jl](https://github.com/JuliaOpt/GLPK.jl)                                   | GPL     | LP, MILP                    |
| [Gurobi](http://gurobi.com)                                                    | [Gurobi.jl](https://github.com/JuliaOpt/Gurobi.jl)                               | Comm.   | LP, MILP, SOCP, MISOCP      |
| [Ipopt](https://projects.coin-or.org/Ipopt)                                    | [Ipopt.jl](https://github.com/JuliaOpt/Ipopt.jl)                                 | EPL     | LP, QP, NLP                 |
| [MOSEK](http://www.mosek.com/)                                                 | [MathOptInterfaceMosek.jl](https://github.com/JuliaOpt/MathOptInterfaceMosek.jl) | Comm.   | LP, MILP, SOCP, MISOCP, SDP |
| [OSQP](https://osqp.org/)                                                      | [OSQP.jl](https://github.com/oxfordcontrol/OSQP.jl)                              | Apache  | LP, QP                      |
| [SCS](https://github.com/cvxgrp/scs>)                                          | [SCS.jl](https://github.com/JuliaOpt/SCS.jl)                                     | MIT     | LP, SOCP, SDP               |


Where:

-   LP = Linear programming
-   QP = Quadratic programming
-   SOCP = Second-order conic programming (including problems with convex quadratic constraints and/or objective)
-   MILP = Mixed-integer linear programming
-   NLP = Nonlinear programming
-   MINLP = Mixed-integer nonlinear programming
-   SDP = Semidefinite programming
-   MISDP = Mixed-integer semidefinite programming

To install Gurobi, for example, and use it with a JuMP model `model`, run:

```julia
import Pkg
Pkg.add("Gurobi")
using JuMP
using Gurobi
model = Model(with_optimizer(Gurobi.Optimizer))
```

Most packages follow the `ModuleName.Optimizer` naming convention, but
exceptions may exist. See the corresponding Julia package README for more
details on how to use the solver.

TODO: Discuss setting solver options.

The following solvers were compatible with JuMP up to release 0.18 but are
not yet compatible with the latest version because they do not implement the
new MathOptInterface API:

- [Artelys Knitro](https://github.com/JuliaOpt/KNITRO.jl)
- [BARON](https://github.com/joehuchette/BARON.jl)
- [Bonmin](https://projects.coin-or.org/Bonmin) and
  [Couenne](https://projects.coin-or.org/Couenne) via
  [AmplNLWriter.jl](https://github.com/JuliaOpt/AmplNLWriter.jl)
- [CDD](https://github.com/JuliaPolyhedra/CDDLib.jl)
- [NLopt](https://github.com/JuliaOpt/NLopt.jl)
- [Pavito](https://github.com/JuliaOpt/Pavito.jl)
- [Pajarito](https://github.com/JuliaOpt/Pajarito.jl)
- [SCIP](https://github.com/SCIP-Interfaces/SCIP.jl)

Solver-specific notes follow below.

### Artelys Knitro

Requires a license. The KNITRO.jl interface currently supports only nonlinear problems.

### BARON

Requires a license. A trial version is available for small problem instances.

### COIN-OR Cbc

Cbc supports "SOS" constraints.

### CPLEX

Requires a working installation of CPLEX with a license (free for faculty
members and graduate teaching assistants). The interface requires using CPLEX as
a shared library, which is unsupported by the CPLEX developers. Special
installation steps are required on Mac OS. CPLEX supports "SOS" constraints.

### ECOS

ECOS can be used by JuMP to solve LPs and SOCPs. ECOS does not support general
quadratic objectives or constraints, only second-order conic constraints
specified by using the `SecondOrderCone` set.

### Gurobi

Requires a working installation of Gurobi with an activated license (free for
academic use). Gurobi supports "SOS" constraints.

### FICO Xpress

Requires a working installation of Xpress with an active license (it is possible
to get a license for academic use, see
[FICO Academic Partner Program](http://subscribe.fico.com/Academic-Partner-Program)).
Supports SOCP and "SOS" constraints.

### MOSEK

Requires a license (free for academic use). The Mosek interface is maintained by
the Mosek team. (Thanks!)

### SCS

SCS can be used by JuMP to solve LPs and SOCPs, and SDPs. SCS is a first order
solver and has low accuracy (``10^{−4}``) by default; see the SCS.jl
documentation for more information.

Installation Guide
==================

JuMP is a package for [Julia](https://julialang.org). To use JuMP, first
[download and install](https://julialang.org/downloads/) Julia or open up
a remote notebook at [JuliaBox](https://www.juliabox.com/) or similar services.

This version of JuMP is compatible with Julia 1.0 and later.

From Julia, JuMP is installed by using the built-in package manager:
```julia
import Pkg
Pkg.add("JuMP")
```

Getting Solvers
---------------

JuMP depends on solvers to solve optimization problems. Most solvers are not
written in Julia, and some require commercial licenses to use, so installation
is often more complex. We list below the currently available solvers.

!!! note
    This list is open for new contributions. See also
    [Interacting with solvers](@ref) and the
    [MathOptInterface docs](http://www.juliaopt.org/MathOptInterface.jl/v0.9.1/)
    for more details on how JuMP interacts with solvers. Please get in touch
    with any questions about connecting new solvers with JuMP.


| Solver                                                                         | Julia Package                                                                    | License | Supports                           |
| ------------------------------------------------------------------------------ | -------------------------------------------------------------------------------- | ------- | ---------------------------------- |
| [Artelys Knitro](https://www.artelys.com/knitro)                               | [KNITRO.jl](https://github.com/JuliaOpt/KNITRO.jl)                               | Comm.   | LP, MILP, SOCP, MISOCP, NLP, MINLP |
| [Cbc](https://projects.coin-or.org/Cbc)                                        | [Cbc.jl](https://github.com/JuliaOpt/Cbc.jl)                                     | EPL     | MILP                               |
| [CDCS](https://github.com/oxfordcontrol/CDCS)                                  | [CDCS.jl](https://github.com/oxfordcontrol/CDCS.jl)                              | GPL     | LP, SOCP, SDP                      |
| [CDD](https://github.com/cddlib/cddlib)                                        | [CDDLib.jl](https://github.com/JuliaPolyhedra/CDDLib.jl)                         | GPL     | LP                                 |
| [Clp](https://projects.coin-or.org/Clp)                                        | [Clp.jl](https://github.com/JuliaOpt/Clp.jl)                                     | EPL     | LP                                 |
| [COSMO](https://github.com/oxfordcontrol/COSMO.jl)                             | [COSMO.jl](https://github.com/oxfordcontrol/COSMO.jl)                            | Apache  | LP, QP, SOCP, SDP                  |
| [CPLEX](http://www-01.ibm.com/software/commerce/optimization/cplex-optimizer/) | [CPLEX.jl](https://github.com/JuliaOpt/CPLEX.jl)                                 | Comm.   | LP, MILP, SOCP, MISOCP             |
| [CSDP](https://projects.coin-or.org/Csdp/)                                     | [CSDP.jl](https://github.com/JuliaOpt/CSDP.jl)                                   | EPL     | LP, SDP                            |
| [ECOS](https://github.com/ifa-ethz/ecos)                                       | [ECOS.jl](https://github.com/JuliaOpt/ECOS.jl)                                   | GPL     | LP, SOCP                           |
| [FICO Xpress](http://www.fico.com/en/products/fico-xpress-optimization-suite)  | [Xpress.jl](https://github.com/JuliaOpt/Xpress.jl)                               | Comm.   | LP, MILP, SOCP, MISOCP             |
| [GLPK](http://www.gnu.org/software/glpk/)                                      | [GLPK.jl](https://github.com/JuliaOpt/GLPK.jl)                                   | GPL     | LP, MILP                           |
| [Gurobi](http://gurobi.com)                                                    | [Gurobi.jl](https://github.com/JuliaOpt/Gurobi.jl)                               | Comm.   | LP, MILP, SOCP, MISOCP             |
| [Ipopt](https://projects.coin-or.org/Ipopt)                                    | [Ipopt.jl](https://github.com/JuliaOpt/Ipopt.jl)                                 | EPL     | LP, QP, NLP                        |
| [Juniper](https://github.com/lanl-ansi/Juniper.jl)                             | [Juniper.jl](https://github.com/lanl-ansi/Juniper.jl)                            | MIT     | MISOCP, MINLP                      |
| [MOSEK](http://www.mosek.com/)                                                 | [MosekTools.jl](https://github.com/JuliaOpt/MosekTools.jl)                       | Comm.   | LP, MILP, SOCP, MISOCP, SDP        |
| [OSQP](https://osqp.org/)                                                      | [OSQP.jl](https://github.com/oxfordcontrol/OSQP.jl)                              | Apache  | LP, QP                             |
| [ProxSDP](https://github.com/mariohsouto/ProxSDP.jl)                           | [ProxSDP.jl](https://github.com/mariohsouto/ProxSDP.jl)                          | MIT     | LP, SOCP, SDP                      |
| [SCIP](https://scip.zib.de/)                                                   | [SCIP.jl](https://github.com/SCIP-Interfaces/SCIP.jl)                            | ZIB     | MILP, MINLP                        |
| [SCS](https://github.com/cvxgrp/scs)                                           | [SCS.jl](https://github.com/JuliaOpt/SCS.jl)                                     | MIT     | LP, SOCP, SDP                      |
| [SDPA](http://sdpa.sourceforge.net/)                                           | [SDPA.jl](https://github.com/JuliaOpt/SDPA.jl), [SDPAFamily.jl](https://github.com/ericphanson/SDPAFamily.jl)                                   | GPL     | LP, SDP                            |
| [SDPT3](https://blog.nus.edu.sg/mattohkc/softwares/sdpt3/)                     | [SDPT3.jl](https://github.com/JuliaOpt/SDPT3.jl)                                 | GPL     | LP, SOCP, SDP                      |
| [SeDuMi](http://sedumi.ie.lehigh.edu/)                                         | [SeDuMi.jl](https://github.com/JuliaOpt/SeDuMi.jl)                               | GPL     | LP, SOCP, SDP                      |
| [Tulip](https://github.com/ds4dm/Tulip.jl)                                     | [Tulip.jl](https://github.com/ds4dm/Tulip.jl)                                    | MPL-2   | LP                                 |

Where:

-   LP = Linear programming
-   QP = Quadratic programming
-   SOCP = Second-order conic programming (including problems with convex quadratic constraints and/or objective)
-   MILP = Mixed-integer linear programming
-   NLP = Nonlinear programming
-   MINLP = Mixed-integer nonlinear programming
-   SDP = Semidefinite programming
-   MISDP = Mixed-integer semidefinite programming

You may also use [AmplNLWriter](https://github.com/JuliaOpt/AmplNLWriter.jl) to
access solvers that support the [nl format](https://en.wikipedia.org/wiki/Nl_(format)).
Such solvers include [Bonmin](https://projects.coin-or.org/Bonmin) and
[Couenne](https://projects.coin-or.org/Couenne). See a more complete list
[here](https://ampl.com/products/solvers/all-solvers-for-ampl/).

To install Gurobi, for example, and use it with a JuMP model `model`, run:

```julia
import Pkg
Pkg.add("Gurobi")
using JuMP
using Gurobi
model = Model(Gurobi.Optimizer)
```

Most packages follow the `ModuleName.Optimizer` naming convention, but
exceptions may exist. See the corresponding Julia package README for more
details on how to use the solver.

Use [`set_parameters`](@ref) to set solver-specific options. Continuing the
example from above,
```julia
set_parameters(model, "Presolve" => 0, "Heuristics" => 0.01)
```
sets Gurobi's
[`Presolve`](https://www.gurobi.com/documentation/8.1/refman/presolve.html#parameter:Presolve)
parameter to zero and
[`Heuristics`](https://www.gurobi.com/documentation/8.1/refman/heuristics.html#parameter:Heuristics)
to 0.01.

The following solvers were compatible with JuMP up to release 0.18 but are
not yet compatible with the latest version because they do not implement the
new MathOptInterface API:

- [Alpine](https://github.com/lanl-ansi/Alpine.jl)
- [BARON](https://github.com/joehuchette/BARON.jl)
- [NLopt](https://github.com/JuliaOpt/NLopt.jl)
- [Pavito](https://github.com/JuliaOpt/Pavito.jl)
- [Pajarito](https://github.com/JuliaOpt/Pajarito.jl)

Solver-specific notes follow below.

### Artelys Knitro

Requires a license.

### BARON

Requires a license. A trial version is available for small problem instances.

### CDD

CDD can solve the problem both using `Float64` and `Rational{BigInt}`
arithmetics. The arithmetic used the type `T` given in `CDDLib.Optimizer{T}`.
Only `CDDLib.Optimizer{Float64}` can be used with JuMP as JuMP inputs the
problem in `Float64` arithmetics. Use [MOI](https://github.com/JuliaOpt/MathOptInterface.jl)
directly for `CDDLib.Optimizer{Rational{BigInt}}`.

### COIN-OR Cbc

Cbc supports "SOS" constraints.

### COSMO

COSMO can solve LPs, QPs, SOCPs and SDPs. It can handle SDPs with
 quadratic objective functions and supports chordal decomposition of large structured
PSD constraints. COSMO is a first order method that performs well on large problems
but has a low accuracy by default (``10^{−4}``).
See the [COSMO.jl](https://oxfordcontrol.github.io/COSMO.jl/stable/)
documentation for more information.

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

Requires a license (free for academic use).
The [Mosek interface](https://github.com/JuliaOpt/Mosek.jl) is maintained by
the Mosek team. (Thanks!)
Note that even if the package implementing MathOptInterface is `MosekTools`,
for consistency the MOI optimizer is called `Mosek.Optimizer` so do the
following to create a model with the Mosek solver:
```julia
julia> using MosekTools
julia> model = Model(Mosek.Optimizer)
```

### ProxSDP

ProxSDP solves general SDP problems by means of a first order proximal algorithm
based on the primal-dual hybrid gradient, also known as Chambolle-Pock method.
The main advantage of ProxSDP over other state-of-the-art solvers is the ability
to exploit the low-rank property inherent to several SDP problems. ProxSDP
is a first order solver and has low accuracy. See the [ProxSDP.jl](https://github.com/mariohsouto/ProxSDP.jl)
documentation for more information.

### SCS

SCS can be used by JuMP to solve LPs and SOCPs, and SDPs. SCS is a first order
solver and has low accuracy (``10^{−4}``) by default; see the SCS.jl
documentation for more information.

### SDPA

SDPA is a second order solver which comes in several variants. The main version has a C++ interface which [SDPA.jl](https://github.com/JuliaOpt/SDPA.jl) uses for efficiently communicating the problem instance to the solver. The three high-precision variants, SDPA-GMP (arbitrary precision), SDPA-QD ("quad-double" precision) and SDPA-DD ("double-double" precision) do not expose a library interface, but can used via [SDPAFamily.jl](https://github.com/ericphanson/SDPAFamily.jl) which writes and reads files to interact with the solver binary.

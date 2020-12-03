# Installation Guide

JuMP is a package for [Julia](https://julialang.org). To use JuMP, first
[download and install](https://julialang.org/downloads/) Julia.

!!! note
    This version of JuMP is compatible with Julia 1.0 and later.

From Julia, JuMP is installed using the built-in package manager:
```julia
import Pkg
Pkg.add("JuMP")
```
## Installing a solver

JuMP depends on solvers to solve optimization problems, and you will need to
install one before you can solve problems with JuMP. The table below lists the
currently available solvers.

Common choices of solvers are [Clp.jl](https://github.com/jump-dev/Clp.jl) for
linear programs, [Cbc.jl](https://github.com/jump-dev/Clp.jl) for mixed-integer
linear programs, and [Ipopt.jl](https://github.com/jump-dev/Clp.jl) for
non-linear programs.

Install a solver using the Julia package manager, replacing `"Clp"` by the
Julia package name as appropriate.
```julia
import Pkg
Pkg.add("Clp")
```

Once installed, you can use Clp as a solver with JuMP as follows, using
[`set_optimizer_attributes`](@ref) to set solver-specific options:
```julia
using JuMP
using Clp
model = Model(Clp.Optimizer)
set_optimizer_attributes(model, "LogLevel" => 1, "PrimalTolerance" => 1e-7)
```

!!! note
    Most packages follow the `ModuleName.Optimizer` naming convention, but
    exceptions may exist. See the README of the Julia package's Github
    repository for more details on how to use a particular solver, including any
    solver-specific options.

### Supported solvers

!!! note
    Most solvers are not written in Julia, and some require commercial licenses
    to use, so installation is often more complex. If a solver has a `Yes` in
    the `Manual Installation` column, the solver requires a manual installation
    step, such as downloading and installing a binary, or obtaining a commercial
    license. Consult the README of the relevant Julia package for more
    information. If the `Manual Installation` column is missing an entry,
    installing the Julia package will download and install any relevant solver
    binaries automatically, and you shouldn't need to do anything other than
    `Pkg.add`.

| Solver                                                                         | Julia Package                                                                    | Manual Installation | License  | Supports           |
| ------------------------------------------------------------------------------ | -------------------------------------------------------------------------------- | ------------------- | -------- | ------------------ |
| [Alpine](https://github.com/lanl-ansi/Alpine.jl)                               | [Alpine.jl](https://github.com/lanl-ansi/Alpine.jl)                              |     | Triad NS | MINLP                              |
| [Artelys Knitro](https://www.artelys.com/knitro)                               | [KNITRO.jl](https://github.com/jump-dev/KNITRO.jl)                               | Yes | Comm.    | LP, MILP, SOCP, MISOCP, NLP, MINLP |
| [BARON](http://minlp.com/baron)                                                | [BARON.jl](https://github.com/joehuchette/BARON.jl)                              | Yes | Comm.    | MINLP                              |
| [Cbc](https://github.com/coin-or/Cbc)                                          | [Cbc.jl](https://github.com/jump-dev/Cbc.jl)                                     |     | EPL      | MILP                               |
| [CDCS](https://github.com/oxfordcontrol/CDCS)                                  | [CDCS.jl](https://github.com/oxfordcontrol/CDCS.jl)                              |     | GPL      | LP, SOCP, SDP                      |
| [CDD](https://github.com/cddlib/cddlib)                                        | [CDDLib.jl](https://github.com/JuliaPolyhedra/CDDLib.jl)                         |     | GPL      | LP                                 |
| [Clp](https://github.com/coin-or/Clp)                                          | [Clp.jl](https://github.com/jump-dev/Clp.jl)                                     |     | EPL      | LP                                 |
| [COSMO](https://github.com/oxfordcontrol/COSMO.jl)                             | [COSMO.jl](https://github.com/oxfordcontrol/COSMO.jl)                            |     | Apache   | LP, QP, SOCP, SDP                  |
| [CPLEX](https://www.ibm.com/analytics/cplex-optimizer/)                        | [CPLEX.jl](https://github.com/jump-dev/CPLEX.jl)                                 | Yes | Comm.    | LP, MILP, SOCP, MISOCP             |
| [CSDP](https://github.com/coin-or/Csdp)                                        | [CSDP.jl](https://github.com/jump-dev/CSDP.jl)                                   |     | EPL      | LP, SDP                            |
| [EAGO](https://github.com/psorlab/EAGO.jl)                                     | [EAGO.jl](https://github.com/psorlab/EAGO.jl)                                    |     | CC BY-NC-SA | NLP                             |
| [ECOS](https://github.com/ifa-ethz/ecos)                                       | [ECOS.jl](https://github.com/jump-dev/ECOS.jl)                                   |     | GPL      | LP, SOCP                           |
| [FICO Xpress](https://www.fico.com/en/products/fico-xpress-optimization-suite) | [Xpress.jl](https://github.com/jump-dev/Xpress.jl)                               | Yes | Comm.    | LP, MILP, SOCP, MISOCP             |
| [GLPK](http://www.gnu.org/software/glpk/)                                      | [GLPK.jl](https://github.com/jump-dev/GLPK.jl)                                   |     | GPL      | LP, MILP                           |
| [Gurobi](https://gurobi.com)                                                   | [Gurobi.jl](https://github.com/jump-dev/Gurobi.jl)                               | Yes | Comm.    | LP, MILP, SOCP, MISOCP             |
| [Hypatia](https://github.com/chriscoey/Hypatia.jl)                             | [Hypatia.jl](https://github.com/chriscoey/Hypatia.jl)                            |     | MIT      | LP, SOCP, SDP                      |
| [Ipopt](https://github.com/coin-or/Ipopt)                                      | [Ipopt.jl](https://github.com/jump-dev/Ipopt.jl)                                 |     | EPL      | LP, QP, NLP                        |
| [Juniper](https://github.com/lanl-ansi/Juniper.jl)                             | [Juniper.jl](https://github.com/lanl-ansi/Juniper.jl)                            |     | MIT      | MISOCP, MINLP                      |
| [MOSEK](https://www.mosek.com/)                                                | [MosekTools.jl](https://github.com/jump-dev/MosekTools.jl)                       | Yes | Comm.    | LP, MILP, SOCP, MISOCP, SDP        |
| [NLopt](https://github.com/stevengj/nlopt)                                     | [NLopt.jl](https://github.com/JuliaOpt/NLopt.jl)                                 |     | GPL      | LP, QP, NLP                        |
| [OSQP](https://osqp.org/)                                                      | [OSQP.jl](https://github.com/oxfordcontrol/OSQP.jl)                              |     | Apache   | LP, QP                             |
| [Pavito](https://github.com/jump-dev/Pavito.jl)                                | [Pavito.jl](https://github.com/jump-dev/Pavito.jl)                               |     | MPL-2    | MICP                               |
| [ProxSDP](https://github.com/mariohsouto/ProxSDP.jl)                           | [ProxSDP.jl](https://github.com/mariohsouto/ProxSDP.jl)                          |     | MIT      | LP, SOCP, SDP                      |
| [SCIP](https://scip.zib.de/)                                                   | [SCIP.jl](https://github.com/SCIP-Interfaces/SCIP.jl)                            | Yes | ZIB      | MILP, MINLP                        |
| [SCS](https://github.com/cvxgrp/scs)                                           | [SCS.jl](https://github.com/jump-dev/SCS.jl)                                     |     | MIT      | LP, SOCP, SDP                      |
| [SDPA](http://sdpa.sourceforge.net/)                                           | [SDPA.jl](https://github.com/jump-dev/SDPA.jl), [SDPAFamily.jl](https://github.com/ericphanson/SDPAFamily.jl) |     | GPL | LP, SDP    |
| [SDPNAL](https://blog.nus.edu.sg/mattohkc/softwares/sdpnalplus/)               | [SDPNAL.jl](https://github.com/jump-dev/SDPNAL.jl)                               | Yes | CC BY-SA | LP, SDP                            |
| [SDPT3](https://blog.nus.edu.sg/mattohkc/softwares/sdpt3/)                     | [SDPT3.jl](https://github.com/jump-dev/SDPT3.jl)                                 | Yes | GPL      | LP, SOCP, SDP                      |
| [SeDuMi](http://sedumi.ie.lehigh.edu/)                                         | [SeDuMi.jl](https://github.com/jump-dev/SeDuMi.jl)                               | Yes | GPL      | LP, SOCP, SDP                      |
| [Tulip](https://github.com/ds4dm/Tulip.jl)                                     | [Tulip.jl](https://github.com/ds4dm/Tulip.jl)                                    |     | MPL-2    | LP                                 |

Where:

-   LP = Linear programming
-   QP = Quadratic programming
-   SOCP = Second-order conic programming (including problems with convex quadratic constraints and/or objective)
-   MILP = Mixed-integer linear programming
-   NLP = Nonlinear programming
-   MINLP = Mixed-integer nonlinear programming
-   SDP = Semidefinite programming
-   MISDP = Mixed-integer semidefinite programming

!!! note
    Developed a solver? This table is open for new contributions! Start by
    making a pull request to edit the [installation.md](https://github.com/jump-dev/JuMP.jl/blob/master/docs/src/installation.md)
    file.

!!! note
    Developing a solver? See [Interacting with solvers](@ref) and the
    [MathOptInterface docs](https://jump.dev/MathOptInterface.jl/stable/)
    for more details on how JuMP interacts with solvers. Please get in touch
    via the [Developer Chatroom](https://jump.dev/pages/governance/#developer-chatroom)
    with any questions about connecting new solvers with JuMP.

## AMPL and GAMS

Use [AmplNLWriter](https://github.com/jump-dev/AmplNLWriter.jl) to access
solvers that support the [nl format](https://en.wikipedia.org/wiki/Nl_(format)).
Such solvers include [Bonmin](https://github.com/coin-or/Bonmin) and
[Couenne](https://github.com/coin-or/Couenne). See a more complete list
[here](https://ampl.com/products/solvers/all-solvers-for-ampl/).

Use [GAMS.jl](https://github.com/GAMS-dev/gams.jl) to access solvers via an
installation of [GAMS](https://www.gams.com). Among them are:
[AlphaECP](https://www.gams.com/latest/docs/S_ALPHAECP.html),
[Antigone](https://www.gams.com/latest/docs/S_ANTIGONE.html),
[BARON](https://www.gams.com/latest/docs/S_BARON.html),
[CONOPT](https://www.gams.com/latest/docs/S_CONOPT.html),
[Couenne](https://www.gams.com/latest/docs/S_COUENNE.html),
[LocalSolver](https://www.gams.com/latest/docs/S_LOCALSOLVER.html),
[PATHNLP](https://www.gams.com/latest/docs/S_PATHNLP.html),
[SHOT](https://www.gams.com/latest/docs/S_SHOT.html),
[SNOPT](https://www.gams.com/latest/docs/S_SNOPT.html),
[SoPlex](https://www.gams.com/latest/docs/S_SOPLEX.html).
See a complete list [here](https://www.gams.com/latest/docs/S_MAIN.html).

!!! note
    [GAMS.jl](https://github.com/GAMS-dev/gams.jl) requires an installation of
    the commercial software [GAMS](https://www.gams.com) for which a
    [free community license](https://www.gams.com/latest/docs/UG_License.html#GAMS_Community_Licenses)
    exists.

## Solver-specific notes
### Artelys Knitro

Requires a license.

### BARON

Requires a license. A trial version is available for small problem instances.

### CDD

CDD can solve the problem both using `Float64` and `Rational{BigInt}`
arithmetics. The arithmetic used the type `T` given in `CDDLib.Optimizer{T}`.
Only `CDDLib.Optimizer{Float64}` can be used with JuMP as JuMP inputs the
problem in `Float64` arithmetics. Use [MOI](https://github.com/jump-dev/MathOptInterface.jl)
directly for `CDDLib.Optimizer{Rational{BigInt}}`.

### COIN-OR Cbc

Cbc supports "SOS" constraints.

### COSMO

COSMO can solve LPs, QPs, SOCPs and SDPs. It can handle SDPs with quadratic
objective functions and supports chordal decomposition of large structured PSD
constraints. COSMO is a first order method that performs well on large problems
but has a low accuracy by default (``10^{−4}``).
See the [COSMO.jl documentation](https://oxfordcontrol.github.io/COSMO.jl/stable/)
for more information.

### CPLEX

Requires a working installation of CPLEX with a license (free for faculty
members and graduate teaching assistants). CPLEX supports "SOS" constraints.

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
[FICO Academic Partner Program](https://fico.com/en/xpress-academic-license)).
Supports SOCP and "SOS" constraints.

### MOSEK

Requires a license (free for academic use). The [Mosek interface](https://github.com/MOSEK/Mosek.jl)
is maintained by the Mosek team. (Thanks!) Note that even if the package
implementing MathOptInterface is `MosekTools`, for consistency the MOI optimizer
is called `Mosek.Optimizer` so do the following to create a model with the Mosek
solver:
```julia
using MosekTools
model = Model(Mosek.Optimizer)
```

### ProxSDP

ProxSDP solves general SDP problems by means of a first order proximal algorithm
based on the primal-dual hybrid gradient, also known as Chambolle-Pock method.
The main advantage of ProxSDP over other state-of-the-art solvers is the ability
to exploit the low-rank property inherent to several SDP problems. ProxSDP is a
first order solver and has low accuracy. See the [ProxSDP.jl](https://github.com/mariohsouto/ProxSDP.jl)
documentation for more information.

### SCS

SCS can be used by JuMP to solve LPs and SOCPs, and SDPs. SCS is a first order
solver and has low accuracy (``10^{−4}``) by default; see the [SCS.jl](https://github.com/jump-dev/SCS.jl)
documentation for more information.

### SDPA

SDPA is a second order solver which comes in several variants. The main version
has a C++ interface which [SDPA.jl](https://github.com/jump-dev/SDPA.jl) uses
for efficiently communicating the problem instance to the solver. The three
high-precision variants, SDPA-GMP (arbitrary precision), SDPA-QD ("quad-double"
precision) and SDPA-DD ("double-double" precision) do not expose a library
interface, but can used via [SDPAFamily.jl](https://github.com/ericphanson/SDPAFamily.jl),
which writes and reads files to interact with the solver binary.

## Previously supported solvers

The following solvers were compatible with JuMP up to release 0.18 but are
not yet compatible with the latest version because they do not implement the
new MathOptInterface API:

- [Pajarito](https://github.com/JuliaOpt/Pajarito.jl)

Please join the [Developer Chatroom](https://jump.dev/pages/governance/#developer-chatroom)
if you have interest in reviving a previously supported solver.

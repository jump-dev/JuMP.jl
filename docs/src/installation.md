# Installation Guide

!!! info
    Installation troubles? Check the [Common installation issues](@ref) section
    below.

JuMP is a package for [Julia](https://julialang.org). To use JuMP, first
[download and install](https://julialang.org/downloads/) Julia.

!!! note
    This version of JuMP is compatible with Julia 1.0 and later.

From Julia, JuMP is installed using the built-in package manager:
```julia
import Pkg
Pkg.add("JuMP")
```

!!! tip
    We recommend you create a Pkg _environment_ for each project you use JuMP
    for, instead of adding lots of packages to the global environment. The
    [Pkg manager documentation](https://julialang.github.io/Pkg.jl/v1/environments/)
    has more information on this topic.

## Installing a solver

JuMP depends on solvers to solve optimization problems, and you will need to
install one before you can solve problems with JuMP. The table below lists the
currently available solvers.

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

Most solvers are not written in Julia, and some require commercial licenses to
use, so installation is often more complex.
  * If a solver has `Manual` in the `Installation` column, the solver requires a
    manual installation step, such as downloading and installing a binary, or
    obtaining a commercial license. Consult the README of the relevant Julia
    package for more information.
  * If the solver has `Manualᴹ` in the `Installation` column, the solver
    requires an installation of [MATLAB](https://www.mathworks.com/products/matlab.html).
  * If the `Installation` column is missing an entry, installing the Julia
    package will download and install any relevant solver binaries
    automatically, and you shouldn't need to do anything other than `Pkg.add`.

Solvers with a missing entry in the `Julia Package` column are written in Julia.
The link in the `Solver` column is the corresponding Julia package.

| Solver                                                                         | Julia Package                                                                    | Installation | License | Supports             |
| ------------------------------------------------------------------------------ | -------------------------------------------------------------------------------- | ------------ | ------- | ---------------------|
| [Alpine.jl](https://github.com/lanl-ansi/Alpine.jl)                            |                                                                                  |        | Triad NS | (MI)NLP                   |
| [Artelys Knitro](https://www.artelys.com/knitro)                               | [KNITRO.jl](https://github.com/jump-dev/KNITRO.jl)                               | Manual | Comm.    | (MI)LP, (MI)SOCP, (MI)NLP |
| [BARON](http://minlp.com/baron)                                                | [BARON.jl](https://github.com/joehuchette/BARON.jl)                              | Manual | Comm.    | (MI)NLP                   |
| [Bonmin](http://github.com/coin-or/Bonmin)                                     | [AmplNLWriter.jl](https://github.com/jump-dev/AmplNLWriter.jl)                   |        | EPL      | (MI)NLP                   |
| [Cbc](https://github.com/coin-or/Cbc)                                          | [Cbc.jl](https://github.com/jump-dev/Cbc.jl)                                     |        | EPL      | (MI)LP                    |
| [CDCS](https://github.com/oxfordcontrol/CDCS)                                  | [CDCS.jl](https://github.com/oxfordcontrol/CDCS.jl)                              | Manualᴹ | GPL     | LP, SOCP, SDP             |
| [CDD](https://github.com/cddlib/cddlib)                                        | [CDDLib.jl](https://github.com/JuliaPolyhedra/CDDLib.jl)                         |        | GPL      | LP                        |
| [Clp](https://github.com/coin-or/Clp)                                          | [Clp.jl](https://github.com/jump-dev/Clp.jl)                                     |        | EPL      | LP                        |
| [COSMO.jl](https://github.com/oxfordcontrol/COSMO.jl)                          |                                                                                  |        | Apache   | LP, QP, SOCP, SDP         |
| [Couenne](http://github.com/coin-or/Couenne)                                   | [AmplNLWriter.jl](https://github.com/jump-dev/AmplNLWriter.jl)                   |        | EPL      | (MI)NLP                   |
| [CPLEX](https://www.ibm.com/analytics/cplex-optimizer/)                        | [CPLEX.jl](https://github.com/jump-dev/CPLEX.jl)                                 | Manual | Comm.    | (MI)LP, (MI)SOCP          |
| [CSDP](https://github.com/coin-or/Csdp)                                        | [CSDP.jl](https://github.com/jump-dev/CSDP.jl)                                   |        | EPL      | LP, SDP                   |
| [EAGO.jl](https://github.com/psorlab/EAGO.jl)                                  |                                                                                  |        | MIT | NLP                    |
| [ECOS](https://github.com/ifa-ethz/ecos)                                       | [ECOS.jl](https://github.com/jump-dev/ECOS.jl)                                   |        | GPL      | LP, SOCP                  |
| [FICO Xpress](https://www.fico.com/en/products/fico-xpress-optimization-suite) | [Xpress.jl](https://github.com/jump-dev/Xpress.jl)                               | Manual | Comm.    | (MI)LP, (MI)SOCP          |
| [GLPK](http://www.gnu.org/software/glpk/)                                      | [GLPK.jl](https://github.com/jump-dev/GLPK.jl)                                   |        | GPL      | (MI)LP                    |
| [Gurobi](https://gurobi.com)                                                   | [Gurobi.jl](https://github.com/jump-dev/Gurobi.jl)                               | Manual | Comm.    | (MI)LP, (MI)SOCP          |
| [HiGHS](https://github.com/ERGO-Code/HiGHS)                                    | [HiGHS.jl](https://github.com/jump-dev/HiGHS.jl)                                 |        |MIT       | LP                        |
| [Hypatia.jl](https://github.com/chriscoey/Hypatia.jl)                          |                                                                                  |        | MIT      | LP, SOCP, SDP             |
| [Ipopt](https://github.com/coin-or/Ipopt)                                      | [Ipopt.jl](https://github.com/jump-dev/Ipopt.jl)                                 |        | EPL      | LP, QP, NLP               |
| [Juniper.jl](https://github.com/lanl-ansi/Juniper.jl)                          |                                                                                  |        | MIT      | (MI)SOCP, (MI)NLP         |
| [MOSEK](https://www.mosek.com/)                                                | [MosekTools.jl](https://github.com/jump-dev/MosekTools.jl)                       | Manual | Comm.    | (MI)LP, (MI)SOCP, SDP     |
| [NLopt](https://github.com/stevengj/nlopt)                                     | [NLopt.jl](https://github.com/JuliaOpt/NLopt.jl)                                 |        | GPL      | LP, QP, NLP               |
| [OSQP](https://osqp.org/)                                                      | [OSQP.jl](https://github.com/oxfordcontrol/OSQP.jl)                              |        | Apache   | LP, QP                    |
| [PATH](http://pages.cs.wisc.edu/~ferris/path.html)                             | [PATHSolver.jl](https://github.com/chkwon/PATHSolver.jl)                 |        | MIT      | MCP                       |
| [Pavito.jl](https://github.com/jump-dev/Pavito.jl)                             |                                                                                  |        | MPL-2    | (MI)NLP                   |
| [ProxSDP.jl](https://github.com/mariohsouto/ProxSDP.jl)                        |                                                                                  |        | MIT      | LP, SOCP, SDP             |
| [SCIP](https://scipopt.org/)                                                   | [SCIP.jl](https://github.com/scipopt/SCIP.jl)                            |        | ZIB      | (MI)LP, (MI)NLP           |
| [SCS](https://github.com/cvxgrp/scs)                                           | [SCS.jl](https://github.com/jump-dev/SCS.jl)                                     |        | MIT      | LP, SOCP, SDP             |
| [SDPA](http://sdpa.sourceforge.net/)                                           | [SDPA.jl](https://github.com/jump-dev/SDPA.jl), [SDPAFamily.jl](https://github.com/ericphanson/SDPAFamily.jl) |  | GPL | LP, SDP |
| [SDPNAL](https://blog.nus.edu.sg/mattohkc/softwares/sdpnalplus/)               | [SDPNAL.jl](https://github.com/jump-dev/SDPNAL.jl)                               | Manualᴹ | CC BY-SA | LP, SDP                  |
| [SDPT3](https://blog.nus.edu.sg/mattohkc/softwares/sdpt3/)                     | [SDPT3.jl](https://github.com/jump-dev/SDPT3.jl)                                 | Manualᴹ | GPL      | LP, SOCP, SDP            |
| [SeDuMi](http://sedumi.ie.lehigh.edu/)                                         | [SeDuMi.jl](https://github.com/jump-dev/SeDuMi.jl)                               | Manualᴹ | GPL      | LP, SOCP, SDP            |
| [Tulip.jl](https://github.com/ds4dm/Tulip.jl)                                  |                                                                                  |        | MPL-2     | LP                       |

Where:
- LP = Linear programming
- QP = Quadratic programming
- SOCP = Second-order conic programming (including problems with convex
  quadratic constraints and/or objective)
- MCP = Mixed-complementarity programming
- NLP = Nonlinear programming
- SDP = Semidefinite programming
- (MI)XXX = Mixed-integer equivalent of problem type `XXX`

!!! note
    Developed a solver or solver wrapper? This table is open for new
    contributions! Start by making a pull request to edit the [installation.md](https://github.com/jump-dev/JuMP.jl/blob/master/docs/src/installation.md)
    file.

!!! note
    Developing a solver or solver wrapper? See [Models](@ref jump_models) and the
    [MathOptInterface docs](https://jump.dev/MathOptInterface.jl/stable/) for
    more details on how JuMP interacts with solvers. Please get in touch via the
    [Developer Chatroom](https://jump.dev/pages/governance/#developer-chatroom)
    with any questions about connecting new solvers with JuMP.

### Solver-specific notes

* Artelys Knitro

Requires a license.

* BARON

Requires a license. A trial version is available for small problem instances.

* CDD

CDD can solve the problem both using `Float64` and `Rational{BigInt}`
arithmetics. The arithmetic used the type `T` given in `CDDLib.Optimizer{T}`.
Only `CDDLib.Optimizer{Float64}` can be used with JuMP as JuMP inputs the
problem in `Float64` arithmetics. Use [MOI](https://github.com/jump-dev/MathOptInterface.jl)
directly for `CDDLib.Optimizer{Rational{BigInt}}`.

* COIN-OR Cbc

Cbc supports "SOS" constraints.

* COSMO

COSMO can solve LPs, QPs, SOCPs and SDPs. It can handle SDPs with quadratic
objective functions and supports chordal decomposition of large structured PSD
constraints. COSMO is a first order method that performs well on large problems
but has a low accuracy by default (``10^{−4}``).
See the [COSMO.jl documentation](https://oxfordcontrol.github.io/COSMO.jl/stable/)
for more information.

* CPLEX

Requires a working installation of CPLEX with a license (free for faculty
members and graduate teaching assistants). CPLEX supports "SOS" constraints.

* ECOS

ECOS can be used by JuMP to solve LPs and SOCPs. ECOS does not support general
quadratic objectives or constraints, only second-order conic constraints
specified by using the `SecondOrderCone` set.

* Gurobi

Requires a working installation of Gurobi with an activated license (free for
academic use). Gurobi supports "SOS" constraints.

* FICO Xpress

Requires a working installation of Xpress with an active license (it is possible
to get a license for academic use, see
[FICO Academic Partner Program](https://fico.com/en/xpress-academic-license)).
Supports SOCP and "SOS" constraints.

* MOSEK

Requires a license (free for academic use). The [Mosek interface](https://github.com/MOSEK/Mosek.jl)
is maintained by the Mosek team. (Thanks!) Note that even if the package
implementing MathOptInterface is `MosekTools`, for consistency the MOI optimizer
is called `Mosek.Optimizer` so do the following to create a model with the Mosek
solver:
```julia
using MosekTools
model = Model(Mosek.Optimizer)
```

* ProxSDP

ProxSDP solves general SDP problems by means of a first order proximal algorithm
based on the primal-dual hybrid gradient, also known as Chambolle-Pock method.
The main advantage of ProxSDP over other state-of-the-art solvers is the ability
to exploit the low-rank property inherent to several SDP problems. ProxSDP is a
first order solver and has low accuracy. See the [ProxSDP.jl](https://github.com/mariohsouto/ProxSDP.jl)
documentation for more information.

* SCS

SCS can be used by JuMP to solve LPs and SOCPs, and SDPs. SCS is a first order
solver and has low accuracy (``10^{−4}``) by default; see the [SCS.jl](https://github.com/jump-dev/SCS.jl)
documentation for more information.

* SDPA

SDPA is a second order solver which comes in several variants. The main version
has a C++ interface which [SDPA.jl](https://github.com/jump-dev/SDPA.jl) uses
for efficiently communicating the problem instance to the solver. The three
high-precision variants, SDPA-GMP (arbitrary precision), SDPA-QD ("quad-double"
precision) and SDPA-DD ("double-double" precision) do not expose a library
interface, but can used via [SDPAFamily.jl](https://github.com/ericphanson/SDPAFamily.jl),
which writes and reads files to interact with the solver binary.

## AMPL-based solvers

Use [AmplNLWriter](https://github.com/jump-dev/AmplNLWriter.jl) to access
solvers that support the [nl format](https://en.wikipedia.org/wiki/Nl_(format)).

Some solvers, such as [Bonmin](https://github.com/coin-or/Bonmin) and
[Couenne](https://github.com/coin-or/Couenne) can be installed via the Julia
package manager. Others need to be manually installed.

Consult the AMPL documentation for a [complete list of supported solvers](https://ampl.com/products/solvers/all-solvers-for-ampl/).

## GAMS-based solvers

Use [GAMS.jl](https://github.com/GAMS-dev/gams.jl) to access solvers available
through [GAMS](https://www.gams.com). Such solvers include:
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

## Previously supported solvers

The following solvers were compatible with JuMP up to release 0.18 but are
not yet compatible with the latest version because they do not implement the
new MathOptInterface API:

- [Pajarito](https://github.com/JuliaOpt/Pajarito.jl)

Please join the [Developer Chatroom](https://jump.dev/pages/governance/#developer-chatroom)
if you have interest in reviving a previously supported solver.

## Common installation issues

!!! tip
    When in doubt, run `import Pkg; Pkg.update()` to see if updating your
    packages fixes the issue. Remember you will need to exit Julia and start a
    new session for the changes to take effect.


### Check the version of your packages

Each package is versioned with a [three-part number](https://semver.org) of the
form `vX.Y.Z`. You can check which versions you have installed with
`import Pkg; Pkg.status()`.

This should almost always be the most-recent release. You can check the releases
of a package by going to the relevant Github page, and navigating to the
"releases" page. For example, the list of JuMP releases is available at:
[https://github.com/jump-dev/JuMP.jl/releases](https://github.com/jump-dev/JuMP.jl/releases).

If you post on the [community forum](https://discourse.julialang.org/c/domain/opt/13),
please include the output of `Pkg.status()`!

### Unsatisfiable requirements detected

Did you get an error like `Unsatisfiable requirements detected for package JuMP`?
The Pkg documentation has a [section on how to understand and manage these conflicts](https://julialang.github.io/Pkg.jl/v1/managing-packages/#conflicts).

### Installing new packages can make JuMP downgrade to an earlier version

Another common complaint is that after adding a new package, code that
previously worked no longer works.

This usually happens because the new package is not compatible with the latest
version of JuMP. Therefore, the package manager rolls-back JuMP to an earlier
version! Here's an example.

First, we add JuMP:
```julia
(jump_example) pkg> add JuMP
  Resolving package versions...
Updating `~/jump_example/Project.toml`
  [4076af6c] + JuMP v0.21.5
Updating `~/jump_example/Manifest.toml`
  ... lines omitted ...
```
The `+ JuMP v0.21.5` line indicates that JuMP has been added at version
`0.21.5`. However, watch what happens when we add [JuMPeR](https://github.com/iainnz/JuMPeR.jl):
```julia
(jump_example) pkg> add JuMPeR
  Resolving package versions...
Updating `~/jump_example/Project.toml`
  [4076af6c] ↓ JuMP v0.21.5 ⇒ v0.18.6
  [707a9f91] + JuMPeR v0.6.0
Updating `~/jump_example/Manifest.toml`
  ... lines omitted ...
```
JuMPeR gets added at version `0.6.0` (`+ JuMPeR v0.6.0`), but JuMP gets
downgraded from `0.21.5` to `0.18.6` (`↓ JuMP v0.21.5 ⇒ v0.18.6`)! The reason
for this is that JuMPeR doesn't support a version of JuMP newer than `0.18.6`.

!!! tip
    Pay careful attention to the output of the package manager when adding new
    packages, especially when you see a package being downgraded!

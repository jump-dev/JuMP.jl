# Installation Guide

This guide explains how to install Julia and  JuMP. If you have installation
troubles, read the [Common installation issues](@ref) section below.

## Install Julia

JuMP is a package for [Julia](https://julialang.org). To use JuMP, first
[download and install](https://julialang.org/downloads/) Julia.

!!! tip
    If you are new to Julia, read our [Getting started with Julia](@ref)
    tutorial.

### Choosing a version

You can install the "Current stable release" or the "Long-term support (LTS)
release."

 * The "Current stable release" is the latest release of Julia. It has access to
   newer features, and is likely faster.
 * The "Long-term support release" is an older version of Julia that has
   continued to receive bug and security fixes. However, it may not have the
   latest features or performance improvements.

For most users, you should install the "Current stable release," and whenever
Julia releases a new version of the current stable release, you should update
your version of Julia. Note that any code you write on one version of the
current stable release will continue to work on all subsequent releases.

For users in restricted software environments (for example, your enterprise IT
controls what software you can install), you may be better off installing the
long-term support release because you will not have to update Julia as
frequently.

## Install JuMP

JuMP is installed using the built-in Julia package manager. Launch Julia, and
then enter the following at the `julia>` prompt:
```julia
julia> import Pkg

julia> Pkg.add("JuMP")
```

!!! tip
    We recommend you create a Pkg _environment_ for each project you use JuMP
    for, instead of adding lots of packages to the global environment. The
    [Pkg manager documentation](https://julialang.github.io/Pkg.jl/v1/environments/)
    has more information on this topic.

When we release a new version of JuMP, you can update with:
```julia
julia> import Pkg

julia> Pkg.update("JuMP")
```

## Install a solver

JuMP depends on solvers to solve optimization problems. Therefore, you will need
to install one before you can solve problems with JuMP.

Install a solver using the Julia package manager, replacing `"HiGHS"` by the
Julia package name as appropriate.
```julia
julia> import Pkg

julia> Pkg.add("HiGHS")
```

Once installed, you can use HiGHS as a solver with JuMP as follows, using
[`set_attribute`](@ref) to set solver-specific options:
```julia
julia> using JuMP

julia> using HiGHS

julia> model = Model(HiGHS.Optimizer);

julia> set_attribute(model, "output_flag" => false)

julia> set_attribute(model, "primal_feasibility_tolerance" => 1e-8)
```

!!! note
    Most packages follow the `ModuleName.Optimizer` naming convention, but
    exceptions may exist. See the README of the Julia package's GitHub
    repository for more details on how to use a particular solver, including any
    solver-specific options.

## Supported solvers

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
| [Clarabel.jl](https://github.com/oxfordcontrol/Clarabel.jl)                    |                                                                                  |        | Apache   | LP, QP, SOCP, SDP         |
| [Clp](https://github.com/coin-or/Clp)                                          | [Clp.jl](https://github.com/jump-dev/Clp.jl)                                     |        | EPL      | LP                        |
| [COPT](https://www.shanshu.ai/copt)                                            | [COPT.jl](https://github.com/COPT-Public/COPT.jl)                                |        | Comm.    | (MI)LP, SOCP, SDP         |
| [COSMO.jl](https://github.com/oxfordcontrol/COSMO.jl)                          |                                                                                  |        | Apache   | LP, QP, SOCP, SDP         |
| [Couenne](http://github.com/coin-or/Couenne)                                   | [AmplNLWriter.jl](https://github.com/jump-dev/AmplNLWriter.jl)                   |        | EPL      | (MI)NLP                   |
| [CPLEX](https://www.ibm.com/analytics/cplex-optimizer/)                        | [CPLEX.jl](https://github.com/jump-dev/CPLEX.jl)                                 | Manual | Comm.    | (MI)LP, (MI)SOCP          |
| [CSDP](https://github.com/coin-or/Csdp)                                        | [CSDP.jl](https://github.com/jump-dev/CSDP.jl)                                   |        | EPL      | LP, SDP                   |
| [DAQP](https://github.com/darnstrom/daqp)                                      | [DAQP.jl](https://github.com/darnstrom/DAQP.jl)                                  |        | MIT      | (Mixed-binary) QP         |
| [DSDP](http://www.mcs.anl.gov/hs/software/DSDP/)                               | [DSDP.jl](https://github.com/jump-dev/DSDP.jl)                                   |        | DSDP     | LP, SDP                   |
| [EAGO.jl](https://github.com/psorlab/EAGO.jl)                                  |                                                                                  |        | MIT      | (MI)NLP                   |
| [ECOS](https://github.com/ifa-ethz/ecos)                                       | [ECOS.jl](https://github.com/jump-dev/ECOS.jl)                                   |        | GPL      | LP, SOCP                  |
| [FICO Xpress](https://www.fico.com/en/products/fico-xpress-optimization-suite) | [Xpress.jl](https://github.com/jump-dev/Xpress.jl)                               | Manual | Comm.    | (MI)LP, (MI)SOCP          |
| [GLPK](http://www.gnu.org/software/glpk/)                                      | [GLPK.jl](https://github.com/jump-dev/GLPK.jl)                                   |        | GPL      | (MI)LP                    |
| [Gurobi](https://gurobi.com)                                                   | [Gurobi.jl](https://github.com/jump-dev/Gurobi.jl)                               | Manual | Comm.    | (MI)LP, (MI)SOCP          |
| [HiGHS](https://github.com/ERGO-Code/HiGHS)                                    | [HiGHS.jl](https://github.com/jump-dev/HiGHS.jl)                                 |        | MIT      | (MI)LP, QP                |
| [Hypatia.jl](https://github.com/chriscoey/Hypatia.jl)                          |                                                                                  |        | MIT      | LP, SOCP, SDP             |
| [Ipopt](https://github.com/coin-or/Ipopt)                                      | [Ipopt.jl](https://github.com/jump-dev/Ipopt.jl)                                 |        | EPL      | LP, QP, NLP               |
| [Juniper.jl](https://github.com/lanl-ansi/Juniper.jl)                          |                                                                                  |        | MIT      | (MI)SOCP, (MI)NLP         |
| [Loraine.jl](https://github.com/kocvara/Loraine.jl)                            |                                                                                  |        | MIT      | LP, SDP                   |
| [MadNLP.jl](https://github.com/sshin23/MadNLP.jl)                              |                                                                                  |        | MIT      | LP, QP, NLP               |
| [MAiNGO](https://git.rwth-aachen.de/avt-svt/public/maingo)                     | [MAiNGO.jl](https://github.com/MAiNGO-github/MAiNGO.jl)                          |        | EPL 2.0  |(MI)NLP                    |
| [Manopt.jl](https://github.com/JuliaManifolds/Manopt.jl)                       |                                                                                  |        | MIT      | NLP                       |
| [MiniZinc](https://www.minizinc.org/)                                          | [MiniZinc.jl](https://github.com/jump-dev/MiniZinc.jl)                           | Manual | MPL-2    | CP-SAT                    |
| [Minotaur](https://github.com/coin-or/minotaur)                                | [AmplNLWriter.jl](https://github.com/jump-dev/AmplNLWriter.jl)                   | Manual | BSD-like | (MI)NLP                   |
| [MOSEK](https://www.mosek.com/)                                                | [MosekTools.jl](https://github.com/jump-dev/MosekTools.jl)                       | Manual | Comm.    | (MI)LP, (MI)SOCP, SDP     |
| [NLopt](https://github.com/stevengj/nlopt)                                     | [NLopt.jl](https://github.com/JuliaOpt/NLopt.jl)                                 |        | GPL      | LP, QP, NLP               |
| [Octeract](https://octeract.gg)                                                | [AmplNLWriter.jl](https://github.com/jump-dev/AmplNLWriter.jl)                   |        | Comm.    | (MI)NLP                   |
| [Optim.jl](https://github.com/JuliaNLSolvers/Optim.jl/)                        |                                                                                  |        | MIT      | NLP                       |
| [OSQP](https://osqp.org/)                                                      | [OSQP.jl](https://github.com/oxfordcontrol/OSQP.jl)                              |        | Apache   | LP, QP                    |
| [PATH](http://pages.cs.wisc.edu/~ferris/path.html)                             | [PATHSolver.jl](https://github.com/chkwon/PATHSolver.jl)                         |        | MIT      | MCP                       |
| [Pajarito.jl](https://github.com/jump-dev/Pajarito.jl)                         |                                                                                  |        | MPL-2    | (MI)NLP, (MI)SOCP, (MI)SDP |
| [Pavito.jl](https://github.com/jump-dev/Pavito.jl)                             |                                                                                  |        | MPL-2    | (MI)NLP                   |
| [Penbmi](http://www.penopt.com/penbmi.html)                                    | [Penopt.jl](https://github.com/jump-dev/Penopt.jl/)                              |        | Comm.    | Bilinear SDP              |
| [Percival.jl](https://github.com/JuliaSmoothOptimizers/Percival.jl/)           |                                                                                  |        | MPL-2    | NLP                       |
| [PolyJuMP.KKT](https://github.com/jump-dev/PolyJuMP.jl)                        | [PolyJuMP.jl](https://github.com/jump-dev/PolyJuMP.jl)                           |        | MIT      | NLP                       |
| [PolyJuMP.QCQP](https://github.com/jump-dev/PolyJuMP.jl)                       | [PolyJuMP.jl](https://github.com/jump-dev/PolyJuMP.jl)                           |        | MIT      | NLP                       |
| [ProxSDP.jl](https://github.com/mariohsouto/ProxSDP.jl)                        |                                                                                  |        | MIT      | LP, SOCP, SDP             |
| [RAPOSa](https://raposa.usc.es/)                                               | [AmplNLWriter.jl](https://github.com/jump-dev/AmplNLWriter.jl)                   | Manual | RAPOSa   | (MI)NLP                   |
| [SCIP](https://scipopt.org/)                                                   | [SCIP.jl](https://github.com/scipopt/SCIP.jl)                                    |        | Apache   | (MI)LP, (MI)NLP           |
| [SCS](https://github.com/cvxgrp/scs)                                           | [SCS.jl](https://github.com/jump-dev/SCS.jl)                                     |        | MIT      | LP, QP, SOCP, SDP         |
| [SDPA](http://sdpa.sourceforge.net/)                                           | [SDPA.jl](https://github.com/jump-dev/SDPA.jl), [SDPAFamily.jl](https://github.com/ericphanson/SDPAFamily.jl) |  | GPL | LP, SDP |
| [SDPLR](https://github.com/sburer/sdplr)                                       | [SDPLR.jl](https://github.com/jump-dev/SDPLR.jl)                                 |        | GPL      | LP, SDP                   |
| [SDPNAL](https://blog.nus.edu.sg/mattohkc/softwares/sdpnalplus/)               | [SDPNAL.jl](https://github.com/jump-dev/SDPNAL.jl)                               | Manualᴹ | CC BY-SA | LP, SDP                  |
| [SDPT3](https://blog.nus.edu.sg/mattohkc/softwares/sdpt3/)                     | [SDPT3.jl](https://github.com/jump-dev/SDPT3.jl)                                 | Manualᴹ | GPL      | LP, SOCP, SDP            |
| [SeDuMi](http://sedumi.ie.lehigh.edu/)                                         | [SeDuMi.jl](https://github.com/jump-dev/SeDuMi.jl)                               | Manualᴹ | GPL      | LP, SOCP, SDP            |
| [StatusSwitchingQP.jl](https://github.com/PharosAbad/StatusSwitchingQP.jl)     |                                                                                  |        | MIT       | LP, QP                   |
| [Tulip.jl](https://github.com/ds4dm/Tulip.jl)                                  |                                                                                  |        | MPL-2     | LP                       |

Where:
- LP = Linear programming
- QP = Quadratic programming
- SOCP = Second-order conic programming (including problems with convex
  quadratic constraints or objective)
- MCP = Mixed-complementarity programming
- NLP = Nonlinear programming
- SDP = Semidefinite programming
- (MI)XXX = Mixed-integer equivalent of problem type `XXX`
- CP-SAT = Constraint programming and Boolean satisfiability

!!! note
    Developed a solver or solver wrapper? This table is open for new
    contributions. Edit the [installation.md](https://github.com/jump-dev/JuMP.jl/blob/master/docs/src/installation.md)
    file, and use the checklist [Adding a new solver to the documentation](@ref)
    when opening the pull request.

!!! note
    Developing a solver or solver wrapper? See [Models](@ref jump_models) and the
    [MathOptInterface docs](https://jump.dev/MathOptInterface.jl/stable/) for
    more details on how JuMP interacts with solvers. Please get in touch via the
    [Developer Chatroom](https://jump.dev/pages/governance/#developer-chatroom)
    with any questions about connecting new solvers with JuMP.

## AMPL-based solvers

Use [AmplNLWriter](https://github.com/jump-dev/AmplNLWriter.jl) to access
solvers that support the [NL format](https://en.wikipedia.org/wiki/Nl_(format)).

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

## NEOS-based solvers

Use [NEOSServer.jl](https://github.com/odow/NEOSServer.jl) to access solvers
available through the [NEOS Server](https://neos-server.org).

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
of a package by going to the relevant GitHub page, and navigating to the
"releases" page. For example, the list of JuMP releases is available at:
[https://github.com/jump-dev/JuMP.jl/releases](https://github.com/jump-dev/JuMP.jl/releases).

If you post on the [community forum](https://jump.dev/forum), please include the
output of `Pkg.status()`.

### Unsatisfiable requirements detected

Did you get an error like `Unsatisfiable requirements detected for package JuMP`?
The Pkg documentation has a [section on how to understand and manage these conflicts](https://julialang.github.io/Pkg.jl/v1/managing-packages/#conflicts).

### Installing new packages can make JuMP downgrade to an earlier version

Another common complaint is that after adding a new package, code that
previously worked no longer works.

This usually happens because the new package is not compatible with the latest
version of JuMP. Therefore, the package manager rolls-back JuMP to an earlier
version. Here's an example.

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
    packages, especially when you see a package being downgraded.

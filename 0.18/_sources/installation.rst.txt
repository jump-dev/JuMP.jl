
.. include:: warn.rst

.. _jump-installation:

------------------
Installation Guide
------------------

This guide will briefly guide you through installing Julia, JuMP and[a] solver[s] of your choice.

Getting Julia
^^^^^^^^^^^^^

At the time of writing this documentation the latest release of Julia is version ``0.6``, which is the version required by JuMP. You can easily build from source on OS X and Linux, but the binaries will work well for most people.

Download links and more detailed instructions are available on the `Julia website <http://julialang.org>`_.

Getting JuMP
^^^^^^^^^^^^

Once you've installed Julia, installing JuMP is simple. Julia has a git-based package system. To use it, open Julia in interactive mode (i.e. ``julia`` at the command line) and use the package manager::

    julia> Pkg.add("JuMP")

This command checks `METADATA.jl <https://github.com/JuliaLang/METADATA.jl>`_ to determine what the most recent version of JuMP is and then downloads it from its repository on GitHub.

To start using JuMP (after installing a solver), it should be imported into the local scope::

    julia> using JuMP

Getting Solvers
^^^^^^^^^^^^^^^

Solver support in Julia is currently provided by writing a solver-specific package that provides a very thin wrapper around the solver's C interface and providing a standard interface that JuMP can call. If you are interested in providing an interface to your solver, please get in touch. The table below lists the currently supported solvers and their capabilities.



.. _jump-solvertable:

+----------------------------------------------------------------------------------+---------------------------------------------------------------------------------+-----------------------------+-------------+----+------+------+-----+-------+-----+
| Solver                                                                           | Julia Package                                                                   | ``solver=``                 | License     | LP | SOCP | MILP | NLP | MINLP | SDP |
+==================================================================================+=================================================================================+=============================+=============+====+======+======+=====+=======+=====+
| `Artelys Knitro <http://artelys.com/en/optimization-tools/knitro>`_              | `KNITRO.jl <https://github.com/JuliaOpt/KNITRO.jl>`_                            | ``KnitroSolver()``          |  Comm.      |    |      |      |  X  |   X   |     |
+----------------------------------------------------------------------------------+---------------------------------------------------------------------------------+-----------------------------+-------------+----+------+------+-----+-------+-----+
| `BARON <http://archimedes.cheme.cmu.edu/?q=baron>`_                              | `BARON.jl <https://github.com/joehuchette/BARON.jl>`_                           |  ``BaronSolver()``          |  Comm.      |    |      |      |  X  |   X   |     |
+----------------------------------------------------------------------------------+---------------------------------------------------------------------------------+-----------------------------+-------------+----+------+------+-----+-------+-----+
| `Bonmin <https://projects.coin-or.org/Bonmin>`_                                  | `AmplNLWriter.jl <https://github.com/JuliaOpt/AmplNLWriter.jl>`_                | ``AmplNLWriter(...)`` *     |  EPL        | X  |      |  X   |  X  |   X   |     |
+                                                                                  +---------------------------------------------------------------------------------+-----------------------------+             +    +      +      +     +       +     +
|                                                                                  | `CoinOptServices.jl <https://github.com/JuliaOpt/CoinOptServices.jl>`_          | ``OsilBonminSolver()``      |             |    |      |      |     |       |     |
+----------------------------------------------------------------------------------+---------------------------------------------------------------------------------+-----------------------------+-------------+----+------+------+-----+-------+-----+
| `Cbc <https://projects.coin-or.org/Cbc>`_                                        | `Cbc.jl <https://github.com/JuliaOpt/Cbc.jl>`_                                  | ``CbcSolver()``             |  EPL        |    |      |  X   |     |       |     |
+----------------------------------------------------------------------------------+---------------------------------------------------------------------------------+-----------------------------+-------------+----+------+------+-----+-------+-----+
| `Clp <https://projects.coin-or.org/Clp>`_                                        | `Clp.jl <https://github.com/JuliaOpt/Clp.jl>`_                                  | ``ClpSolver()``             |  EPL        | X  |      |      |     |       |     |
+----------------------------------------------------------------------------------+---------------------------------------------------------------------------------+-----------------------------+-------------+----+------+------+-----+-------+-----+
|  `Couenne <https://projects.coin-or.org/Couenne>`_                               | `AmplNLWriter.jl <https://github.com/JuliaOpt/AmplNLWriter.jl>`_                | ``AmplNLWriter(...)`` **    |  EPL        | X  |      |  X   |  X  |   X   |     |
+                                                                                  +---------------------------------------------------------------------------------+-----------------------------+             +    +      +      +     +       +     +
|                                                                                  | `CoinOptServices.jl <https://github.com/JuliaOpt/CoinOptServices.jl>`_          | ``OsilCouenneSolver()``     |             |    |      |      |     |       |     |
+----------------------------------------------------------------------------------+---------------------------------------------------------------------------------+-----------------------------+-------------+----+------+------+-----+-------+-----+
| `CPLEX <http://www-01.ibm.com/software/commerce/optimization/cplex-optimizer/>`_ | `CPLEX.jl <https://github.com/JuliaOpt/CPLEX.jl>`_                              | ``CplexSolver()``           |  Comm.      | X  |  X   |  X   |     |       |     |
+----------------------------------------------------------------------------------+---------------------------------------------------------------------------------+-----------------------------+-------------+----+------+------+-----+-------+-----+
| `ECOS <https://github.com/ifa-ethz/ecos>`_                                       | `ECOS.jl <https://github.com/JuliaOpt/ECOS.jl>`_                                |  ``ECOSSolver()``           |  GPL        | X  |  X   |      |     |       |     |
+----------------------------------------------------------------------------------+---------------------------------------------------------------------------------+-----------------------------+-------------+----+------+------+-----+-------+-----+
| `FICO Xpress <http://www.fico.com/en/products/fico-xpress-optimization-suite>`_  | `Xpress.jl <https://github.com/JuliaOpt/Xpress.jl>`_                            |  ``Xpress.XpressSolver()``  |  Comm.      | X  |   X  |  X   |     |       |     |
+----------------------------------------------------------------------------------+---------------------------------------------------------------------------------+-----------------------------+-------------+----+------+------+-----+-------+-----+
| `GLPK <http://www.gnu.org/software/glpk/>`_                                      | `GLPKMath... <https://github.com/JuliaOpt/GLPKMathProgInterface.jl>`_           |  ``GLPKSolver[LP|MIP]()``   |  GPL        | X  |      |  X   |     |       |     |
+----------------------------------------------------------------------------------+---------------------------------------------------------------------------------+-----------------------------+-------------+----+------+------+-----+-------+-----+
| `Gurobi <http://gurobi.com>`_                                                    | `Gurobi.jl <https://github.com/JuliaOpt/Gurobi.jl>`_                            | ``GurobiSolver()``          |  Comm.      | X  |   X  |  X   |     |       |     |
+----------------------------------------------------------------------------------+---------------------------------------------------------------------------------+-----------------------------+-------------+----+------+------+-----+-------+-----+
| `Ipopt <https://projects.coin-or.org/Ipopt>`_                                    | `Ipopt.jl <https://github.com/JuliaOpt/Ipopt.jl>`_                              | ``IpoptSolver()``           |  EPL        | X  |      |      |  X  |       |     |
+----------------------------------------------------------------------------------+---------------------------------------------------------------------------------+-----------------------------+-------------+----+------+------+-----+-------+-----+
| `MOSEK <http://www.mosek.com/>`_                                                 | `Mosek.jl <https://github.com/JuliaOpt/Mosek.jl>`_                              | ``MosekSolver()``           |  Comm.      | X  |   X  |  X   |  X  |       | X   |
+----------------------------------------------------------------------------------+---------------------------------------------------------------------------------+-----------------------------+-------------+----+------+------+-----+-------+-----+
| `NLopt <http://ab-initio.mit.edu/wiki/index.php/NLopt>`_                         | `NLopt.jl <https://github.com/JuliaOpt/NLopt.jl>`_                              | ``NLoptSolver()``           |  LGPL       |    |      |      |  X  |       |     |
+----------------------------------------------------------------------------------+---------------------------------------------------------------------------------+-----------------------------+-------------+----+------+------+-----+-------+-----+
| `SCS <https://github.com/cvxgrp/scs>`_                                           | `SCS.jl <https://github.com/JuliaOpt/SCS.jl>`_                                  |  ``SCSSolver()``            |  MIT        | X  |  X   |      |     |       | X   |
+----------------------------------------------------------------------------------+---------------------------------------------------------------------------------+-----------------------------+-------------+----+------+------+-----+-------+-----+

Where:

- LP = Linear programming
- SOCP = Second-order conic programming (including problems with convex quadratic constraints and/or objective)
- MILP = Mixed-integer linear programming
- NLP = Nonlinear programming
- MINLP = Mixed-integer nonlinear programming
- SDP = Semidefinite programming
`*` = ``AmplNLWriter(CoinOptServices.bonmin)``
`**` = ``AmplNLWriter(CoinOptServices.couenne)``

To install Gurobi, for example, and use it with a JuMP model ``m``, run::

    Pkg.add("Gurobi")
    using JuMP
    using Gurobi

    m = Model(solver=GurobiSolver())

Setting solver options is discussed in the :ref:`Model <ref-model>` section.

Solver-specific notes follow below.

Artelys Knitro
++++++++++++++

Requires a license. The KNITRO.jl interface currently supports only nonlinear problems.

BARON
+++++

Requires a license. A trial version is available for small problem instances.

COIN-OR Clp and Cbc
+++++++++++++++++++

Binaries for Clp and Cbc are provided on OS X and Windows (32- and 64-bit) by default. On Linux, they will be compiled from source (be sure to have a C++ compiler installed). Cbc supports "SOS" constraints but does *not* support MIP callbacks.


CPLEX
+++++

Requires a working installation of CPLEX with a license (free for faculty members and graduate teaching assistants). The interface requires using CPLEX as a shared library, which is unsupported by the CPLEX developers. Special installation steps are required on OS X. CPLEX supports MIP callbacks and "SOS" constraints.


ECOS
++++

ECOS can be used by JuMP to solve LPs and SOCPs. ECOS does not support general quadratic objectives or constraints, only second-order conic constraints specified by using ``norm`` or the quadratic form ``x'x <= y^2``.

FICO Xpress
+++++++++++

Requires a working installation of Xpress with an active license (it is possible to get license for academic use, see `FICO Academic Partner Program <http://subscribe.fico.com/Academic-Partner-Program>`_). Supports SOCP and "SOS" constraints. Callbacks are not yet supported.

.. warning::
   If you are using 64-bit Xpress, you must use 64-bit Julia (and similarly with 32-bit Xpress).

.. warning::
  ``XpressSolver`` is not exported by default, thus one must use  ``Xpress.XpressSolver`` or import it explicitly.

GLPK
++++

GLPK binaries are provided on OS X and Windows (32- and 64-bit) by default. On Linux, it will be compiled from source. Note that ``GLPKSolverLP`` should be used for continuous problems and ``GLPKSolverMIP`` for problems with integer variables. GLPK supports MIP callbacks but does not support "SOS" constraints.

Gurobi
++++++

Requires a working installation of Gurobi with an activated license (free for academic use). Gurobi supports MIP callbacks and "SOS" constraints.

.. warning::
   If you are using 64-bit Gurobi, you must use 64-bit Julia (and similarly with 32-bit Gurobi).

Ipopt
+++++

Ipopt binaries are provided on OS X and Windows (32- and 64-bit) by default. On Linux, it will be compiled from source.
The default installation of Ipopt uses the open-source MUMPS library for sparse linear algebra.
Significant speedups can be obtained by manually compiling Ipopt to use proprietary sparse linear algebra libraries instead.
Julia can be pointed to use a custom version of Ipopt; we suggest posting to the `julia-opt <https://groups.google.com/forum/#!forum/julia-opt>`_ mailing list with your platform details for guidance on how to do this.


MOSEK
+++++

Requires a license (free for academic use). Mosek does not support the MIP callbacks used in JuMP.
For nonlinear optimization, Mosek supports only convex problems.
The Mosek interface is maintained by the Mosek team. (Thanks!)

NLopt
+++++

NLopt supports only nonlinear models. An algorithm must be specified as an option when using ``NLoptSolver``. NLopt is not recommended for large-scale models, because it does not currently exploit sparsity of derivative matrices.

SCS
+++

SCS can be used by JuMP to solve LPs and SOCPs, and SDPs. SCS is a first order solver and has low accuracy (:math:`10^{-4}`) by default; see the SCS.jl documentation for more information.

COIN-OR Bonmin and Couenne
++++++++++++++++++++++++++

Binaries of Bonmin and Couenne are provided on OS X and Windows (32- and 64-bit) by the [CoinOptServices.jl](https://github.com/JuliaOpt/CoinOptServices.jl) package. On Linux, they will be compiled from source. Once installed, they can be called either via `.osil` files using `OsilBonminSolver` and `OsilCouenneSolver` from [CoinOptServices.jl](https://github.com/JuliaOpt/CoinOptServices.jl), or via `.nl` files using `AmplNLWriter(CoinOptServices.bonmin)` and `AmplNLWriter(CoinOptServices.couenne)` from [AmplNLWriter.jl](https://github.com/JuliaOpt/AmplNLWriter.jl).
We recommend using the ``.nl`` format option, which is currently more stable and has better performance for derivative computations.
Since both Bonmin and Couenne use Ipopt for continuous subproblems, the same MUMPS sparse linear algebra performance caveat applies.

Other AMPL-compatible solvers
+++++++++++++++++++++++++++++

Any other solver not listed above that can be called from `AMPL <http://ampl.com/products/solvers/all-solvers-for-ampl/>`_ can be used by JuMP through the
`AmplNLWriter.jl <https://github.com/JuliaOpt/AmplNLWriter.jl>`_ package. The first argument to ``AmplNLSolver``
can be used to specify a solver executable name.

For example, `SCIP <http://scip.zib.de/>`_ is a powerful noncommercial mixed-integer programming solver. To use SCIP within JuMP, you must first download and `compile SCIP with support for AMPL <http://zverovich.net/2012/08/07/using-scip-with-ampl.html>`_. Then you may use ``AmplNLSolver("/path/to/scipampl")`` where ``scipampl`` is the executable produced from the compilation process.

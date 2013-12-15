.. _jump-installation:

------------------
Installation Guide
------------------

This guide will briefly guide you through installing Julia, JuMP and[a] solver[s] of your choice.

Getting Julia
^^^^^^^^^^^^^

At the time of writing this documentation the latest release of Julia is version ``0.2``, which is the version required by JuMP. You can easily build from source on OSX and Linux, but the binaries will work well for most people.

Download links and more detailed instructions are available on the `Julia website <http://julialang.org>`_.

Getting JuMP
^^^^^^^^^^^^

Once you've installed Julia, installing JuMP is simple. Julia has a git-based package system. To use it, open Julia in interactive mode (i.e. ``julia`` at the command line) and use the package manager::

    julia> Pkg.add("JuMP")

This command checks `METADATA.jl <https://github.com/JuliaLang/METADATA.jl/tree/devel>`_ to determine what the most recent version of JuMP is and then downloads it from its repository on GitHub.

Getting Solvers
^^^^^^^^^^^^^^^

Solver support in Julia is currently provided by writing a solver-specific package that provides a very thin wrapper around the solver's C interface and providing a standard interface that JuMP can call. If you are interested in providing an interface to your solver, please get in touch. We currently have interfaces for COIN-OR, Gurobi, GNU GLPK, and CPLEX.

COIN-OR Clp and Cbc
+++++++++++++++++++

Support for Clp and Cbc is provided via `CoinMP <https://projects.coin-or.org/CoinMP>`_ and the `Cbc.jl <https://github.com/mlubin/Cbc.jl>`_ and `Clp.jl <https://github.com/mlubin/Clp.jl>`_ packages. You can install these solvers through the package manager::

    julia> Pkg.add("Cbc")
    julia> Pkg.add("Clp")  # Clp depends on Cbc, so installing Clp first
                           # will install both.

Regarding the CoinMP binary itself, installation will differ by platform:

* Linux - Only option is to build from source, which will happen automatically.
* OS X - Downloads binary via the `Homebrew.jl <https://github.com/staticfloat/Homebrew.jl>`_ package.
* Windows - **Only 32-bit versions of Julia are supported by the COIN solvers at this time**. The 32-bit version of Julia can be used on 64-bit Windows with no issues. Binary download. Will require `Visual Studio 2012 redistributable <http://www.microsoft.com/en-us/download/details.aspx?id=30679>`_ if not already installed.

Clp and Cbc, if available, are the default choice of solver in JuMP. 

Gurobi
++++++

`Gurobi <http://gurobi.com>`_ is an excellent high-performance commercial solver. It supports quadratic objectives and constraints, and is currently the only solver supported by Julia/JuMP with that functionality. Install Gurobi as you normally would and then add the `Gurobi.jl <https://github.com/lindahua/Gurobi.jl>`_ package::

    julia> Pkg.add("Gurobi")

.. warning::
   If you are using 64-bit Gurobi, you must use 64-bit Julia (and similarly with 32-bit Gurobi).
  
The Gurobi package README contains examples of how to use Gurobi within JuMP.

GLPK
++++

Installing `GLPK <https://github.com/carlobaldassi/GLPK.jl>`_ is a bit more involved than can be covered here - see the `documentation <https://gplkjl.readthedocs.org/en/latest/glpk.html>`_ for more information.

CPLEX
+++++

`CPLEX <http://www-01.ibm.com/software/commerce/optimization/cplex-optimizer/>` is a leading commercial solver. An experimental interface is available via the `CPLEXLink <https://github.com/joehuchette/CPLEXLink.jl>` package. Note that this interface requires using CPLEX as a shared library, which is unsupported by the CPLEX developers.




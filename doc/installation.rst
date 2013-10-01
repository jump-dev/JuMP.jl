.. _jump-installation:

------------------
Installation Guide
------------------

This guide will briefly guide you through installing Julia, JuMP and[a] solver[s] of your choice.

Getting Julia
^^^^^^^^^^^^^

At the time of writing this documention Julia is just about to release version ``0.2``. You should download this version instead of the previous version (``0.1.2``) which is not supported by JuMP or the majority of packages. You can easily build from source on OSX and Linux, but the binaries will work well for most people.

Download links and more detailed instructions are available on the `Julia website <http://julialang.org>`_.

Getting JuMP
^^^^^^^^^^^^

Once you've installed Julia, installing JuMP is simple. Julia has a git-based package system. To use it, open Julia in interactive mode (i.e. ``julia`` at the command line) and use the package manager::

    julia> Pkg.add("JuMP")

This command checks `METADATA.jl <https://github.com/JuliaLang/METADATA.jl/tree/devel>`_ to determine what the most recent version of JuMP is and then downloads it from its repository on GitHub.

Getting Solvers
^^^^^^^^^^^^^^^

Solver support in Julia is currently provided by writing a solver-specific package that provides a very thin wrapper around the solver's C interface and providing a standard interface that JuMP can call. If you are interested in providing an interface to your solver, please get in touch. We currently have interfaces for COIN-OR, Gurobi, and GNU GLPK.

COIN-OR Clp and Cbc
+++++++++++++++++++

Support for Clp and Cbc is provided via `CoinMP <https://projects.coin-or.org/CoinMP>`_ and the `Cbc.jl <https://github.com/mlubin/Cbc.jl>`_ and `Clp.jl <https://github.com/mlubin/Clp.jl>`_ packages. You can install these solvers through the package manager::

    julia> Pkg.add("Cbc")
    julia> Pkg.add("Clp")  # Clp depends on Cbc, so installing Clp first
                           # will install both.

Regarding the CoinMP binary itself, installation will differ by platform:

* Linux - Only option is to build from source, which will happen automatically.
* OSX - Downloads binary via the `Homebrew.jl <https://github.com/staticfloat/Homebrew.jl>`_ package.
* Windows - **Only 32-bit versions of Julia are supported by the COIN solvers at this time**. The 32-bit version of Julia can be used on 64-bit Windows with no issues. Binary download. Will require `Visual Studio 2012 redistributable <http://www.microsoft.com/en-us/download/details.aspx?id=30679>`_ if not already installed.

Clp and Cbc, if available, are the default choice of solver in JuMP. To manually select Clp and/or Cbc as your LP/MIP solver, at the top of your code include the following::

    using Cbc  # For MIPs
    using Clp  # For LPs

and at some point, before you call ``solve(model)``::

    setMIPSolver(:Cbc)  # For MIPs
    setLPSolver(:Clp)   # For LPs

Gurobi
++++++

`Gurobi <http://gurobi.com>`_ is an excellent high-performance commerical solver. It supports quadratic objectives and constraints, and is currently the only solver supported by Julia/JuMP with that functionality. Install Gurobi as you normally would and then add the `Gurobi.jl <https://github.com/lindahua/Gurobi.jl>`_ package::

    julia> Pkg.add("Gurobi")

.. note::
   Previously versions of Gurobi.jl required adding environmental variables but it should now "just work" out-of-the-box.

.. warning::
   If you are using 64-bit Gurobi, you must use 64-bit Julia (and similarily with 32-bit Gurobi).
   
To use Gurobi in your Julia programs::

    using Gurobi
    # ... code...
    # ... building model ...
    setMIPSolver(:Gurobi)  # For MIPs
    setLPSolver(:Gurobi)  # For LPs
    # Gurobi will automatically be used for QCQPs if available.

GLPK
++++

Installing `GLPK <https://github.com/carlobaldassi/GLPK.jl>`_ is a bit more involved than can be covered here - see the `documentation <https://gplkjl.readthedocs.org/en/latest/glpk.html>`_ for more information.




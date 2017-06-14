============================================
JuMP --- Julia for Mathematical Optimization
============================================

.. include:: warn.rst

.. module:: JuMP
   :synopsis: Julia for Mathematical Optimization

`JuMP <https://github.com/JuliaOpt/JuMP.jl>`_ is a domain-specific modeling language for
`mathematical optimization <http://en.wikipedia.org/wiki/Mathematical_optimization>`_
embedded in `Julia <http://julialang.org/>`_.
It currently supports a number of open-source and commercial solvers (see below)
for a variety of problem classes, including **linear programming**, **mixed-integer programming**, **second-order conic programming**, **semidefinite programming**, and **nonlinear programming**.
JuMP's features include:

* User friendliness

  * Syntax that mimics natural mathematical expressions.
  * Complete documentation.

* Speed

  * Benchmarking has shown that JuMP can create problems at similar speeds to
    special-purpose modeling languages such as `AMPL <http://www.ampl.com/>`_.
  * JuMP communicates with solvers in memory, avoiding the need to write
    intermediary files.

* Solver independence

  * JuMP uses a generic solver-independent interface provided by the
    `MathProgBase <https://github.com/mlubin/MathProgBase.jl>`_ package, making it easy
    to change between a number of open-source and commercial optimization software packages ("solvers").
  * Currently supported solvers include
    `Artelys Knitro <http://artelys.com/en/optimization-tools/knitro>`_,
    `Bonmin <https://projects.coin-or.org/Bonmin>`_,
    `Cbc <https://projects.coin-or.org/Cbc>`_,
    `Clp <https://projects.coin-or.org/Clp>`_,
    `Couenne <https://projects.coin-or.org/Couenne>`_,
    `CPLEX <http://www-01.ibm.com/software/commerce/optimization/cplex-optimizer/>`_,
    `ECOS <https://github.com/ifa-ethz/ecos>`_,
    `FICO Xpress <http://www.fico.com/en/products/fico-xpress-optimization-suite>`_,
    `GLPK <http://www.gnu.org/software/glpk/>`_,
    `Gurobi <http://www.gurobi.com>`_,
    `Ipopt <https://projects.coin-or.org/Ipopt>`_,
    `MOSEK <http://www.mosek.com/>`_,
    `NLopt <http://ab-initio.mit.edu/wiki/index.php/NLopt>`_, and
    `SCS <https://github.com/cvxgrp/scs>`_.

* Access to advanced algorithmic techniques

  * Including :ref:`efficient LP re-solves <probmod>` and :ref:`callbacks for mixed-integer programming <callbacks>` which previously required using solver-specific and/or low-level C++ libraries.

* Ease of embedding

  * JuMP itself is written purely in Julia. Solvers are the only binary dependencies.
  * Being embedded in a general-purpose programming language makes it easy to solve optimization problems as part of a larger workflow (e.g., inside a simulation, behind a web server, or as a subproblem in a decomposition algorithm).

    * As a trade-off, JuMP's syntax is constrained by the syntax available in Julia.

  * JuMP is `MPL <https://www.mozilla.org/MPL/2.0/>`_ licensed, meaning that it can be embedded in commercial software that complies with the terms of the license.

While neither Julia nor JuMP have reached version 1.0 yet, the releases are stable enough for everyday use and are being used in a number of research projects and neat applications by a growing community of users who are early adopters. JuMP remains under active development, and we welcome your feedback, suggestions, and bug reports.

Installing JuMP
---------------

If you are familiar with Julia you can get started quickly by using the
package manager to install JuMP::

    julia> Pkg.add("JuMP")

And a solver, e.g.::

    julia> Pkg.add("Clp")  # Will install Cbc as well

Then read the :ref:`quick-start` and/or see a :ref:`simple-example`.
The subsequent sections detail the complete functionality of JuMP.

Contents
--------

.. toctree::
    :maxdepth: 2

    installation.rst
    quickstart.rst
    refmodel.rst
    refvariable.rst
    refexpr.rst
    probmod.rst
    callbacks.rst
    nlp.rst

-----------
Citing JuMP
-----------



If you find JuMP useful in your work, we kindly request that you cite the following `paper <http://dx.doi.org/10.1137/15M1020575>`_:

.. code-block:: none

    @article{DunningHuchetteLubin2017,
    author = {Iain Dunning and Joey Huchette and Miles Lubin},
    title = {JuMP: A Modeling Language for Mathematical Optimization},
    journal = {SIAM Review},
    volume = {59},
    number = {2},
    pages = {295-320},
    year = {2017},
    doi = {10.1137/15M1020575},
    }

A preprint of this paper is freely available on `arXiv <http://arxiv.org/abs/1508.01982>`_.

For an earlier work where we presented a prototype implementation of JuMP, see
`here <http://dx.doi.org/10.1287/ijoc.2014.0623>`_:

.. code-block:: none

    @article{LubinDunningIJOC,
    author = {Miles Lubin and Iain Dunning},
    title = {Computing in Operations Research Using Julia},
    journal = {INFORMS Journal on Computing},
    volume = {27},
    number = {2},
    pages = {238-248},
    year = {2015},
    doi = {10.1287/ijoc.2014.0623},
    }

A preprint of this paper is also `freely available <http://arxiv.org/abs/1312.1431>`_.

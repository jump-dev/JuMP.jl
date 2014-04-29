===========================================
JuMP --- Julia for Mathematical Programming
===========================================

.. module:: JuMP
   :synopsis: Julia for Mathematical Programming

`JuMP <https://github.com/JuliaOpt/JuMP.jl>`_ is a domain-specific modeling language for 
`mathematical programming <http://en.wikipedia.org/wiki/Mathematical_optimization>`_ 
embedded in `Julia <http://julialang.org/>`_. It currently supports a number
of open-source and commercial solvers (`Clp <https://projects.coin-or.org/Clp>`_,
`Cbc <https://projects.coin-or.org/Cbc>`_, `GLPK <http://www.gnu.org/software/glpk/>`_, `Gurobi <http://www.gurobi.com>`_, `MOSEK <http://www.mosek.com/>`_, and `CPLEX <http://www-01.ibm.com/software/commerce/optimization/cplex-optimizer/>`_) via a generic solver-independent 
interface provided by the `MathProgBase <https://github.com/mlubin/MathProgBase.jl>`_
package. 

One the best features of JuMP is its **speed** - benchmarking has shown
that it can create problems at similar speeds to special-purpose modeling
languages such as `AMPL <http://www.ampl.com/>`_ while maintaining the expressiveness
of a generic high-level programming language. JuMP communicates with solvers in-memory, 
avoiding the need to write intermediary files and enabling access to **advanced
features** such as :ref:`efficient LP re-solves <probmod>` and :ref:`callbacks for mixed-integer programming <callbacks>`.

JuMP has recently enabled support for nonlinear programming for functions that can be expressed in closed algebraic form. JuMP computes exact sparse second-order derivatives needed by efficient interior-point solvers.

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
    nlp.rst
    probmod.rst
    callbacks.rst

-----------
Citing JuMP
-----------

Further discussion of the design of JuMP in the context of existing domain-specific languages for mathematical programming, together with extensive benchmarks, is given in [1]_. If you find JuMP useful in your work, we request that you cite this paper.


.. [1] M. Lubin and I. Dunning, "Computing in Operations Research using Julia", 2013. `arXiv:1312.1431 <http://arxiv.org/abs/1312.1431>`_

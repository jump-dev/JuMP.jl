===========================================
JuMP --- Julia for Mathematical Programming
===========================================

.. module:: JuMP
   :synopsis: Julia for Mathematical Programming

JuMP is a domain-specific modeling language for 
`mathematical programming <http://en.wikipedia.org/wiki/Mathematical_optimization>`_ 
embedded in `Julia <http://julialang.org/>`_. It currently supports a number
of open-source and commerical solvers (`COIN Clp <https://projects.coin-or.org/Clp>`_,
`COIN Cbc <https://projects.coin-or.org/Cbc>`_, `GNU GLPK <http://www.gnu.org/software/glpk/>`_,
and `Gurobi <http://www.gurobi.com>`_) via a generic solver-independent 
interface provided by the `MathProgBase <https://github.com/mlubin/MathProgBase.jl>`_
package. One the best features of JuMP is its speed - benchmarking has shown
that it can create problems at similar speeds to special-purpose modeling
languages such as `AMPL <http://www.ampl.com/>`_ while maintaing the expressiveness
of a generic high-level programming language.

If you are familiar with Julia you can get started quickly by using the
package manager to install JuMP::

    julia> Pkg.add("JuMP")

And a solver, e.g.::

    julia> Pkg.add("Clp")  # Will install Cbc as well

Then read the :ref:`quick-start` and/or see a :ref:`simple-example`. We also
have details of the functions and types defined by JuMP.
If you are new to Julia or want more details, read on to the next section.

.. Full installation guide, with solvers
.. include:: installation.rst

.. Outlines key functionality of JuMP
.. include:: quickstart.rst

.. Beginning of Reference section
.. First up, Model
.. include:: refmodel.rst
.. Variable
.. include:: refvariable.rst
.. Expressions, Constraints
.. include:: refexpr.rst

----------
References
----------

Further discussion of the design of JuMP in the context of existing domain-specific languages for mathematical programming, together with extensive benchmarks, is given in [1]_.   


.. [1] M. Lubin and I. Dunning, "Computing in Operations Research using Julia", May 2013, under review. `Preprint <http://www.optimization-online.org/DB_HTML/2013/05/3883.html>`_

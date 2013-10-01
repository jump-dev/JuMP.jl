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

.. Walks through a simple example
.. include:: example.rst

.. Beginning of Reference section
.. First up, Model
.. include:: refmodel.rst
.. Variable
.. include:: refvariable.rst


.. _jump-function-list:

Quadratic Objectives
^^^^^^^^^^^^^^^^^^^^

There is preliminary support for convex quadratic objectives. Currently the
only supported solver is ``Gurobi``; it must be set as the ``lpsolver`` or 
``mipsolver`` when solving QPs or mixed-integer QPs, respectively. The 
``@setObjective`` macro does not yet support quadratic terms, but you may
use instead the (slower) operator overloading functionality and the 
``setObjective`` function::

    MathProgBase.setlpsolver(:Gurobi)
    m = Model(:Min)
    @defVar(m, 0 <= x <= 2 )
    @defVar(m, 0 <= y <= 30 )

    setObjective(m, x*x+ 2x*y + y*y )
    @addConstraint(m, x + y >= 1 )
      
    print(m)

    status = solve(m)

Quadratic Constraints
^^^^^^^^^^^^^^^^^^^^^

There is preliminary support for convex quadratic constraints. Currently the 
only supported solver is ``Gurobi``; it must be set as the ``lpsolver`` or
``mipsolver`` when solving QC programs. The ``@addConstraint`` macro does not 
yet support quadratic expressions, but you may instead use the (slower) 
operator overloading functionality via the ``addConstraint`` function::

    MathProgBase.setlpsolver(:Gurobi)
    m = Model(:Min)
    @defVar(m, -1 <= x <= 1)
    @defVar(m, -1 <= y <= 1)

    @setObjective(m, x + y)
    addConstraint(m, x*x + y*y <= 1)

    print(m)

    status = solve(m)

    
----------------------------------------------------
Expressions, constraints, and the objective function
----------------------------------------------------

Macros vs operator overloading @addConstraint duals @setObjective



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

Then read the :ref:`quick-start`. If you are new to Julia or want more details,
read on.

.. Full installation guide, with solvers
.. include:: installation.rst

.. include:: quickstart.rst

----------------------------------------------------
Expressions, constraints, and the objective function
----------------------------------------------------

Macros vs operator overloading @addConstraint duals @setObjective


.. _jump-variables:

---------
Variables
---------


Variables, also known as columns or decision variables, are the results of the optimization.

Creation
^^^^^^^^

The primary way to create variables is with the ``@defVar`` macro. The first argument will always be a ``Model``. In the examples below we assume ``m`` is already defined. The second argument is an expression that declares the variable name and optionally allows specification of lower and upper bounds. For example::

    @defVar(m, x )              # No bounds
    @defVar(m, x >= lb )        # Lower bound only (note: 'lb <= x' is not valid)
    @defVar(m, x <= ub )        # Upper bound only
    @defVar(m, lb <= x <= ub )  # Lower and upper bounds

All these variations create a new local variable, in this case ``x``. Integer and binary restrictions can optionally be specified with a third argument, ``Int`` or ``Bin``.

Arrays of variables
^^^^^^^^^^^^^^^^^^^

To create arrays of variables we append brackets to the variable name. For example::

    @defVar(m, x[1:M,1:N] >= 0 )

will create an ``M`` by ``N`` array of variables. Both ranges and arbitrary iterable sets are supported as index sets. Currently we only support ranges of the form ``a:b`` where ``a`` is an explicit integer, not a variable. Using ranges will generally be faster than using arbitrary symbols. You can mix both ranges and lists of symbols, as in the following example::

    s = ["Green","Blue"]
    @defVar(m, x[-10:10,s] , Int)
    x[-4,"Green"]

An interesting feature is that bounds can depend on variable indices::

    @defVar(m, x[i=1:10] >= i )

Variables can be constructed manually, one-by-one, using

::

    x = Variable(model, lower, upper, category)

but this is not considered idiomatic.

Modification
^^^^^^^^^^^^

Bounds
++++++
* ``setLower(x, lower)``, ``getLower(x)`` - Set/get the lower bound of a variable.
* ``setUpper(x, upper)``, ``getUpper(x)`` - Set/get the upper bound of a variable.

Values at solution
++++++++++++++++++
* ``getValue(x)`` - Get the value of this variable in the solution.
* ``getDual(x)`` - Get the reduced cost of this variable in the solution.

Names
+++++
Variables can have internal names that can be used for writing models to file. This is currently not done for performance reasons, but may be added if there is demand.
* ``setName(x, newName)``, ``getName(x)`` - Set/get the variable's internal name.



.. _jump-function-list:

--------------------------------
List of functions in JuMP module
--------------------------------

.. function:: setName(x, newName)

    Assigns a name to the the variable ``x``.



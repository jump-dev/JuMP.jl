.. _quick-start:


.. include:: warn.rst
-----------------
Quick Start Guide
-----------------

This quick start guide will introduce the main concepts of JuMP.
If you are familiar with another modeling language embedded in a high-level
language such as PuLP (Python) or a solver-specific interface you will find
most of this familiar, with the exception of *macros*. A deep understanding
of macros is not essential, but if you would like to know more please see
the `Julia documentation <http://docs.julialang.org/en/latest/manual/metaprogramming/>`_.
If you are coming from an AMPL or similar background, you may find some of
the concepts novel but the general appearance will still be familiar.


Creating a Model
^^^^^^^^^^^^^^^^

**Models** are Julia objects. They are created by calling the constructor::

    m = Model()

All variables and constraints are associated with a ``Model`` object. For
a list of all functions related to ``Model``, including how to change the
default solver and set solver parameters, see :ref:`ref-model`.


Defining Variables
^^^^^^^^^^^^^^^^^^

**Variables** are also Julia objects, and are defined using the ``@variable``
macro. The first argument will always be the ``Model`` to associate this
variable with. In the examples below we assume ``m`` is already defined.
The second argument is an expression that declares the variable name and
optionally allows specification of lower and upper bounds. For example::

    @variable(m, x )              # No bounds
    @variable(m, x >= lb )        # Lower bound only (note: 'lb <= x' is not valid)
    @variable(m, x <= ub )        # Upper bound only
    @variable(m, lb <= x <= ub )  # Lower and upper bounds

All these variations introduce a new variable ``x`` in the local scope.
The names of your variables must be valid Julia variable names.
For information about common operations on variables, e.g. changing their
bounds, see the :ref:`ref-variable` section.

**Integer** and **binary** restrictions can optionally be specified with a
third argument, ``Int`` or ``Bin``.

To create arrays of variables we append brackets to the variable name.
For example::

    @variable(m, x[1:M,1:N] >= 0 )

will create an ``M`` by ``N`` array of variables. Both ranges and arbitrary
iterable sets are supported as index sets. Currently we only support ranges
of the form ``a:b`` where ``a`` is an explicit integer, not a variable.
Using ranges will generally be faster than using arbitrary symbols. You can
mix both ranges and lists of symbols, as in the following example::

    s = ["Green", "Blue"]
    @variable(m, x[-10:10,s], Int )
    # e.g. x[-4, "Green"]

Finally, bounds can depend on variable indices::

    @variable(m, x[i=1:10] >= i )


Objective and Constraints
^^^^^^^^^^^^^^^^^^^^^^^^^

JuMP allows users to use a natural notation to describe linear expressions. To add constraints, use the ``@constraint()`` and ``@objective()``
macros, e.g.::

    @constraint(m, x[i] - s[i] <= 0)  # Other options: == and >=
    @constraint(m, sum(x[i] for i=1:numLocation) == 1)
    @objective(m, Max, 5x + 22y + (x+y)/2) # or Min

.. note::
    The ``sense`` passed to ``@objective`` must be a `symbol <http://docs.julialang.org/en/latest/manual/metaprogramming/#symbols>`_ type: ``:Min`` or ``:Max``, although the macro accepts ``:Min`` and ``:Max``, as well as ``Min`` and ``Max`` (without the colon) directly.

The ``sum()`` syntax directly follows Julia's own generator expression syntax. You may use conditions within sums, e.g.::

    sum(expression for i = I1, j = I2 if cond)

which is equivalent to::

    a = zero(AffExpr)
    for i = I1
        for j = I2
            ...
            if cond
                a += expression
            end
            ...
        end
    end

.. note::
    JuMP previously used a special curly brace syntax for ``sum{}``, ``prod{}``, and ``norm2{}``. This has been entirely replaced by ``sum()``, ``prod()``, and ``norm()`` since Julia 0.5. The curly brace syntax is deprecated and will be removed in a future release.


.. Walks through a simple example
.. include:: example.rst

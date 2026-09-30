.. _quick-start:

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

**Variables** are also Julia objects, and are defined using the ``@defVar``
macro. The first argument will always be the ``Model`` to associate this
variable with. In the examples below we assume ``m`` is already defined.
The second argument is an expression that declares the variable name and
optionally allows specification of lower and upper bounds. For example::

    @defVar(m, x )              # No bounds
    @defVar(m, x >= lb )        # Lower bound only (note: 'lb <= x' is not valid)
    @defVar(m, x <= ub )        # Upper bound only
    @defVar(m, lb <= x <= ub )  # Lower and upper bounds

All these variations introduce a new variable ``x`` in the local scope.
The names of your variables must be valid Julia variable names.
For information about common operations on variables, e.g. changing their
bounds, see the :ref:`ref-variable` section.

**Integer** and **binary** restrictions can optionally be specified with a
third argument, ``Int`` or ``Bin``.

To create arrays of variables we append brackets to the variable name.
For example::

    @defVar(m, x[1:M,1:N] >= 0 )

will create an ``M`` by ``N`` array of variables. Both ranges and arbitrary
iterable sets are supported as index sets. Currently we only support ranges
of the form ``a:b`` where ``a`` is an explicit integer, not a variable.
Using ranges will generally be faster than using arbitrary symbols. You can
mix both ranges and lists of symbols, as in the following example::

    s = ["Green", "Blue"]
    @defVar(m, x[-10:10,s], Int )
    # e.g. x[-4, "Green"]

Finally, bounds can depend on variable indices::

    @defVar(m, x[i=1:10] >= i )


Objective and Constraints
^^^^^^^^^^^^^^^^^^^^^^^^^

JuMP allows users to use a natural notation to describe linear expressions.
There are two ways to do so. The first is very similar to other modeling
languages and has no restrictions. The second utilizes Julia's powerful
metaprogramming features to get excellent performance even for large problems,
but has some restrictions on how they can be used.

To add constraints in the first way, use the ``addConstraint()`` and ``setObjective()``
functions, e.g.::

    addConstraint(m, y + z == 4)  # Other options: <= and >=
    setObjective(m, :Max, 5x + 22y + (x+y)/2) # or :Min

The second way is visually very similar, and uses the ``@addConstraint`` and ``@setObjective``
macros, e.g.::

    @addConstraint(m, x[i] - s[i] <= 0)
    @setObjective(m, Max, sum{x[i], i=1:numLocation} )

.. note::
    The ``sense`` passed to ``setObjective`` must be a `symbol <http://docs.julialang.org/en/latest/manual/metaprogramming/#symbols>`_ type: ``:Min`` or ``:Max``.
    The ``@setObjective`` macro accepts ``:Min`` and ``:Max``, as well as ``Min`` and ``Max`` (without the colon) directly.

You may have noticed a special ``sum{}`` operator above. This is defined only for
the second kind of function. The syntax is of the form (where ``IX`` can be any iterable)::

	sum{expression, i = I1, j = I2, ...}

which is equivalent to::

    a = zero(AffExpr)  # Create a new empty affine expression
    for i = I1
        for j = I2
            ...
            a += expression
            ...
        end
    end

You can also put a condition in::

    sum{expression, i = I1, j = I2, ...; cond}

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


.. Walks through a simple example
.. include:: example.rst

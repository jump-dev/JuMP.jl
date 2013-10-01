.. _ref-variable:


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



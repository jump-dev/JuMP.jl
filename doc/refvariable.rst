.. _ref-variable:

---------
Variables
---------

Variables, also known as columns or decision variables, are the results of the optimization.

Constructors
^^^^^^^^^^^^

The primary way to create variables is with the ``@defVar`` macro.
The first argument will always be a ``Model``. In the examples below we assume
``m`` is already defined. The second argument is an expression that declares
the variable name and optionally allows specification of lower and upper bounds.
Adding variables "column-wise", e.g., as in column generation, is supported as well;
see the syntax discussed in the :ref:`probmod` section.

::

    @defVar(m, x )              # No bounds
    @defVar(m, x >= lb )        # Lower bound only (note: 'lb <= x' is not valid)
    @defVar(m, x <= ub )        # Upper bound only
    @defVar(m, lb <= x <= ub )  # Lower and upper bounds

All these variations create a new local variable, in this case ``x``. 
The names of your variables must be valid Julia variable names.
Integer and binary restrictions can optionally be specified with a third argument, ``Int`` or ``Bin``.
For advanced users, ``SemiCont`` and ``SemiInt`` may be used to create
`semicontinuous <http://orinanobworld.blogspot.com/2011/03/semicontinuous-variables.html>`_ or
`semi-integer <http://www.gams.com/mccarl/mccarlhtml/semi-integer_variables.htm>`_ variables,
respectively.

To create arrays of variables we append brackets to the variable name.

::

    @defVar(m, x[1:M,1:N] >= 0 )

will create an ``M`` by ``N`` array of variables. Both ranges and arbitrary
iterable sets are supported as index sets. Currently we only support ranges
of the form ``a:b`` where ``a`` is an explicit integer, not a variable. Using
ranges will generally be faster than using arbitrary symbols. You can mix both
ranges and lists of symbols, as in the following example::

    s = ["Green","Blue"]
    @defVar(m, x[-10:10,s] , Int)
    x[-4,"Green"]

Bounds can depend on variable indices::

    @defVar(m, x[i=1:10] >= i )

And indices can have dependencies on preceding indices (e.g. "triangular indexing")::

    @defVar(m, x[i=1:10;j=i:10] >= 0)

Note the dependency must be on preceding indices, going from left to right. That is,
``@defVar(m, x[i=j:10,i=1:10] >= 0)`` is not valid JuMP code.

Finally, variables can be constructed manually, one-by-one::

    x = Variable(m::Model, lower::Number, upper::Number, category::Symbol, name::String)
    x = Variable(m::Model, lower::Number, upper::Number, category::Symbol)

where ``category`` is one of ``:Cont``, ``:Int``, ``:Bin``, ``:SemiCont``, and ``:SemiInt``.
This form of constructing variables is not considered idiomatic JuMP code.

.. note::
    ``@defVar`` is equivalent to a simple assignment ``x = ...`` in Julia and therefore redefines variables without warning. The following code may lead to unexpected results::
    
    @defVar(m, x[1:10,1:10])
    @defVar(m, x[1:5])

    After the second line, the Julia variable ``x`` refers to a set of variables indexed
    by the range ``1:5``.
    The reference to the first set of variables has been lost, although they will remain
    in the model.

Methods
^^^^^^^

**Bounds**

* ``setLower(x::Variable, lower)``, ``getLower(x::Variable)`` - Set/get the lower bound of a variable.
* ``setUpper(x::Variable, upper)``, ``getUpper(x::Variable)`` - Set/get the upper bound of a variable.


**Helper functions**

* ``sum(x)`` - Operates on arrays of variables, efficiently produces an affine expression. Available in macros.
* ``dot(x, coeffs)`` - Performs a generalized "dot product" for arrays of variables and coefficients up to three dimensions, or equivalently the sum of the elements of the Hadamard product. Available in macros, and also as ``dot(coeffs, x)``.


**Values**

* ``getValue(x)`` - Get the value of this variable in the solution. If ``x`` is a single variable, this will simply return a number. 
  If ``x`` is indexable then it will return an indexable dictionary of values. When the model is unbounded, ``getValue`` will 
  instead return the corresponding components of an unbounded ray, if available from the solver.
* ``setValue(x,v)`` - Provide an initial value ``v`` for this variable that can be used by supporting MILP solvers. If ``v`` is ``NaN``, the solver may attempt to fill in this value to construct a feasible solution.
* ``getDual(x)`` - Get the reduced cost of this variable in the solution. Similar behavior to ``getValue`` for indexable variables.

.. note::
    The ``getValue`` function always returns a floating-point value, even when a variable is constrained to take integer values, as most solvers only guarantee integrality up to a particular numerical tolerance. The built-in ``iround`` function should be used to obtain integer values, e.g., by calling ``iround(getValue(x))``. 


**Names**

Variables (in the sense of columns) can have internal names (different from the Julia variable name) that can be used for writing models to file. This feature is disabled for performance reasons, but will be added if there is demand or a special use case.

* ``setName(x::Variable, newName)``, ``getName(x::Variable)`` - Set/get the variable's internal name.

.. _ref-model:

------
Models
------

Constructor
^^^^^^^^^^^

``Model`` is a type defined by JuMP. It has one constructor that takes one
argument, the objective sense. The objective sense can be either ``:Min``
or ``:Max``::

    m = Model(:Min)
    m = Model(:Max)

All variables and constraints are associated with a ``Model`` objects.

Methods
^^^^^^^

**General**

* ``getNumVars(m::Model)`` - returns the number of variables associated with the ``Model m``
* ``getNumConstraints(m::Model)`` - returns the number of constraints associated with the ``Model m``

**Objective**

* ``getObjective(m::Model)`` - returns the objective function as a ``QuadExpr``.
* ``setObjective(m::Model, a::AffExpr)``, ``setObjective(m::Model, q::QuadExpr)`` - sets the objective function to ``a`` and ``q`` respectively.
* ``getObjectiveSense(m::Model)`` - returns objective sense, either ``:Min`` or ``:Max``.
* ``setObjectiveSense(m::Model, newSense::Symbol)`` - sets the objective sense (``newSense`` is either ``:Min`` or ``:Max``).
* ``getObjectiveValue(m::Model)`` - returns objective value after a call to ``solve``.

**Output**

* ``writeLP(m::Model, filename::String)`` - write the model to ``filename`` in the LP file format.
* ``writeMPS(m::Model, filename::String)`` - write the model to ``filename`` in the MPS file format.


Quadratic Objectives
^^^^^^^^^^^^^^^^^^^^

Quadratic objectives are supported by JuMP but currently the only supported
solver is ``Gurobi``. The other issue is that the ``@setObjective`` macro
**does not yet support quadratic terms**, but you may use instead the (slower)
``setObjective`` function::

    m = Model(:Min)
    @defVar(m, 0 <= x <= 2 )
    @defVar(m, 0 <= y <= 30 )

    setObjective(m, x*x+ 2x*y + y*y )  # Cannot use macro
    @addConstraint(m, x + y >= 1 )
      
    print(m)

    status = solve(m)

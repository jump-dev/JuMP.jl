.. _ref-model:

------
Models
------

Constructor
^^^^^^^^^^^

``Model`` is a type defined by JuMP. All variables and constraints are 
associated with a ``Model`` object. It has a constructor that has no 
required arguments::

    m = Model()

The constructor also accepts an optional keyword argument, ``solver``,
which can be used to change the default solver behavior.

``solver`` must be an ``AbstractMathProgSolver`` object, which is constructed as follows::

    solver = solvername(Option1=Value1, Option2=Value2, ...)

where ``solvername`` is one of the supported solvers. See the :ref:`solver table <jump-solvertable>` for the list of available solvers and corresponding parameter names.  All options are solver-dependent; see corresponding solver packages for more information. 

.. note::
    Be sure that the solver provided supports the problem class of the model. For example ``ClpSolver`` and ``GLPKSolverLP`` support only linear programming problems. ``CbcSolver`` and ``GLPKSolverMIP`` support only mixed-integer programming problems.

As an example, we can create a ``Model`` object that will use GLPK's
exact solver for LPs as follows::
    
    m = Model(solver = GLPKSolverLP(method=:Exact))


Methods
^^^^^^^

**General**

* ``getNumVars(m::Model)`` - returns the number of variables associated with the ``Model m``.
* ``getNumConstraints(m::Model)`` - returns the number of constraints associated with the ``Model m``.

**Objective**

* ``getObjective(m::Model)`` - returns the objective function as a ``QuadExpr``.
* ``setObjective(m::Model, sense::Symbol, a::AffExpr)``, ``setObjective(m::Model, sense::Symbol, q::QuadExpr)`` - sets the objective function to ``a`` and ``q`` respectively, with given objective sense, which must be either ``:Min`` or ``:Max``.
* ``getObjectiveSense(m::Model)`` - returns objective sense, either ``:Min`` or ``:Max``.
* ``setObjectiveSense(m::Model, newSense::Symbol)`` - sets the objective sense (``newSense`` is either ``:Min`` or ``:Max``).
* ``getObjectiveValue(m::Model)`` - returns objective value after a call to ``solve``.

**Output**

* ``writeLP(m::Model, filename::String)`` - write the model to ``filename`` in the LP file format.
* ``writeMPS(m::Model, filename::String)`` - write the model to ``filename`` in the MPS file format.


Quadratic Objectives
^^^^^^^^^^^^^^^^^^^^

Quadratic objectives are supported by JuMP using a solver which implements the
corresponding extensions of the MathProgBase interface, but currently the 
``@setObjective`` macro **does not yet support quadratic terms**. Instead you
may use the (slower) ``setObjective`` function::

    m = Model()
    @defVar(m, 0 <= x <= 2 )
    @defVar(m, 0 <= y <= 30 )

    setObjective(m, :Min, x*x+ 2x*y + y*y )  # Cannot use macro
    @addConstraint(m, x + y >= 1 )
      
    print(m)

    status = solve(m)


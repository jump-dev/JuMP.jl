.. _ref-model:

------
Models
------

Constructor
^^^^^^^^^^^

``Model`` is a type defined by JuMP. All variables and constraints are 
associated with a ``Model`` object. It has a constructor that takes 
one required argument, the objective sense. The objective sense must be 
one of the symbols ``:Min`` or ``:Max``::

    m = Model(:Min)
    m = Model(:Max)

The constructor also accepts two optional keyword arguments, ``lpsolver``,
and ``mipsolver``, which can be used to change the default solver behavior.

``lpsolver`` must be an ``LPSolver`` object, which is constructed as follows::

    solver = LPSolver(solvername, Option1=Value1, Option2=Value2, ...)

where ``solvername`` is one of the suppored solvers (``:Clp``, ``:GLPK``, and ``:Gurobi``). All options are solver-dependent; see corresponding solver packages for more information. 

``mipsolver`` must be a ``MIPSolver`` object, which is built similarly to ``LPSolver``. The currently supported solvers are ``:Cbc``, ``:GLPK``, and ``:Gurobi``.

.. note::
    Currently, the ``mipsolver`` solver is used for any problem with integer variables present. The ``lpsolver`` solver is used for all other problems, including those with continuous variables and **quadratic objectives and/or constraints**.

As an example, we can create a ``Model`` object that will use GLPK's
exact solver for LPs as follows::
    
    m = Model(:Min, lpsolver = LPSolver(:GLPK, GLPKmethod=:Exact))


Methods
^^^^^^^

**General**

* ``getNumVars(m::Model)`` - returns the number of variables associated with the ``Model m``.
* ``getNumConstraints(m::Model)`` - returns the number of constraints associated with the ``Model m``.

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

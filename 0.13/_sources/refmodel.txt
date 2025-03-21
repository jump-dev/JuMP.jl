.. _ref-model:

.. include:: warn.rst

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

* ``MathProgBase.numvar(m::Model)`` - returns the number of variables associated with the ``Model m``.
* ``MathProgBase.numlinconstr(m::Model)`` - returns the number of linear constraints associated with the ``Model m``.
* ``MathProgBase.numquadconstr(m::Model)`` - returns the number of quadratic constraints associated with the ``Model m``.
* ``MathProgBase.numconstr(m::Model)`` - returns the total number of constraints associated with the ``Model m``.
* ``getsolvetime(m::Model)`` - returns the solve time reported by the solver if it is implemented.
* ``getnodecount(m::Model)`` - returns the number of explored branch-and-bound nodes, if it is implemented.
* ``getobjbound(m::Model)`` - returns the best known bound on the optimal objective value. This is used, for example, when a branch-and-bound method is stopped before finishing.
* ``getobjgap(m::Model)`` - returns the final relative optimality gap as optimization terminated. That is, it returns :math:`\frac{|b-f|}{|f|}`, where :math:`b` is the best bound and :math:`f` is the best feasible objective value.
* ``getrawsolver(m::Model)`` - returns an object that may be used to access a solver-specific API.
* ``getsimplexiter(m::Model)`` - returns the cumulative number of simplex iterations during the optimization process. In particular, for a MIP it returns the total simplex iterations for all nodes.
* ``getbarrieriter(m::Model)`` - returns the cumulative number of barrier iterations during the optimization process.
* ``internalmodel(m::Model)`` - returns the internal low-level ``AbstractMathProgModel`` object which can be used to access any functionality that is not exposed by JuMP. See the MathProgBase `documentation <https://mathprogbasejl.readthedocs.org/en/latest/lowlevel.html>`_.
* ``solve(m::Model; suppress_warnings=false, relaxation=false)`` - solves the model using the selected solver (or a default for the problem class), and takes two optional arguments that are disabled by default. Setting ``suppress_warnings`` to ``true`` will suppress all JuMP-specific output (e.g. warnings about infeasibility and lack of dual information) but will not suppress solver output (which should be done by passing options to the solver). Setting ``relaxation=true`` solves the standard continuous relaxation for the model: that is, integrality is dropped, special ordered set constraints are not enforced, and semi-continuous and semi-integer variables with bounds ``[l,u]`` are replaced with bounds ``[min(l,0),max(u,0)]``.
* ``JuMP.build(m::Model)`` - builds the model in memory at the MathProgBase level without optimizing.
* ``setsolver(m::Model,s::AbstractMathProgSolver)`` - changes the solver which will be used for the next call to ``solve()``, discarding the current internal model if present.
* ``getvariable(m::Model,name::Symbol)`` - returns the variable or group of variables of the given name which were added to the model with ``@variable``. Throws an error if multiple variables were created with the same name.

**Objective**

* ``getobjective(m::Model)`` - returns the objective function as a ``QuadExpr``.
* ``getobjectivesense(m::Model)`` - returns objective sense, either ``:Min`` or ``:Max``.
* ``setobjectivesense(m::Model, newSense::Symbol)`` - sets the objective sense (``newSense`` is either ``:Min`` or ``:Max``).
* ``getobjectivevalue(m::Model)`` - returns objective value after a call to ``solve``.

**Output**

* ``writeLP(m::Model, filename::AbstractString)`` - write the model to ``filename`` in the LP file format.
* ``writeMPS(m::Model, filename::AbstractString)`` - write the model to ``filename`` in the MPS file format.

.. _solvestatus:

Solve status
^^^^^^^^^^^^

The call ``status = solve(m)`` returns a symbol recording the status of the optimization process, as reported by the solver. Typical values are listed in the table below, although the code can take solver-dependent values. For instance, certain solvers prove infeasibility or unboundedness during presolve, but do not report which of the two cases holds. See your solver interface documentation (as linked to in the :ref:`solver table <jump-solvertable>`) for more information.

.. _jump-statustable:

+-----------------+-----------------------------------------+
| Status          | Meaning                                 |
+=================+=========================================+
| ``:Optimal``    | Problem solved to optimality            |
+-----------------+-----------------------------------------+
| ``:Unbounded``  | Problem is unbounded                    |
+-----------------+-----------------------------------------+
| ``:Infeasible`` | Problem is infeasible                   |
+-----------------+-----------------------------------------+
| ``:UserLimit``  | Iteration limit or timeout              |
+-----------------+-----------------------------------------+
| ``:Error``      | Solver exited with an error             |
+-----------------+-----------------------------------------+
| ``:NotSolved``  | Model built in memory but not optimized |
+-----------------+-----------------------------------------+


Quadratic Objectives
^^^^^^^^^^^^^^^^^^^^

Quadratic objectives are supported by JuMP using a solver which implements the
corresponding extensions of the MathProgBase interface. Add them in the same way
you would a linear objective::

    m = Model()
    @variable(m, 0 <= x <= 2 )
    @variable(m, 0 <= y <= 30 )

    @objective(m, Min, x*x+ 2x*y + y*y )
    @constraint(m, x + y >= 1 )

    print(m)

    status = solve(m)

Second-order cone constraints
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Second-order cone constraints of the form :math:`||Ax-b||_2 + a^Tx + c \le 0` can be added directly using the ``norm`` function::

    @constraint(m, norm(A*x) <= 2w - 1)

The special ``norm2{...}`` construct may be used to build up normed expressions with complex indexing operations in much the same way as the ``sum{...}`` construct::

    @constraint(m, norm2{2x[i] - i, i=1:n; c[i] == 1} <= 1)

Accessing the low-level model
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

It is possible to construct the internal low-level model before optimizing. To do this,
call the ``JuMP.build`` function. It is then possible
to obtain this model by using the ``internalmodel`` function. This may be useful when
it is necessary to access some functionality that is not exposed by JuMP. When you are ready to optimize,
simply call ``solve`` in the normal fashion.

.. _probmod:

--------------------
Problem Modification
--------------------

It can be useful to modify models after they have been created and solved, for
example when we are solving many similar models in succession or generating the
model dynamically (e.g. column generation). Additionally it is sometimes
desirable for the solver to re-start from the last solution to reduce running
times for successive solves ("hot-start"). Where available, JuMP exposes this
functionality.

Differences in Solvers
^^^^^^^^^^^^^^^^^^^^^^

Some solvers do not expose the ability to modify a model after creation - the
model must be constructed from scratch each time. JuMP will use the ability to
modify problems exposed by the solver if possible, and will still work even if
the solver does not support this functionality by passing the complete problem
to the solver every time.

Modifying variables
^^^^^^^^^^^^^^^^^^^

As before, variables can be added using the ``@defVar`` macro. To remove a variable,
one can set the bounds on that variable to zero, e.g.::

    setLower(x, 0.0)
    setUpper(x, 0.0)

While bound updates are applied immediately in JuMP, variable bound changes are not
transmitted to the solver until ``solve`` is called again.

To add variables that appear in existing constraints, e.g. in column generation,
there is an alternative form of the ``defVar`` macro::

  @defVar(m, x, objective = objcoef, inconstraints = constrrefs, coefficients = values)
  @defVar(m, x >= lb, objective = objcoef, inconstraints = constrrefs, coefficients = values)
  @defVar(m, x <= ub, objective = objcoef, inconstraints = constrrefs, coefficients = values)
  @defVar(m, lb <= x <= ub, objective = objcoef, inconstraints = constrrefs, coefficients = values)
  @defVar(m, lb <= x <= ub, Int, objective = objcoef, inconstraints = constrrefs, coefficients = values)  # Types are supported

where ``objcoef`` is the coefficient of the variable in the new problem,
``constrrefs`` is a vector of ``ConstraintRef``, and ``values`` is a vector
of numbers. To give an example, consider the following code snippet::

  m = Model()
  @defVar(m, 0 <= x <= 1)
  @defVar(m, 0 <= y <= 1)
  @setObjective(m, Max, 5x + 1y)
  @addConstraint(m, con, x + y <= 1)
  solve(m)  # x = 1, y = 0
  @defVar(m, 0 <= z <= 1, objective = 10.0, inconstraints = [con], coefficients = [1.0])
  # The constraint is now x + y + z <= 1
  # The objective is now 5x + 1y + 10z
  solve(m)  # z = 1

In some situations you may be adding all variables in this way. To do so, first
define a set of empty constraints, e.g. ::

  m = Model()
  @addConstraint(m, con, 0 <= 1)
  @setObjective(m, Max, 0)
  @defVar(m, 0 <= x <= 1, objective = 5, inconstraints = [con], coefficients = [1.0])
  @defVar(m, 0 <= y <= 1, objective = 1, inconstraints = [con], coefficients = [1.0])
  @defVar(m, 0 <= z <= 1, objective = 10, inconstraints = [con], coefficients = [1.0])
  solve(m)

Modifying constraints
^^^^^^^^^^^^^^^^^^^^^

JuMP does not currently support changing constraint coefficients. For less-than
and greater-than constraints, the right-hand-side can be changed, e.g.::

    @addConstraint(m, mycon, x + y <= 4)
    solve(m)
    chgConstrRHS(mycon, 3)  # Now x + y <= 3
    solve(m)  # Hot-start for LPs

Modifying the objective
^^^^^^^^^^^^^^^^^^^^^^^

To change the objective, simply call ``@setObjective`` again - the previous objective
function and sense will be replaced.

Modifying nonlinear models
^^^^^^^^^^^^^^^^^^^^^^^^^^

See :ref:`nonlinear parameters <nonlinearprobmod>`.


.. _ref-expr:

---------------------------
Expressions and Constraints
---------------------------

Constructor
^^^^^^^^^^^

``AffExpr`` is an affine expression type defined by JuMP. It has three fields:
a vector of coefficients, a vector of variables, and a constant. Apart from
a default constructor that takes no arguments, it also has a full constructor that
can be useful if you want to manually build an affine expression::

    aff = AffExpr([x, z], [3.0, 4.0], 2.0)  # 3x + 4z + 2

Note that the coefficients must be floating point numbers. The matching
constraint for ``AffExpr`` is ``LinearConstraint`` which is defined by an
``AffExpr`` and a lower and upper bound. If a solver interface does not
support range constraints, this will automatically translated into two
constraints at solve time. Constructing constraints manually is not an
expected behavior and won't add the constraint to a model automatically.
See below for the correct methods.


There is also ``QuadExpr`` for quadratic expressions type that also provides
a default constructor that takes no arguments and a full constructor. There
are four fields: two vectors of variables, a vector of coefficients, and the
affine part of the expression. This is best explained by example::

    aff = AffExpr([x, z], [3.0, 4.0], 2.0)  # 3x + 4z + 2
    quad = QuadExpr([x,y],[x,z],[3.0,4.0],aff)  # 3x^2 + 4yz + 3x + 4z + 2

The corresponding constraint is ``QuadConstraint``, which is expected to
be a convex quadratic constraint.

Methods
^^^^^^^

* ``@addConstraint(m::Model, con)`` - efficient way to add linear or quadratic constraints.
* ``@addConstraint(m::Model, ref, con)`` - efficient way to add groups of linear or quadratic constraints.
  See Constraint Reference section for details.
* ``addConstraint(m::Model, con)`` - general way to add linear and quadratic
  constraints.
* ``@addConstraints`` - add groups of constraints at once, in the same fashion as @addConstraint. The model must be the first argument, and multiple constraints can be added on multiple lines wrapped in a ``begin ... end`` block. For example::

    @addConstraints(m, begin
      x >= 1
      y - w <= 2
      sum_to_one[i=1:3], z[i] + y == 1
    end)

* ``@defExpr(ref, expr)`` - efficiently builds a linear or quadratic expression but does not add to model immediately. Instead, returns the expression which can then be inserted in other constraints. For example::

    @defExpr(shared, sum{i*x[i], i=1:5})
    @addConstraint(m, shared + y >= 5)
    @addConstraint(m, shared + z <= 10)

The ``ref`` accepts index sets in the same way as ``@defVar``, and those indices can be used in the construction of the expressions::

    @defExpr(expr[i=1:3], i*sum{x[j], j=1:3})

* ``addSOS1(m::Model, coll::Vector{AffExpr})`` - adds special ordered set constraint
  of type 1 (SOS1). Specify the set as a vector of weighted variables, e.g. ``coll = [3x, y, 2z]``.
  Note that solvers expect the weights to be unique. See
  `here <http://lpsolve.sourceforge.net/5.5/SOS.htm>`_ for more details. If there is no inherent
  weighting in your model, an SOS constraint is probably unnecessary.
* ``addSOS2(m::Model, coll::Vector{AffExpr})`` - adds special ordered set constraint
  of type 2 (SOS2). Specify the set as a vector of weighted variables, e.g. ``coll = [3x, y, 2z]``.
  Note that solvers expect the weights to be unique.
  See `here <http://lpsolve.sourceforge.net/5.5/SOS.htm>`_ for more details.
* ``push!(aff::AffExpr, new_coeff::Float64, new_var::Variable)`` - efficient
  way to grow an affine expression by one term. For example, to add ``5x`` to
  an existing expression ``aff``, use ``push!(aff, 5.0, x)``. This is
  significantly more efficient than ``aff += 5.0*x``.
* ``append!(aff::AffExpr, other::AffExpr)`` - efficiently append the terms of
  an affine expression to an existing affine expression. For example, given
  ``aff = 5.0*x`` and ``other = 7.0*y + 3.0*z``, we can grow ``aff`` using
  ``append!(aff, other)`` which results in ``aff`` equaling ``5x + 7y + 3z``.
  This is significantly more efficient than using ``aff += other``.
* ``sum(affs::Array{AffExpr})`` - efficiently sum an array of affine expressions.
* ``getValue(expr)`` - evaluate an ``AffExpr`` or ``QuadExpr``, given the current solution values.

Constraint References
^^^^^^^^^^^^^^^^^^^^^

In order to manipulate constraints after creation, it is necessary to maintain
a reference. The simplest way to do this is to use the special three-argument
named constraint syntax for ``@addConstraint``, which additionally allows you
to create groups of constraints indexed by sets analogously to ``@defVar``.
For example::

    @defVar(m, x[1:3])
    @defVar(m, y[2:2:6])
    @addConstraint(m, xyconstr[i=1:3,j=6:-2:2], x[i] - y[j] == 1)

adds 9 constraints to the model ``m`` of the expected form. The variable ``xyconstr``
is a collection of ``ConstraintRef{LinearConstraint}`` instances indexed
by the ranges ``1:3`` and ``6:-2:2`` (the ordered tuple ``(6,4,2)``), so, for example
``xyconstr[2,4]`` is a reference to the constraint ``x[2] - y[4] == 1``. Indices can
have dependencies on preceding indices, e.g. triangular indexing is allowed::

    @addConstraint(m, triconstr[i=1:3,j=2i:2:6], x[i] - y[j] == 1)

To obtain the dual of a constraint, call ``getDual`` on the constraint reference::

    println(getDual(xyconstr[1,6]))

When an LP model is infeasible, ``getDual`` will return the corresponding component of the
infeasibility ray (Farkas proof), if available from the solver.

Dual information is unavailable for MIPs and has not yet been implemented for quadratic constraints.

For users who prefer to generate constraints in an explicit loop, we also
provide the ``@defConstrRef`` convenience macro, e.g.::

    @defConstrRef constraintName[1:3]

You can then iterate over constraints and store
references in this structure, e.g.::

    @defVar(m, x[1:5] >= 0)
    @defConstrRef myCons[1:5]
    for i = 1:5
      myCons[i] = @addConstraint(m, x[i] >= i)
    end


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

    aff = AffExpr([3.0, 4.0], [x, z], 2.0)  # 3x + 4z + 2

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

    aff = AffExpr([3.0, 4.0], [x, z], 2.0)  # 3x + 4z + 2
    quad = QuadExpr([x,y],[x,z],[3.0,4.0],aff)  # 3x^2 + 4yz + 3x + 4z + 2

The corresponding constraint is ``QuadConstraint``, which is expected to
be a convex quadratic constraint.

Methods
^^^^^^^

* ``@addConstraint(m::Model, con)`` - efficient way to add linear constraints.
  Uses macros and thus does not yet support quadratic constraints.
* ``addConstraint(m::Model, con)`` - general way to add linear and quadratic
  constraints.
* ``push!(aff::AffExpr, new_coeff::Float64, new_var::Variable)`` - efficient
  way to grow an affine expression by one term. For example, to add ``5x`` to
  an existing expression ``aff``, use ``push!(aff, 5.0, x)``. This is
  signficantly more efficient than ``aff += 5.0*x``.
* ``append!(aff::AffExpr, other::AffExpr)`` - efficiently append the terms of
  an affine expression to an existing affine expression. For example, given
  ``aff = 5.0*x`` and ``other = 7.0*y + 3.0*z``, we can grow ``aff`` using
  ``append!(aff, other)`` which results in ``aff`` equaling ``5x + 7y + 3z``.
  This is signficantly more efficient than using ``aff += other``.

Constraint References
^^^^^^^^^^^^^^^^^^^^^

In order to manipulate constraints after creation, it is necessary to maintain
a reference. For linear constraints both ``@addConstraint`` and ``addConstraint``
return an object of type ``ConstraintRef{LinearConstraint}``. To facilitate
the storage of these we provide the convenience macro, e.g.::

    @defConstrRef constraintName[1:3]

That behaves like ``@defVar``. You can then iterate over constraints and store
references in this structure, e.g.::

    @defVar(m, x[1:5] >= 0)
    @defConstrRef myCons[1:5]
    for i = 1:5
      myCons[i] = @addConstraint(m, x[i] >= i)
    end

To obtain the dual of a constraint, call ``getDual`` on the constraint reference::
    
    println(getDual(myCons[1]))

Dual information is unavailable for MIPs and has not yet been implemented for quadratic constraints.

When an LP model is infeasible, ``getDual`` will return the corresponding component of the
infeasibility ray (Farkas proof), if available from the solver.

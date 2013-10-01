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
expected behaviour and won't add the constraint to a model automatically.
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

Both ``@addConstraint`` and ``addConstraint`` return a.


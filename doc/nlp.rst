.. _nonlinear:

------------------
Nonlinear modeling
------------------

JuMP has *experimental* support for general smooth nonlinear (convex and 
nonconvex) optimization problems. JuMP is able to provide exact second-order 
derivatives to solvers. This information can improve solver accuracy and
performance.

Currently, `Ipopt <https://projects.coin-or.org/Ipopt>`_
is the only supported solver. To install Ipopt, run::

    Pkg.add("Ipopt")

Nonlinear objectives and constraints are specified by using the ``@setNLObjective``
and ``@addNLConstraint`` macros. The familiar ``sum{}`` syntax is supported within
these macros, as well as ``prod{}`` which analogously represents the product of
the terms within. Note that the ``@setObjective`` and ``@addConstraint``
macros (and corresponding functions) do *not* currently support nonlinear expressions.
However, a model can contain a mix of linear, quadratic, and nonlinear constraints or
objective functions.  Starting points must be provided by calling ``setValue`` on each
variable. For nonconvex problems, the returned solution is only guaranteed to be
locally optimal. Convexity detection is not currently provided.

For example, we can solve the classical Rosenbrock problem (with a twist) as follows::

    using JuMP
    m = Model()
    @defVar(m, x)
    @defVar(m, y)

    setValue(x, 0.0); setValue(y, 0.0)
    @setNLObjective(m, Min, (1-x)^2 + 100(y-x^2)^2)

    solve(m)
    println("x = ", getValue(x), " y = ", getValue(y))

    # adding a (linear) constraint
    @addConstraint(m, x + y == 10)
    solve(m)
    println("x = ", getValue(x), " y = ", getValue(y))


 

.. _nonlinear:

------------------
Nonlinear Modeling
------------------

JuMP has support for general smooth nonlinear (convex and
nonconvex) optimization problems. JuMP is able to provide exact, sparse second-order
derivatives to solvers. This information can improve solver accuracy and
performance.




Nonlinear objectives and constraints are specified by using the ``@setNLObjective``
and ``@addNLConstraint`` macros. The familiar ``sum{}`` syntax is supported within
these macros, as well as ``prod{}`` which analogously represents the product of
the terms within. Note that the ``@setObjective`` and ``@addConstraint``
macros (and corresponding functions) do *not* currently support nonlinear expressions.
However, a model can contain a mix of linear, quadratic, and nonlinear constraints or
objective functions.  Starting points may be provided by calling ``setValue`` on each
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

Examples: `optimal control <https://github.com/JuliaOpt/JuMP.jl/blob/master/examples/optcontrol.jl>`_, `maximum likelihood estimation <https://github.com/JuliaOpt/JuMP.jl/blob/master/examples/mle.jl>`_, and  `Hock-Schittkowski tests <https://github.com/JuliaOpt/JuMP.jl/tree/master/test/hockschittkowski>`_.

Performance
^^^^^^^^^^^

The execution time when solving a nonlinear programming problem can be divided into two parts, the time spent in the optimization algorithm (the solver) and the time spent evaluating the nonlinear functions and corresponding derivatives. Ipopt explicitly displays these two timings in its output, for example:

.. code-block:: text

    Total CPU secs in IPOPT (w/o function evaluations)   =      7.412
    Total CPU secs in NLP function evaluations           =      2.083
    

For Ipopt in particular, one can improve the performance by installing advanced sparse linear algebra packages, see :ref:`jump-installation`. For other solvers, see their respective documentation for performance tips.

The function evaluation time, on the other hand, is the responsibility of the modeling language. JuMP computes derivatives by using the `ReverseDiffSparse <https://github.com/mlubin/ReverseDiffSparse.jl>`_ package, which implements, in pure Julia, reverse-mode automatic differentiation with state-of-the-art graph coloring methods [1]_ for exploiting sparsity of the Hessian matrix. As a conservative bound, JuMP's performance here currently may be expected to be within a factor of 10 of AMPL's.

.. note::

    JuMP's performance for evaluating derivatives significantly improves if the `ArrayViews <https://github.com/lindahua/ArrayViews.jl>`_ package is installed. This is currently an optional dependency.


.. [1] Gebremdhin et al., "Efficient Computation of Sparse Hessians Using Coloring and Automatic Differentiation", INFORMS Journal on Computing, 21(1), pp. 209-223, 2009.

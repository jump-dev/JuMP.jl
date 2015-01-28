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
objective functions.  Starting points may be provided by using the ``start``
keyword argument to ``@defVar``.
For nonconvex problems, the returned solution is only guaranteed to be
locally optimal. Convexity detection is not currently provided.

For example, we can solve the classical Rosenbrock problem (with a twist) as follows::

    using JuMP
    m = Model()
    @defVar(m, x, start = 0.0)
    @defVar(m, y, start = 0.0)

    @setNLObjective(m, Min, (1-x)^2 + 100(y-x^2)^2)

    solve(m)
    println("x = ", getValue(x), " y = ", getValue(y))

    # adding a (linear) constraint
    @addConstraint(m, x + y == 10)
    solve(m)
    println("x = ", getValue(x), " y = ", getValue(y))

Examples: `optimal control <https://github.com/JuliaOpt/JuMP.jl/blob/master/examples/optcontrol.jl>`_, `maximum likelihood estimation <https://github.com/JuliaOpt/JuMP.jl/blob/master/examples/mle.jl>`_, and  `Hock-Schittkowski tests <https://github.com/JuliaOpt/JuMP.jl/tree/master/test/hockschittkowski>`_.

Performance: Solution time
^^^^^^^^^^^^^^^^^^^^^^^^^^

The execution time when *solving* a nonlinear programming problem can be divided into two parts, the time spent in the optimization algorithm (the solver) and the time spent evaluating the nonlinear functions and corresponding derivatives. Ipopt explicitly displays these two timings in its output, for example:

.. code-block:: text

    Total CPU secs in IPOPT (w/o function evaluations)   =      7.412
    Total CPU secs in NLP function evaluations           =      2.083


For Ipopt in particular, one can improve the performance by installing advanced sparse linear algebra packages, see :ref:`jump-installation`. For other solvers, see their respective documentation for performance tips.

The function evaluation time, on the other hand, is the responsibility of the modeling language. JuMP computes derivatives by using the `ReverseDiffSparse <https://github.com/mlubin/ReverseDiffSparse.jl>`_ package, which implements, in pure Julia, reverse-mode automatic differentiation with state-of-the-art graph coloring methods [1]_ for exploiting sparsity of the Hessian matrix. As a conservative bound, JuMP's performance here currently may be expected to be within a factor of 5 of AMPL's.

.. _nonlinearprobmod:

Performance: Model generation
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^


Before the solution process starts, JuMP processes the model in order to set up callbacks which solvers use to query derivative information. Typically this "model generation" time is not a bottleneck; however, it could be relatively large when solving a sequence of small problems. **It is possible to skip model generation** after the first call to ``solve()`` if the only differences between models are changes in input data; if constraints are added, the model will be rebuilt from scratch.

Unlike some other modeling languages, JuMP does not currently have a special syntax for "parameters", that is, symbolic constants whose values can be changed between solves. Instead, JuMP uses Julia's rules for `variable bindings <http://docs.julialang.org/en/release-0.3/manual/faq/#i-passed-an-argument-x-to-a-function-modified-it-inside-that-function-but-on-the-outside-the-variable-x-is-still-unchanged-why>`_. Internally, when a variable appears in ``@addNLConstraint`` or ``@setNLObjective``, JuMP saves the object that the variable is bound to. If this object is later modified, the corresponding changes will be reflected in JuMP's function evaluations. For example::

    using JuMP
    m = Model()
    @defVar(m, 0.5 <= x <=  2)
    @defVar(m, 0.0 <= y <= 30)
    @setObjective(m, Min, (x+y)^2)
    param = [1.0]
    @addNLConstraint(m, x + y >= param[1])
    solve(m)
    # Optimal objective is 1.0

    # modify the value saved by JuMP
    param[1] = 10.0
    solve(m)
    # optimal objective is 10.0^2

Note that we used a vector ``param``, which is a mutable object. On the other hand, the following code would *not* result in any modifications to the JuMP model::

    param = 1.0
    @addNLConstraint(m, x + y >= param)
    param = 10.0

The line ``param = 10.0`` changes ``param`` to reference a new value in the local scope, but does not affect the value referenced by JuMP.

This variable binding trick for quick model regeneration does not apply to the macros ``@addConstraint`` and ``@setObjective`` for linear and quadratic expressions; see :ref:`probmod` for modifying linear models. We hope to treat in-place model modifications in a more uniform manner in future releases.

.. [1] Gebremdhin et al., "Efficient Computation of Sparse Hessians Using Coloring and Automatic Differentiation", INFORMS Journal on Computing, 21(1), pp. 209-223, 2009.

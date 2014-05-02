.. _nonlinear:

------------------
Nonlinear Modeling
------------------

.. note::
    Nonlinear modeling functionality should be considered as a feature preview. It has not been as widely tested as the (mixed-integer) linear, quadratic, and conic capabilities of JuMP.

JuMP has support for general smooth nonlinear (convex and
nonconvex) optimization problems. JuMP is able to provide exact, sparse second-order
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


Solver Options
^^^^^^^^^^^^^^

Solution options for Ipopt can be passed by using the optional ``IpoptOptions`` parameter to ``solve()``::

    ...
    solve(m, IpoptOptions=[("tol",1e-6)])
    ...

to set the relative convergence tolerance to ``1e-6``. A full list of options is available `here <http://www.coin-or.org/Ipopt/documentation/node41.html>`_. This syntax is temporary and will be replaced in the future once nonlinear modeling extensions are implemented at the solver-independent MathProgBase level.

Performance
^^^^^^^^^^^

The execution time when solving a nonlinear programming problem can be divided into two parts, the time spent in the optimization algorithm (the solver) and the time spent evaluating the nonlinear functions and corresponding derivatives. Ipopt explicitly displays these two timings in its output, for example:

.. code-block:: text

    Total CPU secs in IPOPT (w/o function evaluations)   =      7.412
    Total CPU secs in NLP function evaluations           =      2.083
    

The performance of the solver itself greatly depends on the linear algebra libraries used. By default, the Ipopt binaries installed by the Ipopt.jl package use the open-source MUMPS library for sparse linear algebra. Significant speedups can be obtained by manually compiling Ipopt to use proprietary sparse linear algebra libraries instead. Julia can be pointed to use a custom version of Ipopt; we suggest posting to the `julia-opt <https://groups.google.com/forum/#!forum/julia-opt>`_ mailing list with your platform details for guidance on how to do this.

The function evaluation time is the responsibility of the modeling language. JuMP computes derivatives by using the `ReverseDiffSparse <https://github.com/mlubin/ReverseDiffSparse.jl>`_ package, which implements, in pure Julia, reverse-mode automatic differentiation with state-of-the-art graph coloring methods [1]_ for exploiting sparsity of the Hessian matrix. As a conservative bound, JuMP's performance here currently may be expected to be within a factor of 10 of AMPL's.

.. note::

    JuMP's performance for evaluating derivatives significantly improves if the `ArrayViews <https://github.com/lindahua/ArrayViews.jl>`_ package is installed. This is currently an optional dependency because it requires the 0.3 development version of Julia.


.. [1] Gebremdhin et al., "Efficient Computation of Sparse Hessians Using Coloring and Automatic Differentiation", INFORMS Journal on Computing, 21(1), pp. 209-223, 2009.

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
If a starting value is not provided for a variable, it will be set to the projection
of zero onto the interval defined by the variable bounds.
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

Syntax notes
^^^^^^^^^^^^

The syntax accepted in nonlinear expressions is more restricted than
the syntax for linear and quadratic expressions. We note some important points below.

- All expressions must be simple scalar operations. You cannot use ``dot``,
  matrix-vector products, vector slices, etc. Translate vector operations
  into explicit ``sum{}`` operations or use the ``AffExpr`` plus auxiliary variable
  trick described below.
- There is no operator overloading provided to build up nonlinear expressions.
  For example, if ``x`` is a JuMP variable, the code ``3x`` will return an
  ``AffExpr`` object that can be used inside of future expressions and
  linear constraints.
  However, the code ``sin(x)`` is an error. All nonlinear expressions must
  be inside of macros.
- As a corollary, user-defined functions may not be used within nonlinear
  expressions. See the example below::

    myfunction(a,b) = exp(a)*b
    @defVar(m, x); @defVar(m, y)
    @setNLObjective(m, Min, myfunction(x,y)) # ERROR
    @setNLObjective(m, Min, exp(x)*y) # Okay

- ``AffExpr`` and ``QuadExpr`` objects cannot currently be used inside nonlinear
  expressions. Instead, introduce auxiliary variables, e.g.::

    myexpr = dot(c,x) + 3y # where x and y are variables
    @defVar(m, aux)
    @addConstraint(m, aux == myexpr)
    @setNLObjective(m, Min, sin(aux))
- You can declare embeddable nonlinear expressions with ``@defNLExpr``. For example::

    @defNLExpr(myexpr[i=1:n], sin(x[i]))
    @addNLConstraint(m, myconstr[i=1:n], myexpr[i] <= 0.5)

.. note::
    There is currently no validity check on indices used for parametric nonlinear expressions. Therefore the following example will not raise an error::

        @defNLExpr(myexpr[i=1:2], i * sin(x))
        @addNLConstraint(m, myexpr[0] <= 0.5)

    You will still receive an error in the case that invalid indices are used to access data or variables within nonlinear expressions.
    We do not recommend depending on this behavior, since it may change in a future JuMP release.

Performance: Solution time
^^^^^^^^^^^^^^^^^^^^^^^^^^

The execution time when *solving* a nonlinear programming problem can be divided into two parts, the time spent in the optimization algorithm (the solver) and the time spent evaluating the nonlinear functions and corresponding derivatives. Ipopt explicitly displays these two timings in its output, for example:

.. code-block:: text

    Total CPU secs in IPOPT (w/o function evaluations)   =      7.412
    Total CPU secs in NLP function evaluations           =      2.083


For Ipopt in particular, one can improve the performance by installing advanced sparse linear algebra packages, see :ref:`jump-installation`. For other solvers, see their respective documentation for performance tips.

The function evaluation time, on the other hand, is the responsibility of the modeling language. JuMP computes derivatives by using the `ReverseDiffSparse <https://github.com/mlubin/ReverseDiffSparse.jl>`_ package, which implements, in pure Julia, reverse-mode automatic differentiation with graph coloring methods for exploiting sparsity of the Hessian matrix [1]_. As a conservative bound, JuMP's performance here currently may be expected to be within a factor of 5 of AMPL's.

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

Querying derivatives from a JuMP model
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

For some advanced use cases, one may want to directly query the derivatives
of a JuMP model instead of handing the problem off to a solver.
Internally, JuMP implements the ``AbstractNLPEvaluator`` interface from
`MathProgBase <http://mathprogbasejl.readthedocs.org/en/latest/nlp.html>`_.
To obtain an NLP evaluator object from a JuMP model, use ``JuMPNLPEvaluator``.
The ``getLinearIndex`` method maps from JuMP variables to the variable
indices at the MathProgBase level.

For example::

    m = Model()
    @defVar(m, x)
    @defVar(m, y)

    @setNLObjective(m, Min, sin(x) + sin(y))
    values = zeros(2)
    values[getLinearIndex(x)] = 2.0
    values[getLinearIndex(y)] = 3.0

    d = JuMPNLPEvaluator(m)
    MathProgBase.initialize(d, [:Grad])
    objval = MathProgBase.eval_f(d, values) # == sin(2.0) + sin(3.0)

    ∇f = zeros(2)
    MathProgBase.eval_grad_f(d, ∇f, values)
    # ∇f[getLinearIndex(x)] == cos(2.0)
    # ∇f[getLinearIndex(y)] == cos(3.0)

The ordering of constraints in a JuMP model corresponds to the following ordering
at the MathProgBase nonlinear abstraction layer. There are three groups of constraints:
linear, quadratic, and nonlinear. Linear and quadratic constraints, to be recognized
as such, must be added with the ``@addConstraint`` macros. All constraints added with
the ``@addNLConstraint`` macros are treated as nonlinear constraints.
Linear constraints are ordered first, then quadratic, then nonlinear.
The ``getLinearIndex`` method applied to a constraint reference object
returns the index of the constraint *within its corresponding constraint class*.
For example::

    m = Model()
    @defVar(m, x)
    @addConstraint(m, cons1, x^2 <= 1)
    @addConstraint(m, cons2, x + 1 == 3)
    @addNLConstraint(m, cons3, x + 5 == 10)

    typeof(cons1) # ConstraintRef{GenericQuadConstraint{GenericQuadExpr{Float64,Variable}}} indicates a quadratic constraint
    typeof(cons2) # ConstraintRef{GenericRangeConstraint{GenericAffExpr{Float64,Variable}}} indicates a linear constraint
    typeof(cons3) # ConstraintRef{GenericRangeConstraint{SymbolicOutput}} indicates a nonlinear constraint
    getLinearIndex(cons1) == getLinearIndex(cons2) == getLinearIndex(cons3) == 1

When querying derivatives, ``cons2`` will appear first, because it is the first linear constraint, then ``cons1``, because it is the first quadratic constraint, then ``cons3``, because it is the first nonlinear constraint. Note that for one-sided nonlinear constraints, JuMP subtracts any values on the right-hand side when computing expression. In other words, one-sided linear constraints are always transformed to have a right-hand side of zero.

This method of querying derivatives directly from a JuMP model is convenient for
interacting with the model in a structured way, e.g., for accessing derivatives of
specific variables. For example, in statistical maximum likelihood estimation problems,
one is often interested in the Hessian matrix at the optimal solution,
which can be queried using the ``JuMPNLPEvaluator``.

However, the examples above are *not* a convenient way to access the NLP `standard-form <http://mathprogbasejl.readthedocs.org/en/latest/nlp.html>`_ representation of a JuMP model, because there is no direct way to access the vector of constraint upper and lower bounds. In this case, you should implement an ``AbstractMathProgSolver`` and corresponding ``AbstractMathProgModel`` type following the MathProgBase nonlinear interface, collecting the problem data through the ``MathProgBase.loadnonlinearproblem!`` `method <http://mathprogbasejl.readthedocs.org/en/latest/nlp.html#loadnonlinearproblem!>`_. This approach has the advantage of being nominally independent of JuMP itself in terms of problem format. You may use the ``buildInternalModel`` method to ask JuMP to populate the "solver" without calling ``optimize!``.

.. [1] Dunning, Huchette, and Lubin, "JuMP: A Modeling Language for Mathematical Optimization", `arXiv <http://arxiv.org/abs/1508.01982>`_.

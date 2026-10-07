.. include:: warn.rst

.. _nonlinear:

------------------
Nonlinear Modeling
------------------

JuMP has support for general smooth nonlinear (convex and
nonconvex) optimization problems. JuMP is able to provide exact, sparse second-order
derivatives to solvers. This information can improve solver accuracy and
performance.




Nonlinear objectives and constraints are specified by using the ``@NLobjective``
and ``@NLconstraint`` macros. The familiar ``sum()`` syntax is supported within
these macros, as well as ``prod()`` which analogously represents the product of
the terms within. Note that the ``@objective`` and ``@constraint``
macros (and corresponding functions) do *not* currently support nonlinear expressions.
However, a model can contain a mix of linear, quadratic, and nonlinear constraints or
objective functions.  Starting points may be provided by using the ``start``
keyword argument to ``@variable``.
If a starting value is not provided for a variable, it will be set to the projection
of zero onto the interval defined by the variable bounds.
For nonconvex problems, the returned solution is only guaranteed to be
locally optimal. Convexity detection is not currently provided.

For example, we can solve the classical Rosenbrock problem (with a twist) as follows::

    using JuMP
    m = Model()
    @variable(m, x, start = 0.0)
    @variable(m, y, start = 0.0)

    @NLobjective(m, Min, (1-x)^2 + 100(y-x^2)^2)

    solve(m)
    println("x = ", getvalue(x), " y = ", getvalue(y))

    # adding a (linear) constraint
    @constraint(m, x + y == 10)
    solve(m)
    println("x = ", getvalue(x), " y = ", getvalue(y))

Examples: `optimal control <https://github.com/JuliaOpt/JuMP.jl/blob/master/examples/optcontrol.jl>`_, `maximum likelihood estimation <https://github.com/JuliaOpt/JuMP.jl/blob/master/examples/mle.jl>`_, and  `Hock-Schittkowski tests <https://github.com/JuliaOpt/JuMP.jl/tree/master/test/hockschittkowski>`_.

Syntax notes
^^^^^^^^^^^^

The syntax accepted in nonlinear expressions is more restricted than
the syntax for linear and quadratic expressions. We note some important points below.

- All expressions must be simple scalar operations. You cannot use ``dot``,
  matrix-vector products, vector slices, etc. Translate vector operations
  into explicit ``sum()`` operations or use the ``AffExpr`` plus auxiliary variable
  trick described below.
- There is no operator overloading provided to build up nonlinear expressions.
  For example, if ``x`` is a JuMP variable, the code ``3x`` will return an
  ``AffExpr`` object that can be used inside of future expressions and
  linear constraints.
  However, the code ``sin(x)`` is an error. All nonlinear expressions must
  be inside of macros.
- :ref:`userfunctions` may be used within nonlinear
  expressions only after they are registered.
  For example::

    myfunction(a,b) = exp(a)*b
    @variable(m, x); @variable(m, y)
    @NLobjective(m, Min, myfunction(x,y)) # ERROR. Needs JuMP.register() first.
    @NLobjective(m, Min, exp(x)*y) # Okay

- ``AffExpr`` and ``QuadExpr`` objects cannot currently be used inside nonlinear
  expressions. Instead, introduce auxiliary variables, e.g.::

    myexpr = dot(c,x) + 3y # where x and y are variables
    @variable(m, aux)
    @constraint(m, aux == myexpr)
    @NLobjective(m, Min, sin(aux))
- You can declare embeddable nonlinear expressions with ``@NLexpression``. For example::

    @NLexpression(m, myexpr[i=1:n], sin(x[i]))
    @NLconstraint(m, myconstr[i=1:n], myexpr[i] <= 0.5)

- Anonymous syntax is supported in ``@NLexpression`` and ``@NLconstraint``::

    myexpr = @NLexpression(m, [i=1:n], sin(x[i]))
    myconstr = @NLconstraint(m, [i=1:n], myexpr[i] <= 0.5)

.. _nonlinearprobmod:

Nonlinear Parameters
^^^^^^^^^^^^^^^^^^^^

For nonlinear models only, JuMP offers a syntax for explicit "parameter" objects
which can be used to modify a model in-place just by updating the value of
the parameter.
Nonlinear parameters are declared by using the ``@NLparameter`` macro and may
be indexed by arbitrary sets analogously to JuMP variables and expressions.
The initial value of the parameter must be provided
on the right-hand side of the ``==`` sign as seen below::

    @NLparameter(m, x == 10)
    @NLparameter(m, y[i=1:10] == my_data[i]) # set of parameters indexed from 1 to 10

You may use ``getvalue`` and ``setvalue`` to query or update the value of a parameter::

    getvalue(x) # 10, from above
    setvalue(y[4], 54.3) # y[4] now holds the value 54.3

Nonlinear parameters can be used *within nonlinear expressions* only::

    @variable(m, z)
    @objective(m, Max, x*z)       # error: x is a nonlinear parameter
    @NLobjective(m, Max, x*z)     # ok
    @expression(m, my_expr, x*z^2)      # error: x is a nonlinear parameter
    @NLexpression(m, my_nl_expr, x*z^2) # ok

Nonlinear parameters are useful when solving nonlinear models in a sequence::

    m = Model()
    @variable(m, z)
    @NLparameter(m, x == 1.0)
    @NLobjective(m, Min, (z-x)^2)
    solve(m)
    getvalue(z) # equals 1.0

    # Now, update the value of x to solve a different problem
    setvalue(x, 5.0)
    solve(m)
    getvalue(z) # equals 5.0

Using nonlinear parameters can be faster than creating a new model from scratch
with updated data because JuMP is able to avoid repeating a number of steps
in processing the model before handing it off to the solver.

.. _userfunctions:

User-defined functions
^^^^^^^^^^^^^^^^^^^^^^

JuMP's library of recognized univariate functions is derived from the `Calculus.jl <https://github.com/johnmyleswhite/Calculus.jl>`_ package. If you encounter a standard special function not currently supported by JuMP, consider contributing to the `list of derivative rules <https://github.com/johnmyleswhite/Calculus.jl/blob/cb42f3699177449a42bdc3461c8aea8777aa8c39/src/differentiate.jl#L115>`_ there. In addition to this built-in list of functions, it is possible to register custom (*user-defined*) nonlinear functions to use within nonlinear expressions. JuMP does not support black-box optimization, so all user-defined functions must provide derivatives in some form. Fortunately, JuMP supports **automatic differentiation of user-defined functions**, a feature to our knowledge not available in any comparable modeling systems.

.. note::
    Automatic differentiation is *not* finite differencing. JuMP's automatically computed derivatives are not subject to approximation error.

JuMP uses `ForwardDiff.jl <https://github.com/JuliaDiff/ForwardDiff.jl>`_ to perform automatic differentiation; see the ForwardDiff.jl `documentation <http://www.juliadiff.org/ForwardDiff.jl/v0.5.0/user/limitations.html>`_ for a description of how to write a function suitable for automatic differentiation. The general guideline is to write code that is generic with respect to the number type; don't assume that the input to the function is ``Float64``. To register a user-defined function with derivatives computed by automatic differentiation, use the ``JuMP.register`` method as in the following example::

    mysquare(x) = x^2
    myf(x,y) = (x-1)^2+(y-2)^2

    m = Model()

    JuMP.register(m, :myf, 2, myf, autodiff=true)
    JuMP.register(m, :mysquare, 1, mysquare, autodiff=true)

    @variable(m, x[1:2] >= 0.5)
    @NLobjective(m, Min, myf(x[1],mysquare(x[2])))

The above code creates a JuMP model with the objective function ``(x[1]-1)^2 + (x[2]^2-2)^2``. The first argument to ``JuMP.register`` the model for which the functions are registered. The second argument is a Julia symbol object which serves as the name of the user-defined function in JuMP expressions; the JuMP name need not be the same as the name of the corresponding Julia method. The third argument specifies how many arguments the function takes. The fourth argument is the name of the Julia method which computes the function, and ``autodiff=true`` instructs JuMP to compute exact gradients automatically.

.. note::
    All arguments to user-defined functions are scalars, not vectors. To define a function which takes a large number of arguments, you may use the splatting syntax ``f(x...) = ...``.

Forward-mode automatic differentiation as implemented by ForwardDiff.jl has a computational cost that scales linearly with the number of input dimensions. As such, it is not the most efficient way to compute gradients of user-defined functions if the number of input arguments is large. In this case, users may want to provide their own routines for evaluating gradients. The more general syntax for ``JuMP.register`` which accepts user-provided derivative evaluation routines is::

    JuMP.register(m::Model, s::Symbol, dimension::Integer, f::Function, ∇f::Function, ∇²f::Function)

The input differs for functions which take a single input argument and functions which take more than one. For univariate functions, the derivative evaluation routines should return a number which represents the first and second-order derivatives respectively. For multivariate functions, the derivative evaluation routines will be passed a gradient vector which they must explicitly fill. Second-order derivatives of multivariate functions are not currently supported; this argument should be omitted. The following example sets up the same optimization problem as before, but now we explicitly provide evaluation routines for the user-defined functions::

    mysquare(x) = x^2
    mysquare_prime(x) = 2x
    mysquare_primeprime(x) = 2

    myf(x,y) = (x-1)^2+(y-2)^2
    function ∇f(g,x,y)
        g[1] = 2*(x-1)
        g[2] = 2*(y-2)
    end

    m = Model()

    JuMP.register(m, :myf, 2, myf, ∇f)
    JuMP.register(m, :mysquare, 1, mysquare, mysquare_prime, mysquare_primeprime)

    @variable(m, x[1:2] >= 0.5)
    @NLobjective(m, Min, myf(x[1],mysquare(x[2])))

Factors affecting solution time
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The execution time when solving a nonlinear programming problem can be divided into two parts, the time spent in the optimization algorithm (the solver) and the time spent evaluating the nonlinear functions and corresponding derivatives. Ipopt explicitly displays these two timings in its output, for example:

.. code-block:: text

    Total CPU secs in IPOPT (w/o function evaluations)   =      7.412
    Total CPU secs in NLP function evaluations           =      2.083


For Ipopt in particular, one can improve the performance by installing advanced sparse linear algebra packages, see :ref:`jump-installation`. For other solvers, see their respective documentation for performance tips.

The function evaluation time, on the other hand, is the responsibility of the modeling language. JuMP computes derivatives by using the `ReverseDiffSparse <https://github.com/mlubin/ReverseDiffSparse.jl>`_ package, which implements, in pure Julia, reverse-mode automatic differentiation with graph coloring methods for exploiting sparsity of the Hessian matrix [1]_. As a conservative bound, JuMP's performance here currently may be expected to be within a factor of 5 of AMPL's.



Querying derivatives from a JuMP model
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

For some advanced use cases, one may want to directly query the derivatives
of a JuMP model instead of handing the problem off to a solver.
Internally, JuMP implements the ``AbstractNLPEvaluator`` interface from
`MathProgBase <http://mathprogbasejl.readthedocs.org/en/latest/nlp.html>`_.
To obtain an NLP evaluator object from a JuMP model, use ``JuMP.NLPEvaluator``.
The ``linearindex`` method maps from JuMP variables to the variable
indices at the MathProgBase level.

For example::

    m = Model()
    @variable(m, x)
    @variable(m, y)

    @NLobjective(m, Min, sin(x) + sin(y))
    values = zeros(2)
    values[linearindex(x)] = 2.0
    values[linearindex(y)] = 3.0

    d = JuMP.NLPEvaluator(m)
    MathProgBase.initialize(d, [:Grad])
    objval = MathProgBase.eval_f(d, values) # == sin(2.0) + sin(3.0)

    ∇f = zeros(2)
    MathProgBase.eval_grad_f(d, ∇f, values)
    # ∇f[linearindex(x)] == cos(2.0)
    # ∇f[linearindex(y)] == cos(3.0)

The ordering of constraints in a JuMP model corresponds to the following ordering
at the MathProgBase nonlinear abstraction layer. There are three groups of constraints:
linear, quadratic, and nonlinear. Linear and quadratic constraints, to be recognized
as such, must be added with the ``@constraint`` macros. All constraints added with
the ``@NLconstraint`` macros are treated as nonlinear constraints.
Linear constraints are ordered first, then quadratic, then nonlinear.
The ``linearindex`` method applied to a constraint reference object
returns the index of the constraint *within its corresponding constraint class*.
For example::

    m = Model()
    @variable(m, x)
    @constraint(m, cons1, x^2 <= 1)
    @constraint(m, cons2, x + 1 == 3)
    @NLconstraint(m, cons3, x + 5 == 10)

    typeof(cons1) # JuMP.ConstraintRef{JuMP.Model,JuMP.GenericQuadConstraint{JuMP.GenericQuadExpr{Float64,JuMP.Variable}}} indicates a quadratic constraint
    typeof(cons2) # JuMP.ConstraintRef{JuMP.Model,JuMP.GenericRangeConstraint{JuMP.GenericAffExpr{Float64,JuMP.Variable}}} indicates a linear constraint
    typeof(cons3) # JuMP.ConstraintRef{JuMP.Model,JuMP.GenericRangeConstraint{JuMP.NonlinearExprData}} indicates a nonlinear constraint
    linearindex(cons1) == linearindex(cons2) == linearindex(cons3) == 1

When querying derivatives, ``cons2`` will appear first, because it is the first linear constraint, then ``cons1``, because it is the first quadratic constraint, then ``cons3``, because it is the first nonlinear constraint. Note that for one-sided nonlinear constraints, JuMP subtracts any values on the right-hand side when computing expressions. In other words, one-sided nonlinear constraints are always transformed to have a right-hand side of zero.

The ``JuMP.constraintbounds(m::Model)`` method returns the lower and upper bounds
of all the constraints in the model, concatenated in the order discussed above.

This method of querying derivatives directly from a JuMP model is convenient for
interacting with the model in a structured way, e.g., for accessing derivatives of
specific variables. For example, in statistical maximum likelihood estimation problems,
one is often interested in the Hessian matrix at the optimal solution,
which can be queried using the ``JuMP.NLPEvaluator``.

If you are writing a "solver", we *highly encourage* use of the `MathProgBase nonlinear interface <http://mathprogbasejl.readthedocs.org/en/latest/nlp.html>`_ over querying derivatives using the above methods. These methods are provided for convenience but do not fully integrate with JuMP's solver infrastructure. In particular, they do not allow users to specify your solver to the ``Model()`` constructor nor to call it using ``solve()`` nor to populate the solution back into the model. Use of the MathProgBase interface also has the advantage of being independent of JuMP itself; users of MathProgBase solvers are free to implement their own evaluation routines instead of expressing their model in JuMP.  You may use the ``JuMP.build`` method to ask JuMP to populate the "solver" without calling ``optimize!``.

Raw expression input
^^^^^^^^^^^^^^^^^^^^

In addition to the ``@NLobjective`` and ``@NLconstraint`` macros, it is also
possible to provide Julia ``Expr`` objects directly by using
``JuMP.setNLobjective`` and ``JuMP.addNLconstraint``. This input form
may be useful if the expressions are generated programmatically.
JuMP variables should be spliced into the expression object. For example::

    @variable(m, 1 <= x[i=1:4] <= 5)
    JuMP.setNLobjective(m, :Min, :($(x[1])*$(x[4])*($(x[1])+$(x[2])+$(x[3])) + $(x[3])))
    JuMP.addNLconstraint(m, :($(x[1])*$(x[2])*$(x[3])*$(x[4]) >= 25))

    # Equivalent form using traditional JuMP macros:
    @NLobjective(m, Min, x[1]*x[4]*(x[1]+x[2]+x[3]) + x[3])
    @NLconstraint(m, x[1]*x[2]*x[3]*x[4] >= 25)

See the Julia documentation for more examples and description of Julia expressions.


.. [1] Dunning, Huchette, and Lubin, "JuMP: A Modeling Language for Mathematical Optimization", `arXiv <http://arxiv.org/abs/1508.01982>`_.

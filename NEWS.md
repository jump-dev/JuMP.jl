JuMP release notes
==================

Version 0.13.2 (May 16, 2016)
-----------------------------

  * Compatibility update for MathProgBase

Version 0.13.1 (May 3, 2016)
----------------------------

  * Fix broken deprecation for ``registerNLfunction``.

Version 0.13.0 (April 29, 2016)
-------------------------------

  * Most exported methods and macros have been renamed to avoid camelCase. See the list of changes [here](https://github.com/JuliaOpt/JuMP.jl/blob/e53d0db67cde2a4b80d0c1281f4b49eb0128a1f5/src/deprecated.jl#L30). There is a 1-1 mapping from the old names to the new, and it is safe to simply replace the names to update existing models.
  * Specify variable lower/upper bounds in ``@variable`` using the ``lowerbound`` and ``upperbound`` keyword arguments.
  * Change name printed for variable using the ``basename`` keyword argument to ``@variable``.
  * New ``@variables`` macro allows multiline declaration of groups of variables.
  * A number of solver methods previously available only through MathProgBase are now exposed directly in JuMP. The fix was [recorded](https://youtu.be/qF1lZPJ3a5A) live!
  * Compatibility fixes with Julia 0.5.
  * The "end" indexing syntax is no longer supported within JuMPArrays which do not use 1-based indexing until upstream issues are resolved, see [here](https://github.com/JuliaOpt/JuMP.jl/issues/730).

Version 0.12.2 (March 9, 2016)
------------------------------

  * Small fixes for nonlinear optimization

Version 0.12.1 (March 1, 2016)
------------------------------

  * Fix a regression in slicing for JuMPArrays (when not using 1-based indexing)

Version 0.12.0 (February 27, 2016)
----------------------------------

  * The automatic differentiation functionality has been completely rewritten with a number of user-facing changes:
      - ``@defExpr`` and ``@defNLExpr`` now take the model as the first argument. The previous one-argument version of ``@defExpr`` is deprecated; all expressions should be named. E.g., replace ``@defExpr(2x+y)`` with ``@defExpr(jump_model, my_expr, 2x+y)``.
      - JuMP no longer uses Julia's variable binding rules for efficiently re-solving a sequence of nonlinear models. Instead, we have introduced nonlinear parameters. This is a breaking change, so we have added a warning message when we detect models that may depend on the old behavior.
      - Support for user-defined functions integrated within nonlinear JuMP expressions.
  * Replaced iteration over ``AffExpr`` with ``Number``-like scalar iteration; previous iteration behavior is now available via ``linearterms(::AffExpr)``.
  * Stopping the solver via ``throw(CallbackAbort())`` from a callback no longer triggers an exception. Instead, ``solve()`` returns ``UserLimit`` status.
  * ``getDual()`` now works for conic problems (Thanks @emreyamangil.)

Version 0.11.3 (February 4, 2016)
---------------------------------

  * Bug-fix for problems with quadratic objectives and semidefinite constraints

Version 0.11.2 (January 14, 2016)
---------------------------------

  * Compatibility update for Mosek

Version 0.11.1 (December 1, 2015)
---------------------------------

  * Remove usage of `@compat` in tests.
  * Fix updating quadratic objectives for nonlinear models.

Version 0.11.0 (November 30, 2015)
----------------------------------

  * Julia 0.4.0 is the minimum required version for this release.
  * Fix for scoping semantics of index variables in sum{}. Index variables no longer leak into the surrounding scope.
  * Addition of the ``solve(m::Model, relaxation=true)`` keyword argument to solve the standard continuous realaxation of model ``m``
  * The ``getConstraintBounds()`` method allows access to the lower and upper bounds of all constraints in a (nonlinear) model.
  * Update for breaking changes in MathProgBase

Version 0.10.3 (November 20, 2015)
----------------------------------

  * Fix a rare error when parsing quadratic expressions
  * Fix ``Variable()`` constructor with default arguments
  * Detect unrecognized keywords in ``solve()``

Version 0.10.2 (September 28, 2015)
-----------------------------------

  * Fix for deprecation warnings


Version 0.10.1 (September 3, 2015)
----------------------------------

  * Fixes for ambiguity warnings.
  * Fix for breaking change in precompilation syntax in Julia 0.4-pre

Version 0.10.0 (August 31, 2015)
--------------------------------

  * Support (on Julia 0.4 and later) for conditions in indexing ``@defVar`` and ``@addConstraint`` constructs, e.g. ``@defVar(m, x[i=1:5,j=1:5; i+j >= 3])``
  * Support for vectorized operations on Variables and expressions. See the documentation for details.
  * New ``getVar()`` method to access variables in a model by name
  * Support for semidefinite programming.
  * Dual solutions are now available for general nonlinear problems. You may call ``getDual`` on a reference object for a nonlinear constraint, and ``getDual`` on a variable object for Lagrange multipliers from active bounds.
  * Introduce warnings for two common performance traps: too many calls to ``getValue()`` on a collection of variables and use of the ``+`` operator in a loop to sum expressions.
  * Second-order cone constraints can be written directly with the ``norm()`` and ``norm2{}`` syntax.
  * Implement MathProgBase interface for querying Hessian-vector products.
  * Iteration over ``JuMPContainer``s is deprecated; instead, use the ``keys`` and ``values`` functions, and ``zip(keys(d),values(d))`` for the old behavior.
  * ``@defVar`` returns ``Array{Variable,N}`` when each of ``N`` index sets are of the form ``1:nᵢ``.
  * Module precompilation: on Julia 0.4 and later, ``using JuMP`` is now much faster.

Version 0.9.3 (August 11, 2015)
-------------------------------

  * Fixes for FactCheck testing on julia v0.4.

Version 0.9.2 (June 27, 2015)
-----------------------------

  * Fix bug in @addConstraints.

Version 0.9.1 (April 25, 2015)
------------------------------

  * Fix for Julia 0.4-dev.
  * Small infrastructure improvements for extensions.

Version 0.9.0 (April 18, 2015)
------------------------------

  * Comparison operators for constructing constraints (e.g. ``2x >= 1``) have been deprecated. Instead, construct the constraints explicitly in
    the ``@addConstraint`` macro to add them to the model, or in the ``@LinearConstraint`` macro to create a stand-alone linear constraint instance.
  * ``getValue()`` method implemented to compute the value of a nonlinear subexpression
  * JuMP is now released under the Mozilla Public License version 2.0 (was previously LGPL). MPL is a copyleft license which is less restrictive than LGPL, especially for embedding JuMP within other applications.
  * A number of performance improvements in ReverseDiffSparse for computing derivatives.
  * ``MathProgBase.getsolvetime(m)`` now returns the solution time reported by the solver, if available. (Thanks @odow, Oscar Dowson)
  * Formatting fix for LP format output. (Thanks @sbebo, Leonardo Taccari).

Version 0.8.0 (February 17, 2015)
---------------------------------

  * Nonlinear subexpressions now supported with the ``@defNLExpr`` macro.
  * SCS supported for solving second-order conic problems.
  * ``setXXXCallback`` family deprecated in favor of ``addXXXCallback``.
  * Multiple callbacks of the same type can be registered.
  * Added support for informational callbacks via ``addInfoCallback``.
  * A ``CallbackAbort`` exception can be thrown from callback to safely exit optimization.

Version 0.7.4 (February 4, 2015)
--------------------------------

  * Reduced costs and linear constraint duals are now accessible when quadratic constraints are present.
  * Two-sided nonlinear constraints are supported.
  * Methods for accessing the number of variables and constraints in a model are renamed.
  * New default procedure for setting initial values in nonlinear optimization: project zero onto the variable bounds.
  * Small bug fixes.


Version 0.7.3 (January 14, 2015)
--------------------------------

  * Fix a method ambiguity conflict with Compose.jl (cosmetic fix)

Version 0.7.2 (January 9, 2015)
-------------------------------

  * Fix a bug in ``sum(::JuMPDict)``
  * Added the ``setCategory`` function to change a variables category (e.g. continuous or binary)
  after construction, and ``getCategory`` to retrieve the variable category.

Version 0.7.1 (January 2, 2015)
-------------------------------

  * Fix a bug in parsing linear expressions in macros. Affects only Julia 0.4 and later.

Version 0.7.0 (December 29, 2014)
---------------------------------

### Linear/quadratic/conic programming

  * **Breaking change**: The syntax for column-wise model generation has been changed to use keyword arguments in ``@defVar``.
  * On Julia 0.4 and later, variables and coefficients may be multiplied in any order within macros. That is, variable*coefficient is now valid syntax.
  * ECOS supported for solving second-order conic problems.

### Nonlinear programming

  * Support for skipping model generation when solving a sequence of nonlinear models with changing data.
  * Fix a memory leak when solving a sequence of nonlinear models.
  * The ``@addNLConstraint`` macro now supports the three-argument version to define sets of nonlinear constraints.
  * KNITRO supported as a nonlinear solver.
  * Speed improvements for model generation.
  * The ``@addNLConstraints`` macro supports adding multiple (groups of) constraints at once. Syntax is similar to ``@addConstraints``.
  * Discrete variables allowed in nonlinear problems for solvers which support them (currently only KNITRO).

### General

  * Starting values for variables may now be specified with ``@defVar(m, x, start=value)``.
  * The ``setSolver`` function allows users to change the solver subsequent to model creation.
  * Support for "fixed" variables via the ``@defVar(m, x == 1)`` syntax.
  * Unit tests rewritten to use FactCheck.jl, improved testing across solvers.

Version 0.6.3 (October 19, 2014)
--------------------------------

  * Fix a bug in multiplying two AffExpr objects.


Version 0.6.2 (October 11, 2014)
--------------------------------

  * Further improvements and bug fixes for printing.
  * Fixed a bug in ``@defExpr``.
  * Support for accessing expression graphs through the MathProgBase NLP interface.

Version 0.6.1 (September 19, 2014)
----------------------------------

  * Improvements and bug fixes for printing.

Version 0.6.0 (September 9, 2014)
---------------------------------

  * Julia 0.3.0 is the minimum required version for this release.
  * ``buildInternalModel(m::Model)`` added to build solver-level model in memory without optimizing.
  * Deprecate ``load_model_only`` keyword argument to ``solve``.
  * Add groups of constraints with ``@addConstraints`` macro.
  * Unicode operators now supported, including ``∑`` for ``sum``, ``∏`` for ``prod``, and ``≤``/``≥``
  * Quadratic constraints supported in ``@addConstraint`` macro.
  * Quadratic objectives supported in ``@setObjective`` macro.
  * MathProgBase solver-independent interface replaces Ipopt-specific interface for nonlinear problems
    - **Breaking change**: ``IpoptOptions`` no longer supported to specify solver options, use ``m = Model(solver=IpoptSolver(options...))`` instead.
  * New solver interfaces: ECOS, NLopt, and nonlinear support for MOSEK
  * New option to control whether the lazy constraint callback is executed at each node in the B&B tree or just when feasible solutions are found
  * Add support for semicontinuous and semi-integer variables for those solvers that support them.
  * Add support for index dependencies (e.g. triangular indexing) in ``@defVar``, ``@addConstraint``, and ``@defExpr`` (e.g. ``@defVar(m, x[i=1:10,j=i:10])``).
    - This required some changes to the internal structure of JuMP containers, which may break code that explicitly stored ``JuMPDict`` objects.

Version 0.5.8 (September 24, 2014)
----------------------------------

  * Fix a bug with specifying solvers (affects Julia 0.2 only)

Version 0.5.7 (September 5, 2014)
---------------------------------

  * Fix a bug in printing models

Version 0.5.6 (September 2, 2014)
---------------------------------
  * Add support for semicontinuous and semi-integer variables for those solvers that support them.
    - **Breaking change**: Syntax for ``Variable()`` constructor has changed (use of this interface remains discouraged)
  * Update for breaking changes in MathProgBase

Version 0.5.5 (July 6, 2014)
----------------------------

  * Fix bug with problem modification: adding variables that did not appear in existing constraints or objective.

Version 0.5.4 (June 19, 2014)
----------------------------

  * Update for breaking change in MathProgBase which reduces loading times for ``using JuMP``
  * Fix error when MIPs not solved to optimality


Version 0.5.3 (May 21, 2014)
----------------------------

  * Update for breaking change in ReverseDiffSparse

Version 0.5.2 (May 9, 2014)
---------------------------

  * Fix compatibility with Julia 0.3 prerelease


Version 0.5.1 (May 5, 2014)
---------------------------

  * Fix a bug in coefficient handling inside lazy constraints and user cuts

Version 0.5.0 (May 2, 2014)
---------------------------

  * Support for nonlinear optimization with exact, sparse second-order derivatives automatically computed. Ipopt is currently the only solver supported.
  * ``getValue`` for ``AffExpr`` and ``QuadExpr``
  * **Breaking change**: ``getSolverModel`` replaced by ``getInternalModel``, which returns the internal MathProgBase-level model
  * Groups of constraints can be specified with ``@addConstraint`` (see documentation for details). This is not a breaking change.
  * ``dot(::JuMPDict{Variable},::JuMPDict{Variable})`` now returns the corresponding quadratic expression.

Version 0.4.1 (March 24, 2014)
------------------------------

  * Fix bug where change in objective sense was ignored when re-solving a model.
  * Fix issue with handling zero coefficients in AffExpr.

Version 0.4.0 (March 10, 2014)
------------------------------

  * Support for SOS1 and SOS2 constraints.
  * Solver-independent callback for user heuristics.
  * ``dot`` and ``sum`` implemented for ``JuMPDict`` objects. Now you can say ``@addConstraint(m, dot(a,x) <= b)``.
  * Developers: support for extensions to JuMP. See definition of Model in ``src/JuMP.jl`` for more details.
  * Option to construct the low-level model before optimizing.

Version 0.3.2 (February 17, 2014)
---------------------------------

 * Improved model printing
   - Preliminary support for IJulia output

Version 0.3.1 (January 30, 2014)
--------------------------------

 * Documentation updates
   - Support for MOSEK
   - CPLEXLink renamed to CPLEX

Version 0.3.0 (January 21, 2014)
--------------------------------

 * Unbounded/infeasibility rays: getValue() will return the corresponding
   components of an unbounded ray when a model is unbounded, if supported
   by the selected solver. getDual() will return an infeasibility ray (Farkas proof)
   if a model is infeasible and the selected solver supports this feature.
 * Solver-independent callbacks for user generated cuts.
 * Use new interface for solver-independent QCQP.
 * ``setlazycallback`` renamed to ``setLazyCallback`` for consistency.

Version 0.2.0 (December 15, 2013)
---------------------------------

  * **Breaking change**: Objective sense is specified in setObjective
    instead of in the Model constructor.
  * **Breaking change**: ``lpsolver`` and ``mipsolver`` merged into
    single ``solver`` option.
  * Problem modification with efficient LP restarts and MIP warm-starts.
  * Relatedly, column-wise modeling now supported.
  * Solver-independent callbacks supported. Currently we support only
    a "lazy constraint" callback, which works with Gurobi, CPLEX, and GLPK.
    More callbacks coming soon.

Version 0.1.2 (November 16, 2013)
---------------------------------

  * Bug fixes for printing, improved error messages.
  * Allow ``AffExpr`` to be used in macros; e.g.,
    ``ex = y + z; @addConstraint(m, x + 2*ex <= 3)``

Version 0.1.1 (October 23, 2013)
--------------------------------

  * Update for solver specification API changes in MathProgBase.


Version 0.1.0 (October 3, 2013)
-------------------------------

  * Initial public release.

JuMP release notes
==================

Version 0.20.1 (Oct 18, 2019)
-----------------------------

- Add sections on `@variables` and `@constraints` in the documentation (#2062).
- Fixed product of sparse matrices for Julia v1.3 (#2063).
- Added `set_objective_coefficient` to modify the coefficient of a linear term
  of the objective function (#2008).
- Added `set_time_limit_sec`, `unset_time_limit_sec` and `time_limit_sec` to set
  and query the time limit for the solver in seconds (#2053).

Version 0.20.0 (Aug 24, 2019)
-----------------------------

- Documentation updates.
- Numerous bug fixes.
- Better error messages (#1977, #1978, #1997, #2017).
- Performance improvements (#1947, #2032).
- Added LP sensitivity summary functions `lp_objective_perturbation_range`
  and `lp_rhs_perturbation_range` (#1917).
- Added functions `dual_objective_value`, `raw_status` and `set_parameter`.
- Added function `set_objective_coefficient` to modify the coefficient of
  a linear term of the objective (#2008).
- Added functions `set_normalized_rhs`, `normalized_rhs`, and
  `add_to_function_constant` to modify and get the constant part
  of a constraint (#1935, #1960).
- Added functions `set_normalized_coefficient` and `normalized_coefficient`
  to modify and get the coefficient of a linear term of a constraint
  (#1935, #1960).
- Numerous other improvements in MOI 0.9, see the `NEWS.md` file of MOI for more
  details.

Version 0.19.2 (June 8, 2019)
-----------------------------

- Fix a bug in derivatives that could arise in models with nested nonlinear
  subexpressions.

Version 0.19.1 (May 12, 2019)
-----------------------------

- Usability and performance improvements.
- Bug fixes.

Version 0.19.0 (February 15, 2019)
----------------------------------

**JuMP 0.19 contains significant breaking changes.**

Breaking changes:

- JuMP's abstraction layer for communicating with solvers changed from
  [MathProgBase](https://github.com/JuliaOpt/MathProgBase.jl) (MPB) to
  [MathOptInterface](https://github.com/JuliaOpt/MathOptInterface.jl)
  (MOI). MOI addresses many longstanding design issues. (See @mlubin's
  [slides](http://www.juliaopt.org/meetings/bordeaux2018/lubin.pdf) from
  JuMP-dev 2018.) JuMP 0.19 is compatible only with solvers that have been
  updated for MOI. See the
  [installation guide](http://www.juliaopt.org/JuMP.jl/dev/installation/)
  for a list of solvers that have and have not yet been updated.

- Most solvers have been renamed to `PackageName.Optimizer`. For example,
  `GurobiSolver()` is now `Gurobi.Optimizer`.

- Solvers are no longer added to a model via `Model(solver = XXX(kwargs...))`.
  Instead use `Model(with_optimizer(XXX, kwargs...))`. For example, `Model(with_optimizer(Gurobi.Optimizer, OutputFlag=0))`.

- JuMP containers (e.g., the objects returned by `@variable`) have been
  redesigned. `Containers.SparseAxisArray` replaces `JuMPDict`, `JuMPArray` was
  rewritten (inspired by `AxisArrays`) and renamed `Containers.DenseAxisArray`,
  and you can now request a container type with the `container=` keyword to the
  macros. See the corresponding
  [documentation](http://www.juliaopt.org/JuMP.jl/dev/variables/#Variable-containers-1)
  for more details.

- The statuses returned by solvers have changed. See the possible status
  values
  [here](http://www.juliaopt.org/MathOptInterface.jl/stable/apireference.html#Termination-Status-1).
  The MOI statuses are much richer than the MPB statuses and can be used to
  distinguish between previously indistinguishable cases (e.g. did the solver
  have a feasible solution when it stopped because of the time limit?).

- Starting values are separate from result values. Use `value` to query
  the value of a variable in a solution. Use `start_value` and `set_start_value`
  to get and set an initial starting point provided to the solver. The solutions
  from previous solves are no longer automatically set as the starting points
  for the next solve.

- The data structures for affine and quadratic expressions `AffExpr` and
  `QuadExpr` have changed. Internally, terms are stored in dictionaries instead
  of lists. Duplicate coefficients can no longer exist. Accessors and iteration
  methods have changed.

- `JuMPNLPEvaluator` no longer includes the linear and quadratic parts of the
  model in the evaluation calls. These are now handled separately to allow NLP
  solvers that support various types of constraints.

- JuMP solver-independent callbacks have been replaced by solver-specific
  callbacks. See your favorite solver for more details. (See the note below: No
  solver-specific callbacks are implemented yet.)

- The `norm()` syntax is no longer recognized inside macros. Use the
  `SecondOrderCone()` set instead.

- JuMP no longer performs automatic transformation between special quadratic
  forms and second-order cone constraints. Support for these
  constraint classes depends on the solver.

- The symbols `:Min` and `:Max` are no longer used as optimization senses.
  Instead, JuMP uses the `OptimizationSense` enum from MathOptInterface.
  `@objective(model, Max, ...)`, `@objective(model, Min, ...)`,
  `@NLobjective(model, Max, ...)`, and `@objective(model, Min, ...)` remain
  valid, but `@objective(m, :Max, ...)` is no longer accepted.

- The sign conventions for duals has changed in some cases for consistency with
  conic duality (see the
  [documentation](http://www.juliaopt.org/MathOptInterface.jl/v0.6.2/apimanual.html#Duals-1)).
  The `shadow_price` helper method returns duals with signs that match
  conventional LP interpretations of dual values as sensitivities of the
  objective value to relaxations of constraints.

- `@constraintref` is no longer defined. Instead, create the appropriate
  container to hold constraint references manually. For example,
  ```julia
  constraints = Dict() # Optionally, specify types for improved performance.
  for i in 1:N
    constraints[i] = @constraint(model, ...)
  end
  ```

- The `lowerbound`, `upperbound`, and `basename` keyword arguments to the `@variable`
  macro have been renamed to `lower_bound`, `upper_bound`, and `base_name`,
  for consistency with JuMP's new
  [style recommendations](http://www.juliaopt.org/JuMP.jl/dev/style/).

- We rely on broadcasting syntax to apply accessors to collections of
  variables, e.g., `value.(x)` instead of `getvalue(x)` for collections. (Use
  `value(x)` when `x` is a scalar object.)

New features:

- Splatting (like `f(x...)`) is recognized in restricted settings in nonlinear
  expressions.

- Support for deleting constraints and variables.

- The documentation has been completely rewritten using docstrings and
  Documenter.

- Support for modeling mixed conic and quadratic models (e.g., conic models
  with quadratic objectives and bi-linear matrix inequalities).

- Significantly improved support for modeling new types of constraints and for
  extending JuMP's macros.

- Support for providing dual warm starts.

- Improved support for accessing solver-specific attributes (e.g., the
  irreducible inconsistent subsystem).

- Explicit control of whether symmetry-enforcing constraints are added to PSD
  constraints.

- Support for modeling exponential cones.

- Significant improvements in internal code quality and testing.

- Style and naming guidelines.

- Direct mode and manual mode provide explicit control over when copies of a
  model are stored and/or regenerated. See the corresponding
  [documentation](http://www.juliaopt.org/JuMP.jl/dev/solvers/).

There are known regressions from JuMP 0.18 that will be addressed in a future
release (0.19.x or later):

- Performance regressions in model generation
  ([issue](https://github.com/JuliaOpt/JuMP.jl/issues/1403)). Please file an
  issue anyway if you notice a significant performance regression. We have
  plans to address a number of performance issues, but we might not be aware of
  all of them.

- Fast incremental NLP solves are not yet reimplemented
  ([issue](https://github.com/JuliaOpt/JuMP.jl/issues/1185)).

- We do not yet have an implementation of solver-specific callbacks.

- The column generation syntax in `@variable` has been removed (i.e., the
  `objective`, `coefficients`, and `inconstraints` keyword arguments). Support
  for column generation will be re-introduced in a future release.

- The ability to solve the continuous relaxation (i.e. via
  `solve(model; relaxation = true)`) is not yet reimplemented ([issue](https://github.com/JuliaOpt/JuMP.jl/issues/1611)).

Version 0.18.5 (December 1, 2018)
---------------------------------

   * Support views in some derivative evaluation functions.
   * Improved compatibility with PackageCompiler.

Version 0.18.4 (October 8, 2018)
--------------------------------

   * Fix a bug in model printing on Julia 0.7 and 1.0.

Version 0.18.3 (October 1, 2018)
------------------------------

   * Add support for Julia v1.0 (Thanks @ExpandingMan)
   * Fix matrix expressions with quadratic functions (#1508)

Version 0.18.2 (June 10, 2018)
------------------------------

   * Fix a bug in second-order derivatives when expressions are present (#1319)
   * Fix a bug in `@constraintref` (#1330)

Version 0.18.1 (April 9, 2018)
------------------------------

   * Fix for nested tuple destructuring (#1193)
   * Preserve internal model when relaxation=true (#1209)
   * Minor bug fixes and updates for example

Version 0.18.0 (July 27, 2017)
------------------------------

   * Drop support for Julia 0.5.
   * Update for ForwardDiff 0.5.
   * Minor bug fixes.


Version 0.17.1 (June 9, 2017)
-------------------------------

   * Use of `constructconstraint!` in `@SDconstraint`.
   * Minor bug fixes.

Version 0.17.0 (May 27, 2017)
-------------------------------

   * **Breaking change**: Mixing quadratic and conic constraints is no longer supported.
   * **Breaking change**: The ``getvariable`` and ``getconstraint`` functions are replaced by indexing on the corresponding symbol. For instance, to access the variable with name ``x``, one should now write ``m[:x]`` instead of ``getvariable(m, :x)``. As a consequence, creating a variable and constraint with the same name now triggers a warning, and accessing one of them afterwards throws an error. This change is breaking only in the latter case.
   * Addition of the ``getobjectivebound`` function that mirrors the functionality of the MathProgBase ``getobjbound`` function except that it takes into account transformations performed by JuMP.
   * Minor bug fixes.

The following changes are primarily of interest to developers of JuMP extensions:

   * The new syntax ``@constraint(model, expr in Cone)`` creates the constraint ensuring that ``expr`` is inside ``Cone``. The ``Cone`` argument is passed to ``constructconstraint!`` which enables the call to the dispatched to an extension.
   * The ``@variable`` macro now calls ``constructvariable!`` instead of directly calling the ``Variable`` constructor. Extra arguments and keyword arguments passed to ``@variable`` are passed to ``constructvariable!`` which enables the call to be dispatched to an extension.
   * Refactor the internal function ``conicdata`` (used build the MathProgBase conic model) into smaller subfunctions to make these parts reusable by extensions.



Version 0.16.2 (March 28, 2017)
-------------------------------

   * Minor bug fixes and printing tweaks
   * Address deprecation warnings for Julia 0.6

Version 0.16.1 (March 7, 2017)
------------------------------

   * Better support for ``AbstractArray`` in JuMP (Thanks @tkoolen)
   * Minor bug fixes

Version 0.16.0 (February 23, 2017)
---------------------------

   * **Breaking change**: JuMP no longer has a mechanism for selecting solvers by default (the previous mechanism was flawed and incompatible with Julia 0.6). Not specifying a solver before calling ``solve()`` will result in an error.
   * **Breaking change**: User-defined functions are no longer global. The first argument to ``JuMP.register`` is now a JuMP ``Model`` object within whose scope the function will be registered. Calling ``JuMP.register`` without a ``Model`` now produces an error.
   * **Breaking change**: Use the new ``JuMP.fix`` method to fix a variable to a value or to update the value to which a variable is fixed. Calling ``setvalue`` on a fixed variable now results in an error in order to avoid silent behavior changes. (Thanks @joaquimg)
   * Nonlinear expressions now print out similarly to linear/quadratic expressions (useful for debugging!)
   * New ``category`` keyword to ``@variable``. Used for specifying categories of anonymous variables.
   * Compatibility with Julia 0.6-dev.
   * Minor fixes and improvements (Thanks @cossio, @ccoffrin, @blegat)


Version 0.15.1 (January 31, 2017)
---------------------------------

  * Bugfix for ``@LinearConstraints`` and friends


Version 0.15.0 (December 22, 2016)
----------------------------------

  * Julia 0.5.0 is the minimum required version for this release.
  * Document support for BARON solver
  * Enable info callbacks in more states than before, e.g. for recording solutions.
    New ``when`` argument to ``addinfocallback`` ([#814](https://github.com/JuliaOpt/JuMP.jl/pull/814), thanks @yeesian)
  * Improved support for anonymous variables. This includes new warnings for potentially confusing use of the traditional non-anonymous syntax:
    * When multiple variables in a model are given the same name
    * When non-symbols are used as names, e.g., ``@variable(m, x[1][1:N])``
  * Improvements in iterating over JuMP containers ([#836](https://github.com/JuliaOpt/JuMP.jl/pull/836), thanks @IssamT)
  * Support for writing variable names in .lp file output (Thanks @leethargo)
  * Support for querying duals to SDP problems (Thanks @blegat)
  * The comprehension syntax with curly braces ``sum{}``, ``prod{}``, and ``norm2{}`` has been deprecated
    in favor of Julia's native comprehension syntax ``sum()``, ``prod()`` and ``norm()`` as previously announced.
    (For early adopters of the new syntax, ``norm2()`` was renamed to ``norm()`` without deprecation.)
  * Unit tests rewritten to use Base.Test instead of FactCheck
  * Improved support for operations with matrices of JuMP types (Thanks @ExpandingMan)
  * The syntax to halt a solver from inside a callback has changed from ``throw(CallbackAbort())`` to ``return JuMP.StopTheSolver``
  * Minor bug fixes

Version 0.14.2 (December 12, 2016)
----------------------------------

  * Allow singeton anonymous variables (includes bugfix)

Version 0.14.1 (September 12, 2016)
-----------------------------------

  * More consistent handling of states in informational callbacks,
    includes a new ``when`` parameter to ``addinfocallback`` for
    specifying in which state an informational callback should be called.

Version 0.14.0 (August 7, 2016)
-------------------------------

  * Compatibility with Julia 0.5 and ForwardDiff 0.2
  * Support for "anonymous" variables, constraints, expressions, and parameters, e.g.,
    ``x = @variable(m, [1:N])`` instead of ``@variable(m, x[1:N])``
  * Support for retrieving constraints from a model by name via ``getconstraint``
  * ``@NLconstraint`` now returns constraint references (as expected).
  * Support for vectorized expressions within lazy constraints
  * On Julia 0.5, parse new comprehension syntax ``sum(x[i] for i in 1:N if isodd(i))``
    instead of ``sum{ x[i], i in 1:N; isodd(i) }``. The old syntax with curly
    braces will be deprecated in JuMP 0.15.
  * Now possible to provide nonlinear expressions as "raw" Julia ``Expr`` objects
    instead of using JuMP's nonlinear macros. This input format is useful for
    programmatically generated expressions.
  * ``s/Mathematical Programming/Mathematical Optimization/``
  * Support for local cuts (Thanks to @madanim, Mehdi Madani)
  * Document Xpress interface developed by @joaquimg, Joaquim Dias Garcia
  * Minor bug and deprecation fixes (Thanks @odow, @jrevels)


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

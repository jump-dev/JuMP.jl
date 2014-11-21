JuMP release notes
==================

Unreleased
----------

  * The ``@defConstrRef`` macro is deprecated in favor of the three-argument version of ``@addConstraint``.
  * The ``@addNLConstraint`` macro now supports the three-argument version to define sets of nonlinear constraints.
  * Speed improvements for nonlinear model generation.
  * KNITRO supported as a nonlinear solver.
  * On Julia 0.4, variables and coefficients may be multiplied in any order within macros. That is, variable*coefficient is now valid syntax.

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

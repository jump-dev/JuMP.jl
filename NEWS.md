JuMP release notes
==================

Unversioned
-----------

  * Support for SOS1 and SOS2 constraints.

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

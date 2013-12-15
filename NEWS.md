JuMP release notes
==================

Version 0.2.0 (December 15, 2013)
--------------------------

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

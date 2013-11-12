JuMP release notes
==================

Version 0.2.0 (UNRELEASED)
--------------------------

  * Breaking change: Objective sense is specified in setObjective
    instead of in the Model constructor.
  * Breaking change: ``lpsolver`` and ``mipsolver`` merged into
    single ``solver`` option.
  * Problem modification with efficient LP restarts and MIP warm-starts.
  * Relatedly, column-wise modeling now supported.

Version 0.1.1 (October 23, 2013)
--------------------------------

  * Update for solver specification API changes in MathProgBase.


Version 0.1.0 (October 3, 2013)
-------------------------------

  * Initial public release.

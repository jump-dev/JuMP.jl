```@meta
CurrentModule = JuMP
```

# Release notes

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## Version 1.11.1 (May 19, 2023)

### Fixed

 - Fixed a poor error message when `sum(::DenseAxisArray; dims)` was called
   (#3338)
 - Fixed support for dependent sets in the [`@variable`](@ref) macro (#3344)
 - Fixed a performance bug in constraints with sparse symmetric matrices (#3349)

### Other

 - Improved the printing of complex numbers (#3332)
 - When printing, sets which contain constants ending in `.0` now print as
   integers. This follows the behavior of constants in functions (#3341)
 - Added `InfiniteOpt` to the extensions documentation (#3343)
 - Added more documentation for the exponential cone (#3345) (#3347)
 - Added checklists for developers (#3346) (#3355)
 - Fixed test support upcoming Julia nightly (#3351)
 - Fixed `extension-tests.yml` action (#3353)
 - Add more solvers to the documentation (#3359) (#3361) (#3362)

## Version 1.11.0 (May 3, 2023)

### Added

 - Added new methods to [`print_active_bridges`](@ref) for printing a particular
   objective, constraint, or variable (#3316)

### Fixed

 - Fixed tests for MOI v1.14.0 release (#3312)
 - Fixed indexing containers when an axis is `Vector{Any}` that contains a
   `Vector{Any}` element (#3280)
 - Fixed `getindex(::AbstractJuMPScalar)` which is called for an expression like
   `x[]` (#3314)
 - Fixed bug in `set_string_names_on_creation` with a vector of variables
   (#3322)
 - Fixed bug in `memoize` function in nonlinear documentation (#3337)

### Other

 - Fixed typos in the documentation (#3317) (#3318) (#3328)
 - Added a test for the order of setting start values (#3315)
 - Added READMEs of solvers and extensions to the docs (#3309) (#3320) (#3327)
   (#3329) (#3333)
 - Style improvements to `src/variables.jl` (#3324)
 - Clarify that column generation does not find global optimum (#3325)
 - Add a GitHub actions workflow for testing extensions prior to release (#3331)
 - Document the release process for JuMP (#3334)
 - Fix links to discourse and chatroom (#3335)

## Version 1.10.0 (April 3, 2023)

### Added

 - Added [`Nonnegatives`](@ref), [`Nonpositives`](@ref) and [`Zeros`](@ref), and
   support vector-valued inequality syntax in the JuMP macros (#3273)
 - Added special support for `LinearAlgebra.Symmetric` and `LinearAlgebra.Hermitian`
   matrices in [`Zeros`](@ref) constraints (#3281) (#3296)
 - Added [`HermitianMatrixSpace`](@ref) and the `Hermitian` tag for generating a
   matrix of variables that is Hermitian (#3292) (#3293)
 - Added [`Semicontinuous`](@ref) and [`Semiinteger`](@ref) (#3302)
 - Added support for keyword indexing of containers (#3237)

### Fixed

 - Fixed `[compat]` bound for MathOptInterface in `Project.toml` (#3272)

### Other

 - Split out the [Nested optimization problems](@ref) tutorial (#3274)
 - Updated doctests to ensure none have hidden state (#3275) (#3276)
 - Clarified how lazy constraints may revisit points (#3278)
 - Added [P-Norm](@ref) example (#3282)
 - Clarified docs that macros create new bindings (#3284)
 - Fixed threading example (#3283)
 - Added plot to [The minimum distortion problem](@ref) (#3288)
 - Added Google style rules for Vale and fixed warnings (#3285)
 - Added citation for the JuMP 1.0 paper (#3294)
 - Updated package versions in the documentation (#3298)
 - Added comment for the order in which start values must be set (#3303)
 - Improved error message for unrecognized constraint operators (#3311)

## Version 1.9.0 (March 7, 2023)

### Added

 - Added [`get_attribute`](@ref) and [`set_attribute`](@ref). These replace
   [`get_optimizer_attribute`](@ref) and [`set_optimizer_attribute`](@ref),
   although the `_optimizer_` functions remain for backward compatibility. (#3219)
 - Added [`set_start_values`](@ref) for setting all supported start values in a
   model (#3238)
 - Add [`remove_bridge`](@ref) and [`print_active_bridges`](@ref) (#3259)

### Fixed

 - The matrix returned by a variable in [`HermitianPSDCone`](@ref) is now a
   `LinearAlgebra.Hermitian` matrix. This is potentially breaking if you have
   written code to assume the return is a `Matrix`. (#3245) (#3246)
 - Fixed missing support for `Base.isreal` of expressions (#3252)

### Other

 - Fixed a thread safety issue in the [Parallelism](@ref) tutorial (#3240)
   (#3243)
 - Improved the error message when unsupported operators are used in `@NL`
   macros (#3236)
 - Clarified the documentation to say that matrices in [`HermitianPSDCone`](@ref)
   must be `LinearAlgebra.Hermitian` (#3241)
 - Minor style fixes to internal macro code (#3247)
 - Add [Quantum state discrimination](@ref) tutorial (#3250)
 - Improve error message when `begin...end` not passed to plural macros (#3255)
 - Document how to register function with varying number of input arguments
   (#3258)
 - Tidy tests by removing unneeded `JuMP.` prefixes (#3260)
 - Clarified the introduction to the [Complex number support](@ref) tutorial (#3262)
 - Fixed typos in the Documentation (#3263) (#3266) (#3268) (#3269)

## Version 1.8.2 (February 27, 2023)

### Fixed

 - Fixed dot product between complex JuMP expression and number (#3244)

### Other

 - Polish simple SDP examples (#3232)

## Version 1.8.1 (February 23, 2023)

### Fixed

 - Fixed support for `init` in nonlinear generator expressions (#3226)

### Other

 - Use and document `import MathOptInterface as MOI` (#3222)
 - Removed references in documentation to multiobjective optimization being
   unsupported (#3223)
 - Added tutorial on multi-objective portfolio optimization (#3227)
 - Refactored some of the conic tutorials (#3229)
 - Fixed typos in the documentation (#3230)
 - Added tutorial on parallelism (#3231)

## Version 1.8.0 (February 16, 2023)

### Added

 - Added `-->` syntax support for indicator constraints. The old syntax of `=>`
   remains supported (#3207)
 - Added `<-->` syntax for reified constraints. For now, few solvers support
   reified constraints (#3206)
 - Added [`fix_discrete_variables`](@ref). This is most useful for computing the
   dual of a mixed-integer program (#3208)
 - Added support for vector-valued objectives. For details, see the
   [Multi-objective knapsack](@ref) tutorial (#3176)

### Fixed

 - Fixed a bug in [`lp_sensitivity_report`](@ref) by switching to an explicit
   LU factorization of the basis matrix (#3182)
 - Fixed a bug that prevented `[; kwarg]` arguments in macros (#3220)

### Other

 - Minor fixes to the documentation (#3200) (#3201) (#3203) (#3210)
 - Added tutorial [Constraint programming](@ref) (#3202)
 - Added more examples to [Tips and Tricks](@ref conic_tips_and_tricks)
 - Remove `_distance_to_set` in favor of `MOI.Utilities.distance_to_set` (#3209)
 - Improve [The diet problem](@ref) tutorial by adding the variable as a column
   in the dataframe (#3213)
 - Improve [The knapsack problem example](@ref) tutorial (#3216) (#3217)
 - Added the [Ellipsoid approximation](@ref) tutorial (#3218)

## Version 1.7.0 (January 25, 2023)

### Added

 - Added support for `view` of a [`Containers.DenseAxisArray`](@ref) (#3152) (#3180)
 - Added support for containers of variables in [`ComplexPlane`](@ref) (#3184)
 - Added support for `minimum` and `maximum` generators in nonlinear expressions
   (#3189)
 - Added `SnoopPrecompile` statements that reduce the time-to-first-solve in
   Julia 1.9 (#3193) (#3195) (#3196) (#3197)

### Other

 - Large refactoring of the tests (#3166) (#3167) (#3168) (#3169) (#3170) (#3171)
 - Remove unreachable code due to `VERSION` checks (#3172)
 - Document how to test JuMP extensions (#3174)
 - Fix method ambiguities in `Containers` (#3173)
 - Improve error message that is thrown when `=` is used instead of `==` in
   the [`@constraint`](@ref) macro (#3178)
 - Improve the error message when `Bool` is used instead of `Bin` in the
   [`@variable`](@ref) macro (#3180)
 - Update versions of the documentation (#3185)
 - Tidy the import of packages and remove unnecessary prefixes (#3186) (#3187)
 - Refactor `src/JuMP.jl` by moving methods into more relevant files (#3188)
 - Fix docstring of `Model` not appearing in the documentation (#3198)

## Version 1.6.0 (January 1, 2023)

### Added

 - Added a `result` keyword argument to [`solution_summary`](@ref) to allow
   summarizing models with multiple solutions (#3138)
 - Added [`relax_with_penalty!`](@ref), which is a useful tool when debugging
   infeasible models (#3140)
 - Added [`has_start_value`](@ref) (#3157)
 - Added support for [`HermitianPSDCone`](@ref) in constraints (#3154)

### Fixed

 - Fixed promotion of complex expressions (#3150) (#3164)

### Other

 - Added Benders tutorial with in-place resolves (#3145)
 - Added more [Tips and tricks](@ref linear_tips_and_tricks) for linear programs
   (#3144) (#3163)
 - Clarified documentation that `start` can depend on the indices of a
   variable container (#3148)
 - Replace instances of `length` and `size` by the recommended `eachindex` and
   `axes` (#3149)
 - Added a warning explaining why the model is dirty when accessing solution
   results from a modified model (#3156)
 - Clarify documentation that `PSD` ensures a symmetric matrix (#3159)
 - Maintenance of the JuMP test suite (#3146) (#3158) (#3162)

## Version 1.5.0 (December 8, 2022)

### Added

  - Add support for complex-valued variables:
    * [`HermitianPSDCone`](@ref) (#3109)
    * [`ComplexPlane`](@ref) and [`ComplexVariable`](@ref) (#3134)
  - Add support for [`MOI.OptimizerWithAttributes`](@ref) in
    [`set_optimizer_attribute`](@ref) and [`get_optimizer_attribute`](@ref) (#3129)

### Fixed

  - Fixed error message for vectorized interval constraints (#3123)
  - Fixed passing `AbstractString` to [`set_optimizer_attribute`](@ref) (#3127)

### Other

  - Update package versions used in docs (#3119) (#3133) (#3139)
  - Fixed output of diet tutorial (#3120)
  - Explain how to use `Dates.period` in [`set_time_limit_sec`](@ref) (#3121)
  - Update to `JuliaFormatter` v1.0.15 (#3130)
  - Fixed `HTTP` server example in web_app.jl (#3131)
  - Update docs to build with `Documenter#master` (#3094)
  - Add tests for `LinearAlgebra` operations (#3132)
  - Tidy these release notes (#3135)
  - Added documentation for [Complex number support](@ref) (#3141)
  - Removed the "workforce scheduling" and "steelT3" tutorials (#3143)

## Version 1.4.0 (October 29, 2022)

### Added

  - Added [`Containers.rowtable`](@ref) which converts a container into a vector
    of `NamedTuple`s to support the Tables.jl interface. This simplifies
    converting [`Containers.DenseAxisArray`](@ref) and [`Containers.SparseAxisArray`](@ref)
    objects into tabular forms such as a DataFrame (#3104)
  - Added a new method to [`Containers.container`](@ref) so that index names are
    passed to the container (#3088)

### Fixed

  - Fixed a bug in `copy_to(dest::Model, src::MOI.ModelLike)` when `src` has
    nonlinear components (#3101)
  - Fixed the printing of `(-1.0 + 0.0im)` coefficients in complex expressions
    (#3112)
  - Fixed a parsing bug in nonlinear expressions with generator statements that
    contain multiple `for` statements (#3116)

### Other

  - Converted the multi-commodity flow tutorial to use an SQLite database (#3098)
  - Fixed a number of typos in the documentation (#3103) (#3107) (#3018)
  - Improved various style aspects of the PDF documentation (#3095) (#3098)
    (#3102)

## Version 1.3.1 (September 28, 2022)

### Fixed

  - Fixed a performance issue in `relax_integrality` (#3087)
  - Fixed the type stability of operators with `Complex` arguments (#3072)
  - Fixed a bug which added additional `+()` terms to some nonlinear expressions
    (#3091)
  - Fixed potential method ambiguities with `AffExpr` and `QuadExpr` objects
    (#3092)

### Other

  - Added [vale](https://vale.sh) as a linter for the documentation (#3080)
  - Added a tutorial on debugging JuMP models (#3043)
  - Fixed a number of typos in the documentation (#3079) (#3083)
  - Many other small tweaks to the documentation (#3068) (#3073) (#3074) (#3075)
    (#3076) (#3077) (#3078) (#3081) (#3082) (#3084) (#3085) (#3089)

## Version 1.3.0 (September 5, 2022)

### Added

  - Support slicing in `SparseAxisArray` (#3031)

### Fixed

  - Fixed a bug introduced in v1.2.0 that prevented `DenseAxisArray`s with
    `Vector` keys (#3064)

### Other

  - Released the JuMP logos under the CC BY 4.0 license (#3063)
  - Minor tweaks to the documentation (#3054) (#3056) (#3057) (#3060) (#3061)
    (#3065)
  - Improved code coverage of a number of files (#3048) (#3049) (#3050) (#3051)
    (#3052) (#3053) (#3058) (#3059)

## Version 1.2.1 (August 22, 2022)

### Fixed

  - Fixed a bug when parsing two-sided nonlinear constraints (#3045)

## Version 1.2.0 (August 16, 2022)

### Breaking

This is a large minor release because it significantly refactors the internal
code for handling nonlinear programs to use the `MathOptInterface.Nonlinear`
submodule that was introduced in MathOptInterface v1.3.0. As a consequence, the
internal datastructure in `model.nlp_data` has been removed, as has the
`JuMP._Derivatives` submodule. Despite the changes, the public API for nonlinear
programming has not changed, and any code that uses only the public API and that
worked with v1.1.1 will continue to work with v1.2.0.

### Added

  - Added `all_constraints(model; include_variable_in_set_constraints)` which
    simplifies returning a list of all constraint indices in the model.
  - Added the ability to delete nonlinear constraints via
    `delete(::Model, ::NonlinearConstraintRef)`.
  - Added the ability to provide an explicit Hessian for a multivariate
    user-defined function.
  - Added support for querying the primal value of a nonlinear constraint via
    `value(::NonlinearConstraintRef)`

### Fixed

  - Fixed a bug in `Containers.DenseAxisArray` so that it now supports indexing
    with keys that hash to the same value, even if they are different types, for
    example, `Int32` and `Int64`.
  - Fixed a bug printing the model when the solver does not support `MOI.Name`.

### Other

  - Added a constraint programming formulation to the Sudoku tutorial.
  - Added newly supported solvers Pajarito, Clarabel, and COPT to the
    installation table.
  - Fixed a variety of other miscellaneous issues in the documentation.

## Version 1.1.1 (June 14, 2022)

### Other

  - Fixed problem displaying LaTeX in the documentation
  - Minor updates to the style guide
  - Updated to MOI v1.4.0 in the documentation

## Version 1.1.0 (May 25, 2022)

### Added

  - Added `num_constraints(::Model; count_variable_in_set_constraints)` to
    simplify the process of counting the number of constraints in a model
  - Added `VariableRef(::ConstraintRef)` for querying the variable associated
    with a bound or integrality constraint.
  - Added `set_normalized_coefficients` for modifying the variable coefficients
    of a vector-valued constraint.
  - Added `set_string_names_on_creation` to disable creating `String` names for
    variables and constraints. This can improve performance.

### Fixed

  - Fixed a bug passing `nothing` to the `start` keyword of `@variable`

### Other

  - New tutorials:
    - Sensitivity analysis of a linear program
    - Serving web apps
  - Minimal ellipse SDP tutorial refactored and improved
  - Docs updated to the latest version of each package
  - Lots of minor fixes and improvements to the documentation

## Version 1.0.0 (March 24, 2022)

**Read more about this release, along with an acknowledgement of all the
contributors in our [JuMP 1.0.0 is released](https://jump.dev/blog/1.0.0-release/)
blog post.**

### Breaking

  - The previously deprecated functions (v0.23.0, v0.23.1) have been removed.
    Deprecation was to improve consistency of function names:
    - `num_nl_constraints` (see `num_nonlinear_constraints`)
    - `all_nl_constraints` (see `all_nonlinear_constraints`)
    - `add_NL_expression` (see `add_nonlinear_expression`)
    - `set_NL_objective`  (see `set_nonlinear_objective`)
    - `add_NL_constraint` (see `add_nonlinear_constraint`)
    - `nl_expr_string` (see `nonlinear_expr_string`)
    - `nl_constraint_string` (see `nonlinear_constraint_string`)
    - `SymMatrixSpace` (see `SymmetricMatrixSpace`)
  - The unintentionally exported variable `JuMP.op_hint` has been renamed to the
    unexported `JuMP._OP_HINT`

### Fixed

  - Fixed a bug writing .nl files
  - Fixed a bug broadcasting `SparseAxisArray`s

## Version 0.23.2 (March 14, 2022)

### Added

  - Added `relative_gap` to `solution_summary`
  - `register` now throws an informative error if the function is not
    differentiable using ForwardDiff. In some cases, the check in `register`
    will encounter a false negative, and the informative error will be thrown at
    run-time. This usually happens when the function is non-differentiable in a
    subset of the domain.

### Fixed

  - Fixed a scoping issue when extending the `container` keyword of containers

### Other

  - Docs updated to the latest version of each package

## Version 0.23.1 (March 2, 2022)

### Deprecated

  - `nl_expr_string` and `nl_constraint_string` have been renamed to
    `nonlinear_expr_string` and `nonlinear_constraint_string`. The old methods
    still exist with deprecation warnings. This change should impact very few
    users because to call them you must rely on private internals of the
    nonlinear API. Users are encouraged to use `sprint(show, x)` instead, where
    `x` is the nonlinear expression or constraint of interest.

### Added

  - Added support for `Base.abs2(x)` where `x` is a variable or affine
    expression. This is mainly useful for complex-valued constraints.

### Fixed

  - Fixed addition of complex and real affine expressions
  - Fixed arithmetic for Complex-valued quadratic expressions
  - Fixed variable bounds passed as `Rational{Int}(Inf)`
  - Fixed printing of the coefficient `(0 + 1im)`
  - Fixed a bug when `solution_summary` is called prior to `optimize!`

## Version 0.23.0 (February 25, 2022)

**JuMP v0.23.0 is a breaking release. It is also a release-candidate for JuMP
v1.0.0. That is, if no issues are found with the v0.23.0 release, then it will
be re-tagged as v1.0.0.**

### Breaking

  - Julia 1.6 is now the minimum supported version
  - MathOptInterface has been updated to v1.0.0
  - All previously deprecated functionality has been removed
  - `PrintMode`, `REPLMode` and `IJuliaMode` have been removed in favor of the
    MIME types `MIME"text/plain"` and `MIME"text/latex"`. Replace instances of
    `::Type{REPLMode}` with `::MIME"text/plain"`, `REPLMode` with
    `MIME("text/plain")`, `::Type{IJuliaMode}` with `::MIME"text/latex"`, and
    `IJuliaMode` with `MIME("text/latex")`.
  - Functions containing the `nl_` acronym have been renamed to the more
    explicit `nonlinear_`. For example, `num_nl_constraints` is now
    `num_nonlinear_constraints` and `set_NL_objective` is now
    `set_nonlinear_objective`. Calls to the old functions throw an error
    explaining the new name.
  - `SymMatrixSpace` has been renamed to `SymmetricMatrixSpace`

### Added

  - Added `nonlinear_dual_start_value` and `set_nonlinear_dual_start_value`
  - Added preliminary support for `Complex` coefficient types

### Fixed

  - Fixed a bug in `solution_summary`

### Other

  - MILP examples have been migrated from GLPK to HiGHS
  - Fixed various typos
  - Improved section on setting constraint start values

### Troubleshooting problems when updating

If you experience problems when updating, you are likely using previously
deprecated functionality. (By default, Julia does not warn when you use
deprecated features.)

To find the deprecated features you are using, start Julia with `--depwarn=yes`:
```
$ julia --depwarn=yes
```
Then install JuMP v0.22.3:
```julia
julia> using Pkg
julia> pkg"add JuMP@0.22.3"
```
And then run your code. Apply any suggestions, or search the release notes below
for advice on updating a specific deprecated feature.

## Version 0.22.3 (February 10, 2022)

### Fixed

  - Fixed a reproducibility issue in the TSP tutorial
  - Fixed a reproducibility issue in the `max_cut_sdp` tutorial
  - Fixed a bug broadcasting an empty SparseAxisArray

### Other

  - Added a warning and improved documentation for the modify-then-query case
  - Fixed a typo in the docstring of `RotatedSecondOrderCone`
  - Added Aqua.jl as a check for code health
  - Added introductions to each section of the tutorials
  - Improved the column generation and Benders decomposition tutorials
  - Updated documentation to MOI v0.10.8
  - Updated JuliaFormatter to v0.22.2

## Version 0.22.2 (January 10, 2022)

### Added

  - The function `all_nl_constraints` now returns all nonlinear constraints
    in a model
  - `start_value` and `set_start_value` can now be used to get and set the
    primal start for constraint references
  - Plural macros now return a tuple containing the elements that were defined
    instead of `nothing`
  - Anonymous variables are now printed as `_[i]` where `i` is the index of the
    variable instead of `noname`. Calling `name(x)` still returns `""` so this
    is non-breaking.

### Fixed

  - Fixed handling of `min` and `max` in nonlinear expressions
  - CartesianIndex is no longer allowed as a key for DenseAxisArrays.

### Other

  - Improved the performance of GenericAffExpr
  - Added a tutorial on the Travelling Salesperson Problem
  - Added a tutorial on querying the Hessian of a nonlinear program
  - Added documentation on using custom solver binaries.

## Version 0.22.1 (November 29, 2021)

### Added

  * Export `OptimizationSense` enum, with instances: `MIN_SENSE`, `MAX_SENSE`,
    and `FEASIBILITY_SENSE`
  * Add `Base.isempty(::Model)` to match `Base.empty(::Model)`

### Fixed

  * Fix bug in container with tuples as indices
  * Fix bug in `set_time_limit_sec`

### Other

  * Add tutorial "Design patterns for larger models"
  * Remove release notes section from PDF
  * General edits of the documentation and error messages

## Version 0.22.0 (November 10, 2021)

**JuMP v0.22 is a breaking release**

### Breaking

JuMP 0.22 contains a number of breaking changes. However, these should be
invisible for the majority of users. You will mostly encounter these breaking
changes if you: wrote a JuMP extension, accessed `backend(model)`, or called
`@SDconstraint`.

The breaking changes are as follows:

 * MathOptInterface has been updated to v0.10.4. For users who have interacted
   with the MOI backend, this contains a large number of breaking changes. Read
   the [MathOptInterface release notes](https://jump.dev/MathOptInterface.jl/v0.10/release_notes/#v0.10.0-(September-6,-2021))
   for more details.
 * The `bridge_constraints` keyword argument to `Model` and `set_optimizer` has
   been renamed `add_bridges` to reflect that more thing were bridged than just
   constraints.
 * The `backend(model)` field now contains a concrete instance of a
   `MOI.Utilities.CachingOptimizer` instead of one with an abstractly typed
   optimizer field. In most cases, this will lead to improved performance.
   However, calling `set_optimizer` after `backend` invalidates the old
   backend. For example:
   ```julia
   model = Model()
   b = backend(model)
   set_optimizer(model, GLPK.Optimizer)
   @variable(model, x)
   # b is not updated with `x`! Get a new b by calling `backend` again.
   new_b = backend(model)
   ```
 * All usages of `@SDconstraint` are deprecated. The new syntax is
   `@constraint(model, X >= Y, PSDCone())`.
 * Creating a `DenseAxisArray` with a `Number` as an axis will now display a
   warning. This catches a common error in which users write
   `@variable(model, x[length(S)])` instead of
   `@variable(model, x[1:length(S)])`.
 * The `caching_mode` argument to `Model`, for example,
   `Model(caching_mode = MOIU.MANUAL)` mode has been removed. For more control
   over the optimizer, use `direct_model` instead.
 * The previously deprecated `lp_objective_perturbation_range` and
   `lp_rhs_perturbation_range` functions have been removed. Use
   `lp_sensitivity_report` instead.
 * The `.m` fields of `NonlinearExpression` and `NonlinearParameter` have been
   renamed to `.model`.
 * Infinite variable bounds are now ignored. Thus, `@variable(model, x <= Inf)`
   will show `has_upper_bound(x) == false`. Previously, these bounds were passed
   through to the solvers which caused numerical issues for solvers expecting
   finite bounds.
 * The `variable_type` and `constraint_type` functions were removed. This should
   only affect users who previously wrote JuMP extensions. The functions can be
   deleted without consequence.
 * The internal functions `moi_mode`, `moi_bridge_constraints`,
   `moi_add_constraint`, and `moi_add_to_function_constant` are no longer
   exported.
 * The un-used method `Containers.generate_container` has been deleted.
 * The `Containers` API has been refactored, and `_build_ref_sets` is now
   public as `Containers.build_ref_sets`.
 * The `parse_constraint_` methods for extending `@constraint` at parse time
   have been refactored in a breaking way. Consult the Extensions documentation
   for more details and examples.

### Added

 * The `TerminationStatusCode` and `ResultStatusCode` enums are now exported
   by JuMP. Prefer `termination_status(model) == OPTIMAL` instead of
   `== MOI.OPTIMAL`, although the `MOI.` prefix way still works.
 * Copy a `x::DenseAxisArray` to an `Array` by calling `Array(x)`.
 * `NonlinearExpression` is now a subtype of `AbstractJuMPScalar`
 * Constraints such as `@constraint(model, x + 1 in MOI.Integer())` are now
   supported.
 * `primal_feasibility_report` now accepts a function as the first argument.
 * Scalar variables `@variable(model, x[1:2] in MOI.Integer())` creates two
   variables, both of which are constrained to be in the set `MOI.Integer`.
 * Conic constraints can now be specified as inequalities under a different
   partial ordering. So `@constraint(model, x - y in MOI.Nonnegatives())` can
   now be written as `@constraint(model, x >= y, MOI.Nonnegatives())`.
 * Names are now set for vectorized constraints.

### Fixed

 * Fixed a performance issue when `show` was called on a `SparseAxisArray` with
   a large number of elements.
 * Fixed a bug displaying barrier and simplex iterations in `solution_summary`.
 * Fixed a bug by implementing `hash` for `DenseAxisArray` and
   `SparseAxisArray`.
 * Names are now only set if the solver supports them. Previously, this
   prevented solvers such as Ipopt from being used with `direct_model`.
 * `MutableArithmetics.Zero` is converted into a `0.0` before being returned to
   the user. Previously, some calls to `@expression` would return the
   undocumented `MutableArithmetics.Zero()` object. One example is summing over
   an empty set `@expression(model, sum(x[i] for i in 1:0))`. You will now get
   `0.0` instead.
 * `AffExpr` and `QuadExpr` can now be used with `== 0` instead of `iszero`.
   This fixes a number of issues relating to Julia standard libraries such as
   `LinearAlgebra` and `SparseArrays`.
 * Fixed a bug when registering a user-defined function with splatting.

### Other

 * The documentation is now available as a PDF.
 * The documentation now includes a full copy of the MathOptInterface
   documentation to make it easy to link concepts between the docs. (The
   MathOptInterface documentation has also been significantly improved.)
 * The documentation contains a large number of improvements and clarifications
   on a range of topics. Thanks to @sshin23, @DilumAluthge, and @jlwether.
 * The documentation is now built with Julia 1.6 instead of 1.0.
 * Various error messages have been improved to be more readable.

## Version 0.21.10 (September 4, 2021)

### Added

  * Added `add_NL_expression`
  * `add_NL_xxx` functions now support `AffExpr` and `QuadExpr` as terms

### Fixed

  * Fixed a bug in `solution_summary`
  * Fixed a bug in `relax_integrality`

### Other

  * Improved error message in `lp_sensitivity_report`

## Version 0.21.9 (August 1, 2021)

### Added

  * Containers now support arbitrary container types by passing the type to the
    `container` keyword and overloading `Containers.container`.
  * `is_valid` now supports nonlinear constraints
  * Added `unsafe_backend` for querying the inner-most optimizer of a JuMP
    model.
  * Nonlinear parameters now support the plural `@NLparameters` macro.
  * Containers (for example, `DenseAxisArray`) can now be used in vector-valued
    constraints.

### Other

  * Various improvements to the documentation.

## Version 0.21.8 (May 8, 2021)

### Added

  * The `@constraint` macro is now extendable in the same way as `@variable`.
  * `AffExpr` and `QuadExpr` can now be used in nonlinear macros.

### Fixed

  * Fixed a bug in `lp_sensitivity_report`.
  * Fixed an inference issue when creating empty `SparseAxisArray`s.

## Version 0.21.7 (April 12, 2021)

### Added

  * Added `primal_feasibility_report`, which can be used to check whether a
    primal point satisfies primal feasibility.
  * Added `coefficient`, which returns the coefficient associated with a
    variable in affine and quadratic expressions.
  * Added `copy_conflict`, which returns the IIS of an infeasible model.
  * Added `solution_summary`, which returns (and prints) a struct containing a
    summary of the solution.
  * Allow `AbstractVector` in vector constraints instead of just `Vector`.
  * Added `latex_formulation(model)` which returns an object representing the
    latex formulation of a model. Use `print(latex_formulation(model))` to print
    the formulation as a string.
  * User-defined functions in nonlinear expressions are now automatically
    registered to aid quick model prototyping. However, a warning is printed to
    encourage the manual registration.
  * DenseAxisArray's now support broadcasting over multiple arrays.
  * Container indices can now be iterators of `Base.SizeUnknown`.

### Fixed

  * Fixed bug in `rad2deg` and `deg2rad` in nonlinear expressions.
  * Fixed a MethodError bug in `Containers` when forcing container type.
  * Allow partial slicing of a DenseAxisArray, resolving an issue from 2014.
  * Fixed a bug printing variable names in IJulia.
  * Ending an IJulia cell with `model` now prints a summary of the model (like
    in the REPL) not the latex formulation. Use `print(model)` to print the latex
    formulation.
  * Fixed a bug when copying models containing nested arrays.

### Other

  * Tutorials are now part of the documentation, and more refactoring has taken
    place.
  * Added JuliaFormatter added as a code formatter.
  * Added some precompilation statements to reduce initial latency.
  * Various improvements to error messages to make them more helpful.
  * Improved performance of `value(::NonlinearExpression)`.
  * Improved performance of `fix(::VariableRef)`.

## Version 0.21.6 (January 29, 2021)

### Added

  * Added support for skew symmetric variables via
    `@variable(model, X[1:2, 1:2] in SkewSymmetricMatrixSpace())`.
  * `lp_sensitivity_report` has been added which significantly improves the
    performance of querying the sensitivity summary of an LP.
    `lp_objective_perturbation_range` and `lp_rhs_perturbation_range` are
    deprecated.
  * Dual warm-starts are now supported with `set_dual_start_value` and
    `dual_start_value`.
  * `∈` (`\in<tab>`) can now be used in macros instead of `=` or `in`.
  * Use `haskey(model::Model, key::Symbol)` to check if a name `key` is
    registered in a model.
  * Added `unregister(model::Model, key::Symbol)` to unregister a name `key`
    from `model`.
  * Added `callback_node_status` for use in callbacks.
  * Added `print_bridge_graph` to visualize the bridging graph generated by
    MathOptInterface.
  * Improved error message for containers with duplicate indices.

### Fixed

  * Various fixes to pass tests on Julia 1.6.
  * Fixed a bug in the printing of nonlinear expressions in IJulia.
  * Fixed a bug when nonlinear expressions are passed to user-defined functions.
  * Some internal functions that were previously exported are now no longer
    exported.
  * Fixed a bug when relaxing a fixed binary variable.
  * Fixed a `StackOverflowError` that occurred when `SparseAxisArray`s had a
    large number of elements.
  * Removed an unnecessary type assertion in `list_of_constraint_types`.
  * Fixed a bug when copying models with registered expressions.

### Other

  * The documentation has been significantly overhauled. It now has distinct
    sections for the manual, API reference, and examples. The existing examples
    in `/examples` have now been moved to `/docs/src/examples` and rewritten
    using `Literate.jl`, and they are now included in the documentation.
  * JuliaFormatter has been applied to most of the codebase. This will continue
    to roll out over time, as we fix upstream issues in the formatter, and will
    eventually become compulsory.
  * The root cause of a large number of method invalidations has been resolved.
  * We switched continuous integration from Travis and Appveyor to GitHub
    Actions.

## Version 0.21.5 (September 18, 2020)

### Fixed

  * Fixed deprecation warnings
  * Throw `DimensionMismatch` for incompatibly sized functions and sets
  * Unify treatment of `keys(x)` on JuMP containers

## Version 0.21.4 (September 14, 2020)

### Added

  * Add debug info when adding unsupported constraints
  * Add `relax_integrality` for solving continuous relaxation
  * Allow querying constraint conflicts

### Fixed

  * Dispatch on `Real` for `MOI.submit`
  * Implement `copy` for `CustomSet` in tests
  * Don't export private macros
  * Fix invalid assertion in nonlinear
  * Error if constraint has `NaN` right-hand side
  * Improve speed of tests
  * Lots of work modularizing files in `/test`
  * Improve line numbers in macro error messages
  * Print nonlinear subexpressions
  * Various documentation updates
  * Dependency updates:
    * Datastructures 0.18
    * MathOptFormat v0.5
    * Prep for MathOptInterface 0.9.15

## Version 0.21.3 (June 18, 2020)

- Added Special Order Sets (SOS1 and SOS2) to JuMP with default weights to ease
  the creation of such constraints (#2212).
- Added functions `simplex_iterations`, `barrier_iterations` and `node_count`
  (#2201).
- Added function `reduced_cost` (#2205).
- Implemented `callback_value` for affine and quadratic expressions (#2231).
- Support `MutableArithmetics.Zero` in objective and constraints (#2219).
- Documentation improvements:
  * Mention tutorials in the docs (#2223).
  * Update COIN-OR links (#2242).
  * Explicit link to the documentation of `MOI.FileFormats` (#2253).
  * Typo fixes (#2261).
- Containers improvements:
  * Fix `Base.map` for `DenseAxisArray` (#2235).
  * Throw `BoundsError` if number of indices is incorrect for `DenseAxisArray`
    and `SparseAxisArray` (#2240).
- Extensibility improvements:
  * Implement a `set_objective` method fallback that redirects to
    `set_objective_sense` and `set_objective_function` (#2247).
  * Add `parse_constraint` method with arbitrary number of arguments (#2051).
  * Add `parse_constraint_expr` and `parse_constraint_head` (#2228).

## Version 0.21.2 (April 2, 2020)

- Added `relative_gap()` to access `MOI.RelativeGap()` attribute (#2199).
- Documentation fixes:
  * Added link to source for docstrings in the documentation (#2207).
  * Added docstring for `@variables` macro (#2216).
  * Typo fixes (#2177, #2184, #2182).
- Implementation of methods for Base functions:
  * Implemented `Base.empty!` for `JuMP.Model` (#2198).
  * Implemented `Base.conj` for JuMP scalar types (#2209).

### Fixed

  * Fixed sum of expression with scalar product in macro (#2178).
  * Fixed writing of nonlinear models to MathOptFormat (#2181).
  * Fixed construction of empty SparseAxisArray (#2179).
  * Fixed constraint with zero function (#2188).

## Version 0.21.1 (Feb 18, 2020)

- Improved the clarity of the `with_optimizer` deprecation warning.

## Version 0.21.0 (Feb 16, 2020)

### Breaking

- Deprecated `with_optimizer` (#2090, #2084, #2141). You can replace
  `with_optimizer` by either nothing, `optimizer_with_attributes` or a closure:
  * replace `with_optimizer(Ipopt.Optimizer)` by `Ipopt.Optimizer`.
  * replace `with_optimizer(Ipopt.Optimizer, max_cpu_time=60.0)`
    by `optimizer_with_attributes(Ipopt.Optimizer, "max_cpu_time" => 60.0)`.
  * replace `with_optimizer(Gurobi.Optimizer, env)` by `() -> Gurobi.Optimizer(env)`.
  * replace `with_optimizer(Gurobi.Optimizer, env, Presolve=0)`
    by `optimizer_with_attributes(() -> Gurobi.Optimizer(env), "Presolve" => 0)`.

  alternatively to `optimizer_with_attributes`, you can also set the attributes
  separately with `set_optimizer_attribute`.
- Renamed `set_parameter` and `set_parameters` to `set_optimizer_attribute` and
  `set_optimizer_attributes` (#2150).
- Broadcast should now be explicit inside macros. `@SDconstraint(model, x >= 1)`
  and `@constraint(model, x + 1 in SecondOrderCone())` now throw an error
  instead of broadcasting `1` along the dimension of `x` (#2107).
- `@SDconstraint(model, x >= 0)` is now equivalent to `@constraint(model, x in PSDCone())`
  instead of `@constraint(model, (x .- 0) in PSDCone())` (#2107).
- The macros now create the containers with `map` instead of `for` loops,
  as a consequence, containers created by `@expression` can now have any element
  type and containers of constraint references now have concrete element types
  when possible. This fixes a long-standing issue where `@expression` could
  only be used to generate a collection of linear expressions. Now it works for
  quadratic expressions as well (#2070).
- Calling `deepcopy(::AbstractModel)` now throws an error.
- The constraint name is now printed in the model string (#2108).

### Added

- Added support for solver-independent and solver-specific callbacks (#2101).
- Added `write_to_file` and `read_from_file`, supported formats are CBF, LP,
  MathOptFormat, MPS and SDPA (#2114).
- Added support for complementarity constraints (#2132).
- Added support for indicator constraints (#2092).
- Added support for querying multiple solutions with the `result` keyword (#2100).
- Added support for constraining variables on creation (#2128).
- Added method `delete` that deletes a vector of variables at once if it is
  supported by the underlying solver (#2135).
- The arithmetic between JuMP expression has be refactored into the
  MutableArithmetics package (#2107).
- Improved error on complex values in NLP (#1978).
- Added an example of column generation (#2010).

### Fixed

- Incorrect coefficients generated when using Symmetric variables (#2102)

## Version 0.20.1 (Oct 18, 2019)

- Add sections on `@variables` and `@constraints` in the documentation (#2062).
- Fixed product of sparse matrices for Julia v1.3 (#2063).
- Added `set_objective_coefficient` to modify the coefficient of a linear term
  of the objective function (#2008).
- Added `set_time_limit_sec`, `unset_time_limit_sec` and `time_limit_sec` to set
  and query the time limit for the solver in seconds (#2053).

## Version 0.20.0 (Aug 24, 2019)

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

## Version 0.19.2 (June 8, 2019)

- Fix a bug in derivatives that could arise in models with nested nonlinear
  subexpressions.

## Version 0.19.1 (May 12, 2019)

- Usability and performance improvements.
- Bug fixes.

## Version 0.19.0 (February 15, 2019)

**JuMP 0.19 contains significant breaking changes.**

### Breaking

- JuMP's abstraction layer for communicating with solvers changed from
  [MathProgBase](https://github.com/JuliaOpt/MathProgBase.jl) (MPB) to
  [MathOptInterface](https://github.com/JuliaOpt/MathOptInterface.jl)
  (MOI). MOI addresses many longstanding design issues. (See @mlubin's
  [slides](https://www.juliaopt.org/meetings/bordeaux2018/lubin.pdf) from
  JuMP-dev 2018.) JuMP 0.19 is compatible only with solvers that have been
  updated for MOI. See the
  [installation guide](https://www.juliaopt.org/JuMP.jl/dev/installation/)
  for a list of solvers that have and have not yet been updated.

- Most solvers have been renamed to `PackageName.Optimizer`. For example,
  `GurobiSolver()` is now `Gurobi.Optimizer`.

- Solvers are no longer added to a model via `Model(solver = XXX(kwargs...))`.
  Instead use `Model(with_optimizer(XXX, kwargs...))`. For example, `Model(with_optimizer(Gurobi.Optimizer, OutputFlag=0))`.

- JuMP containers (for example, the objects returned by `@variable`) have been
  redesigned. `Containers.SparseAxisArray` replaces `JuMPDict`, `JuMPArray` was
  rewritten (inspired by `AxisArrays`) and renamed `Containers.DenseAxisArray`,
  and you can now request a container type with the `container=` keyword to the
  macros. See the corresponding
  [documentation](https://www.juliaopt.org/JuMP.jl/dev/variables/#Variable-containers-1)
  for more details.

- The statuses returned by solvers have changed. See the possible status
  values
  [here](https://www.juliaopt.org/MathOptInterface.jl/stable/apireference.html#Termination-Status-1).
  The MOI statuses are much richer than the MPB statuses and can be used to
  distinguish between previously indistinguishable cases (for example, did the solver
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
  [documentation](https://www.juliaopt.org/MathOptInterface.jl/v0.6.2/apimanual.html#Duals-1)).
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
  [style recommendations](https://www.juliaopt.org/JuMP.jl/dev/style/).

- We rely on broadcasting syntax to apply accessors to collections of
  variables, for example, `value.(x)` instead of `getvalue(x)` for collections. (Use
  `value(x)` when `x` is a scalar object.)

### Added

- Splatting (like `f(x...)`) is recognized in restricted settings in nonlinear
  expressions.

- Support for deleting constraints and variables.

- The documentation has been completely rewritten using docstrings and
  Documenter.

- Support for modeling mixed conic and quadratic models (for example, conic models
  with quadratic objectives and bi-linear matrix inequalities).

- Significantly improved support for modeling new types of constraints and for
  extending JuMP's macros.

- Support for providing dual warm starts.

- Improved support for accessing solver-specific attributes (for example, the
  irreducible inconsistent subsystem).

- Explicit control of whether symmetry-enforcing constraints are added to PSD
  constraints.

- Support for modeling exponential cones.

- Significant improvements in internal code quality and testing.

- Style and naming guidelines.

- Direct mode and manual mode provide explicit control over when copies of a
  model are stored or regenerated. See the corresponding
  [documentation](https://www.juliaopt.org/JuMP.jl/dev/solvers/).

### Regressions

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

- The column generation syntax in `@variable` has been removed (that is, the
  `objective`, `coefficients`, and `inconstraints` keyword arguments). Support
  for column generation will be re-introduced in a future release.

- The ability to solve the continuous relaxation (that is, via
  `solve(model; relaxation = true)`) is not yet reimplemented ([issue](https://github.com/JuliaOpt/JuMP.jl/issues/1611)).

## Version 0.18.5 (December 1, 2018)

   * Support views in some derivative evaluation functions.
   * Improved compatibility with PackageCompiler.

## Version 0.18.4 (October 8, 2018)

   * Fix a bug in model printing on Julia 0.7 and 1.0.

## Version 0.18.3 (October 1, 2018)

   * Add support for Julia v1.0 (Thanks @ExpandingMan)
   * Fix matrix expressions with quadratic functions (#1508)

## Version 0.18.2 (June 10, 2018)

   * Fix a bug in second-order derivatives when expressions are present (#1319)
   * Fix a bug in `@constraintref` (#1330)

## Version 0.18.1 (April 9, 2018)

   * Fix for nested tuple destructuring (#1193)
   * Preserve internal model when relaxation=true (#1209)
   * Minor bug fixes and updates for example

## Version 0.18.0 (July 27, 2017)

   * Drop support for Julia 0.5.
   * Update for ForwardDiff 0.5.
   * Minor bug fixes.

## Version 0.17.1 (June 9, 2017)

   * Use of `constructconstraint!` in `@SDconstraint`.
   * Minor bug fixes.

## Version 0.17.0 (May 27, 2017)

   * **Breaking change**: Mixing quadratic and conic constraints is no longer supported.
   * **Breaking change**: The `getvariable` and `getconstraint` functions are replaced by indexing on the corresponding symbol. For instance, to access the variable with name `x`, one should now write `m[:x]` instead of `getvariable(m, :x)`. As a consequence, creating a variable and constraint with the same name now triggers a warning, and accessing one of them afterwards throws an error. This change is breaking only in the latter case.
   * Addition of the `getobjectivebound` function that mirrors the functionality of the MathProgBase `getobjbound` function except that it takes into account transformations performed by JuMP.
   * Minor bug fixes.

The following changes are primarily of interest to developers of JuMP extensions:

   * The new syntax `@constraint(model, expr in Cone)` creates the constraint ensuring that `expr` is inside `Cone`. The `Cone` argument is passed to `constructconstraint!` which enables the call to the dispatched to an extension.
   * The `@variable` macro now calls `constructvariable!` instead of directly calling the `Variable` constructor. Extra arguments and keyword arguments passed to `@variable` are passed to `constructvariable!` which enables the call to be dispatched to an extension.
   * Refactor the internal function `conicdata` (used build the MathProgBase conic model) into smaller sub-functions to make these parts reusable by extensions.

## Version 0.16.2 (March 28, 2017)

   * Minor bug fixes and printing tweaks
   * Address deprecation warnings for Julia 0.6

## Version 0.16.1 (March 7, 2017)

   * Better support for `AbstractArray` in JuMP (Thanks @tkoolen)
   * Minor bug fixes

## Version 0.16.0 (February 23, 2017)

   * **Breaking change**: JuMP no longer has a mechanism for selecting solvers by default (the previous mechanism was flawed and incompatible with Julia 0.6). Not specifying a solver before calling `solve()` will result in an error.
   * **Breaking change**: User-defined functions are no longer global. The first argument to `JuMP.register` is now a JuMP `Model` object within whose scope the function will be registered. Calling `JuMP.register` without a `Model` now produces an error.
   * **Breaking change**: Use the new `JuMP.fix` method to fix a variable to a value or to update the value to which a variable is fixed. Calling `setvalue` on a fixed variable now results in an error in order to avoid silent behavior changes. (Thanks @joaquimg)
   * Nonlinear expressions now print out similarly to linear/quadratic expressions (useful for debugging!)
   * New `category` keyword to `@variable`. Used for specifying categories of anonymous variables.
   * Compatibility with Julia 0.6-dev.
   * Minor fixes and improvements (Thanks @cossio, @ccoffrin, @blegat)

## Version 0.15.1 (January 31, 2017)

  * Bugfix for `@LinearConstraints` and friends

## Version 0.15.0 (December 22, 2016)

  * Julia 0.5.0 is the minimum required version for this release.
  * Document support for BARON solver
  * Enable info callbacks in more states than before, for example, for recording solutions.
    New `when` argument to `addinfocallback` ([#814](https://github.com/JuliaOpt/JuMP.jl/pull/814), thanks @yeesian)
  * Improved support for anonymous variables. This includes new warnings for potentially confusing use of the traditional non-anonymous syntax:
    * When multiple variables in a model are given the same name
    * When non-symbols are used as names, for example, `@variable(m, x[1][1:N])`
  * Improvements in iterating over JuMP containers ([#836](https://github.com/JuliaOpt/JuMP.jl/pull/836), thanks @IssamT)
  * Support for writing variable names in .lp file output (Thanks @leethargo)
  * Support for querying duals to SDP problems (Thanks @blegat)
  * The comprehension syntax with curly braces `sum{}`, `prod{}`, and `norm2{}` has been deprecated
    in favor of Julia's native comprehension syntax `sum()`, `prod()` and `norm()` as previously announced.
    (For early adopters of the new syntax, `norm2()` was renamed to `norm()` without deprecation.)
  * Unit tests rewritten to use Base.Test instead of FactCheck
  * Improved support for operations with matrices of JuMP types (Thanks @ExpandingMan)
  * The syntax to halt a solver from inside a callback has changed from `throw(CallbackAbort())` to `return JuMP.StopTheSolver`
  * Minor bug fixes

## Version 0.14.2 (December 12, 2016)

  * Allow singleton anonymous variables (includes bugfix)

## Version 0.14.1 (September 12, 2016)

  * More consistent handling of states in informational callbacks,
    includes a new `when` parameter to `addinfocallback` for
    specifying in which state an informational callback should be called.

## Version 0.14.0 (August 7, 2016)

  * Compatibility with Julia 0.5 and ForwardDiff 0.2
  * Support for "anonymous" variables, constraints, expressions, and parameters, for example,
    `x = @variable(m, [1:N])` instead of `@variable(m, x[1:N])`
  * Support for retrieving constraints from a model by name via `getconstraint`
  * `@NLconstraint` now returns constraint references (as expected).
  * Support for vectorized expressions within lazy constraints
  * On Julia 0.5, parse new comprehension syntax `sum(x[i] for i in 1:N if isodd(i))`
    instead of `sum{ x[i], i in 1:N; isodd(i) }`. The old syntax with curly
    braces will be deprecated in JuMP 0.15.
  * Now possible to provide nonlinear expressions as "raw" Julia `Expr` objects
    instead of using JuMP's nonlinear macros. This input format is useful for
    programmatically generated expressions.
  * `s/Mathematical Programming/Mathematical Optimization/`
  * Support for local cuts (Thanks to @madanim, Mehdi Madani)
  * Document Xpress interface developed by @joaquimg, Joaquim Dias Garcia
  * Minor bug and deprecation fixes (Thanks @odow, @jrevels)

## Version 0.13.2 (May 16, 2016)

  * Compatibility update for MathProgBase

## Version 0.13.1 (May 3, 2016)

  * Fix broken deprecation for `registerNLfunction`.

## Version 0.13.0 (April 29, 2016)

  * Most exported methods and macros have been renamed to avoid camelCase. See the list of changes [here](https://github.com/JuliaOpt/JuMP.jl/blob/e53d0db67cde2a4b80d0c1281f4b49eb0128a1f5/src/deprecated.jl#L30). There is a 1-1 mapping from the old names to the new, and it is safe to simply replace the names to update existing models.
  * Specify variable lower/upper bounds in `@variable` using the `lowerbound` and `upperbound` keyword arguments.
  * Change name printed for variable using the `basename` keyword argument to `@variable`.
  * New `@variables` macro allows multi-line declaration of groups of variables.
  * A number of solver methods previously available only through MathProgBase are now exposed directly in JuMP. The fix was [recorded](https://youtu.be/qF1lZPJ3a5A) live.
  * Compatibility fixes with Julia 0.5.
  * The "end" indexing syntax is no longer supported within JuMPArrays which do not use 1-based indexing until upstream issues are resolved, see [here](https://github.com/JuliaOpt/JuMP.jl/issues/730).

## Version 0.12.2 (March 9, 2016)

  * Small fixes for nonlinear optimization

## Version 0.12.1 (March 1, 2016)

  * Fix a regression in slicing for JuMPArrays (when not using 1-based indexing)

## Version 0.12.0 (February 27, 2016)

  * The automatic differentiation functionality has been completely rewritten with a number of user-facing changes:
      - `@defExpr` and `@defNLExpr` now take the model as the first argument. The previous one-argument version of `@defExpr` is deprecated; all expressions should be named. For example, replace `@defExpr(2x+y)` with `@defExpr(jump_model, my_expr, 2x+y)`.
      - JuMP no longer uses Julia's variable binding rules for efficiently re-solving a sequence of nonlinear models. Instead, we have introduced nonlinear parameters. This is a breaking change, so we have added a warning message when we detect models that may depend on the old behavior.
      - Support for user-defined functions integrated within nonlinear JuMP expressions.
  * Replaced iteration over `AffExpr` with `Number`-like scalar iteration; previous iteration behavior is now available via `linearterms(::AffExpr)`.
  * Stopping the solver via `throw(CallbackAbort())` from a callback no longer triggers an exception. Instead, `solve()` returns `UserLimit` status.
  * `getDual()` now works for conic problems (Thanks @emreyamangil.)

## Version 0.11.3 (February 4, 2016)

  * Bug-fix for problems with quadratic objectives and semidefinite constraints

## Version 0.11.2 (January 14, 2016)

  * Compatibility update for Mosek

## Version 0.11.1 (December 1, 2015)

  * Remove usage of `@compat` in tests.
  * Fix updating quadratic objectives for nonlinear models.

## Version 0.11.0 (November 30, 2015)

  * Julia 0.4.0 is the minimum required version for this release.
  * Fix for scoping semantics of index variables in sum{}. Index variables no longer leak into the surrounding scope.
  * Addition of the `solve(m::Model, relaxation=true)` keyword argument to solve the standard continuous relaxation of model `m`
  * The `getConstraintBounds()` method allows access to the lower and upper bounds of all constraints in a (nonlinear) model.
  * Update for breaking changes in MathProgBase

## Version 0.10.3 (November 20, 2015)

  * Fix a rare error when parsing quadratic expressions
  * Fix `Variable()` constructor with default arguments
  * Detect unrecognized keywords in `solve()`

## Version 0.10.2 (September 28, 2015)

  * Fix for deprecation warnings

## Version 0.10.1 (September 3, 2015)

  * Fixes for ambiguity warnings.
  * Fix for breaking change in precompilation syntax in Julia 0.4-pre

## Version 0.10.0 (August 31, 2015)

  * Support (on Julia 0.4 and later) for conditions in indexing `@defVar` and `@addConstraint` constructs, for example, `@defVar(m, x[i=1:5,j=1:5; i+j >= 3])`
  * Support for vectorized operations on Variables and expressions. See the documentation for details.
  * New `getVar()` method to access variables in a model by name
  * Support for semidefinite programming.
  * Dual solutions are now available for general nonlinear problems. You may call `getDual` on a reference object for a nonlinear constraint, and `getDual` on a variable object for Lagrange multipliers from active bounds.
  * Introduce warnings for two common performance traps: too many calls to `getValue()` on a collection of variables and use of the `+` operator in a loop to sum expressions.
  * Second-order cone constraints can be written directly with the `norm()` and `norm2{}` syntax.
  * Implement MathProgBase interface for querying Hessian-vector products.
  * Iteration over `JuMPContainer`s is deprecated; instead, use the `keys` and `values` functions, and `zip(keys(d),values(d))` for the old behavior.
  * `@defVar` returns `Array{Variable,N}` when each of `N` index sets are of the form `1:nᵢ`.
  * Module precompilation: on Julia 0.4 and later, `using JuMP` is now much faster.

## Version 0.9.3 (August 11, 2015)

  * Fixes for FactCheck testing on julia v0.4.

## Version 0.9.2 (June 27, 2015)

  * Fix bug in @addConstraints.

## Version 0.9.1 (April 25, 2015)

  * Fix for Julia 0.4-dev.
  * Small infrastructure improvements for extensions.

## Version 0.9.0 (April 18, 2015)

  * Comparison operators for constructing constraints (for example, `2x >= 1`) have been deprecated. Instead, construct the constraints explicitly in
    the `@addConstraint` macro to add them to the model, or in the `@LinearConstraint` macro to create a stand-alone linear constraint instance.
  * `getValue()` method implemented to compute the value of a nonlinear subexpression
  * JuMP is now released under the Mozilla Public License version 2.0 (was previously LGPL). MPL is a copyleft license which is less restrictive than LGPL, especially for embedding JuMP within other applications.
  * A number of performance improvements in ReverseDiffSparse for computing derivatives.
  * `MathProgBase.getsolvetime(m)` now returns the solution time reported by the solver, if available. (Thanks @odow, Oscar Dowson)
  * Formatting fix for LP format output. (Thanks @sbebo, Leonardo Taccari).

## Version 0.8.0 (February 17, 2015)

  * Nonlinear subexpressions now supported with the `@defNLExpr` macro.
  * SCS supported for solving second-order conic problems.
  * `setXXXCallback` family deprecated in favor of `addXXXCallback`.
  * Multiple callbacks of the same type can be registered.
  * Added support for informational callbacks via `addInfoCallback`.
  * A `CallbackAbort` exception can be thrown from callback to safely exit optimization.

## Version 0.7.4 (February 4, 2015)

  * Reduced costs and linear constraint duals are now accessible when quadratic constraints are present.
  * Two-sided nonlinear constraints are supported.
  * Methods for accessing the number of variables and constraints in a model are renamed.
  * New default procedure for setting initial values in nonlinear optimization: project zero onto the variable bounds.
  * Small bug fixes.

## Version 0.7.3 (January 14, 2015)

  * Fix a method ambiguity conflict with Compose.jl (cosmetic fix)

## Version 0.7.2 (January 9, 2015)

  * Fix a bug in `sum(::JuMPDict)`
  * Added the `setCategory` function to change a variables category (for example, continuous or binary)
  after construction, and `getCategory` to retrieve the variable category.

## Version 0.7.1 (January 2, 2015)

  * Fix a bug in parsing linear expressions in macros. Affects only Julia 0.4 and later.

## Version 0.7.0 (December 29, 2014)

### Linear/quadratic/conic programming

  * **Breaking change**: The syntax for column-wise model generation has been changed to use keyword arguments in `@defVar`.
  * On Julia 0.4 and later, variables and coefficients may be multiplied in any order within macros. That is, variable*coefficient is now valid syntax.
  * ECOS supported for solving second-order conic problems.

### [Nonlinear programming](@id _nonlinear_programming_release_notes)

  * Support for skipping model generation when solving a sequence of nonlinear models with changing data.
  * Fix a memory leak when solving a sequence of nonlinear models.
  * The `@addNLConstraint` macro now supports the three-argument version to define sets of nonlinear constraints.
  * KNITRO supported as a nonlinear solver.
  * Speed improvements for model generation.
  * The `@addNLConstraints` macro supports adding multiple (groups of) constraints at once. Syntax is similar to `@addConstraints`.
  * Discrete variables allowed in nonlinear problems for solvers which support them (currently only KNITRO).

### General

  * Starting values for variables may now be specified with `@defVar(m, x, start=value)`.
  * The `setSolver` function allows users to change the solver subsequent to model creation.
  * Support for "fixed" variables via the `@defVar(m, x == 1)` syntax.
  * Unit tests rewritten to use FactCheck.jl, improved testing across solvers.

## Version 0.6.3 (October 19, 2014)

  * Fix a bug in multiplying two AffExpr objects.

## Version 0.6.2 (October 11, 2014)

  * Further improvements and bug fixes for printing.
  * Fixed a bug in `@defExpr`.
  * Support for accessing expression graphs through the MathProgBase NLP interface.

## Version 0.6.1 (September 19, 2014)

  * Improvements and bug fixes for printing.

## Version 0.6.0 (September 9, 2014)

  * Julia 0.3.0 is the minimum required version for this release.
  * `buildInternalModel(m::Model)` added to build solver-level model in memory without optimizing.
  * Deprecate `load_model_only` keyword argument to `solve`.
  * Add groups of constraints with `@addConstraints` macro.
  * Unicode operators now supported, including `∑` for `sum`, `∏` for `prod`, and `≤`/`≥`
  * Quadratic constraints supported in `@addConstraint` macro.
  * Quadratic objectives supported in `@setObjective` macro.
  * MathProgBase solver-independent interface replaces Ipopt-specific interface for nonlinear problems
    - **Breaking change**: `IpoptOptions` no longer supported to specify solver options, use `m = Model(solver=IpoptSolver(options...))` instead.
  * New solver interfaces: ECOS, NLopt, and nonlinear support for MOSEK
  * New option to control whether the lazy constraint callback is executed at each node in the B&B tree or just when feasible solutions are found
  * Add support for semicontinuous and semi-integer variables for those solvers that support them.
  * Add support for index dependencies (for example, triangular indexing) in `@defVar`, `@addConstraint`, and `@defExpr` (for example, `@defVar(m, x[i=1:10,j=i:10])`).
    - This required some changes to the internal structure of JuMP containers, which may break code that explicitly stored `JuMPDict` objects.

## Version 0.5.8 (September 24, 2014)

  * Fix a bug with specifying solvers (affects Julia 0.2 only)

## Version 0.5.7 (September 5, 2014)

  * Fix a bug in printing models

## Version 0.5.6 (September 2, 2014)

  * Add support for semicontinuous and semi-integer variables for those solvers that support them.
    - **Breaking change**: Syntax for `Variable()` constructor has changed (use of this interface remains discouraged)
  * Update for breaking changes in MathProgBase

## Version 0.5.5 (July 6, 2014)

  * Fix bug with problem modification: adding variables that did not appear in existing constraints or objective.

## Version 0.5.4 (June 19, 2014)

  * Update for breaking change in MathProgBase which reduces loading times for `using JuMP`
  * Fix error when MIPs not solved to optimality

## Version 0.5.3 (May 21, 2014)

  * Update for breaking change in ReverseDiffSparse

## Version 0.5.2 (May 9, 2014)

  * Fix compatibility with Julia 0.3 prerelease

## Version 0.5.1 (May 5, 2014)

  * Fix a bug in coefficient handling inside lazy constraints and user cuts

## Version 0.5.0 (May 2, 2014)

  * Support for nonlinear optimization with exact, sparse second-order derivatives automatically computed. Ipopt is currently the only solver supported.
  * `getValue` for `AffExpr` and `QuadExpr`
  * **Breaking change**: `getSolverModel` replaced by `getInternalModel`, which returns the internal MathProgBase-level model
  * Groups of constraints can be specified with `@addConstraint` (see documentation for details). This is not a breaking change.
  * `dot(::JuMPDict{Variable},::JuMPDict{Variable})` now returns the corresponding quadratic expression.

## Version 0.4.1 (March 24, 2014)

  * Fix bug where change in objective sense was ignored when re-solving a model.
  * Fix issue with handling zero coefficients in AffExpr.

## Version 0.4.0 (March 10, 2014)

  * Support for SOS1 and SOS2 constraints.
  * Solver-independent callback for user heuristics.
  * `dot` and `sum` implemented for `JuMPDict` objects. Now you can say `@addConstraint(m, dot(a,x) <= b)`.
  * Developers: support for extensions to JuMP. See definition of Model in `src/JuMP.jl` for more details.
  * Option to construct the low-level model before optimizing.

## Version 0.3.2 (February 17, 2014)

 * Improved model printing
   - Preliminary support for IJulia output

## Version 0.3.1 (January 30, 2014)

 * Documentation updates
 * Support for MOSEK
 * CPLEXLink renamed to CPLEX

## Version 0.3.0 (January 21, 2014)

 * Unbounded/infeasibility rays: `getValue()` will return the corresponding
   components of an unbounded ray when a model is unbounded, if supported
   by the selected solver. `getDual()` will return an infeasibility ray (Farkas
   proof) if a model is infeasible and the selected solver supports this
   feature.
 * Solver-independent callbacks for user generated cuts.
 * Use new interface for solver-independent QCQP.
 * `setlazycallback` renamed to `setLazyCallback` for consistency.

## Version 0.2.0 (December 15, 2013)

### Breaking

   * Objective sense is specified in `setObjective` instead of in the `Model`
     constructor.
   * `lpsolver` and `mipsolver` merged into single `solver` option.

### Added

   * Problem modification with efficient LP restarts and MIP warm-starts.
   * Relatedly, column-wise modeling now supported.
   * Solver-independent callbacks supported. Currently we support only
     a "lazy constraint" callback, which works with Gurobi, CPLEX, and GLPK.
     More callbacks coming soon.

## Version 0.1.2 (November 16, 2013)

  * Bug fixes for printing, improved error messages.
  * Allow `AffExpr` to be used in macros; for example,
    `ex = y + z; @addConstraint(m, x + 2*ex <= 3)`

## Version 0.1.1 (October 23, 2013)

  * Update for solver specification API changes in MathProgBase.

## Version 0.1.0 (October 3, 2013)

  * Initial public release.

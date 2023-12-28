# Development roadmap

The JuMP developers have compiled this roadmap document to share their plans and
goals with the JuMP community. Contributions to roadmap issues are especially
invited.

Most of these issues will require changes to both JuMP and MathOptInterface, and
are non-trivial in their implementation. They are in no particular order, but
represent broad themes that we see as areas in which JuMP could be improved.

 - Support nonlinear expressions with vector-valued inputs and outputs. There
   are a few related components:
    - Representing terms like `log(det(X))` as necessary for `Convex.jl`
    - Automatic differentiation of terms with vector inputs and outputs
    - User-defined functions with vector--as opposed to scalar--inputs, which is
      particularly useful for optimal control problems
    - User-defined functions with vector outputs, avoiding the need for
      [User-defined operators with vector outputs](@ref)
 - Add support for modeling with SI units. The [UnitJuMP.jl](https://github.com/trulsf/UnitJuMP.jl)
   extension is a good proof of concept for what this would look like. We want
   to make units a first-class concept in JuMP. See [#1350](https://github.com/jump-dev/JuMP.jl/issues/1350)
   for more details.

## Completed

 - **Done [#3106](https://github.com/jump-dev/JuMP.jl/pull/3106)** Make
   nonlinear programming a first-class citizen. There have been many issues
   and discussions about this: currently nonlinear constraints are handled
   through a `MOI.NLPBlock` and have various limitations and restrictions.
   - [https://github.com/jump-dev/JuMP.jl/issues/1185](https://github.com/jump-dev/JuMP.jl/issues/1185)
   - [https://github.com/jump-dev/JuMP.jl/issues/1198](https://github.com/jump-dev/JuMP.jl/issues/1198)
   - [https://github.com/jump-dev/JuMP.jl/issues/2788](https://github.com/jump-dev/JuMP.jl/issues/2788)
   - [https://github.com/jump-dev/MathOptInterface.jl/issues/846](https://github.com/jump-dev/MathOptInterface.jl/issues/846)
   - [https://github.com/jump-dev/MathOptInterface.jl/issues/1397](https://github.com/jump-dev/MathOptInterface.jl/issues/1397)
 - **Done [#3385](https://github.com/jump-dev/JuMP.jl/pull/3385)** Add support for coefficient types other than `Float64`:
   [https://github.com/jump-dev/JuMP.jl/issues/2025](https://github.com/jump-dev/JuMP.jl/issues/2025)
   Since the very beginning, JuMP has hard-coded the coefficient type as
   `Float64`. This has made it impossible to support solvers which can use other
   types such as `BigFloat` or `Rational{BigInt}`.
 - **Done [#3385](https://github.com/jump-dev/JuMP.jl/pull/3635)** Add support for constraint programming:
   [https://github.com/jump-dev/JuMP.jl/issues/2227](https://github.com/jump-dev/JuMP.jl/issues/2227)
   JuMP has a strong focus on linear, conic and nonlinear optimization problems.
   We want to add better support for constraint programming.
 - **Done [#3176](https://github.com/jump-dev/JuMP.jl/pull/3176)** Add support for multiobjective problems:
   [https://github.com/jump-dev/JuMP.jl/issues/2099](https://github.com/jump-dev/JuMP.jl/issues/2099)
   JuMP is restricted to problems with scalar-valued objectives. We want to
   extend this to vector-valued problems.
 - **Done [#3629](https://github.com/jump-dev/JuMP.jl/pull/3629)** Refactor the
   internal code of JuMP's macros. The code in `src/macros.jl` is
   some of the oldest part of JuMP and is difficult to read, modify, and extend.
   We should overhaul the internals of JuMP's macros---without making
   user-visible breaking changes---to improve their long-term maintainability.

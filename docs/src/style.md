Style guide and design principles
=================================

Style guide
-----------

### JuMP Source Code

 - 80 character line limit. (Except for long URLs in comments.)
 - 4 spaces for indentation.
- Use `lowercase_with_underscore` for variable and function names.
- Types and Modules should be in `CamelCase`.
- Use `UPPERCASE` for constants.
 - No one-letter variable names. This includes using `m = Model()`. The only
exceptions are for indices in loops such as `for i in 1:N`.
 - Start internal functions (i.e., those not meant to be used by a user) with an
 underscore, e.g., `_add_constraint`.
 - Use three backticks for name of internal function docstrings.
 - Variables should never start with an underscore.
 - To instantiate a vector, use `T[]` instead of `Vector{T}()`.
 - Use proper punctuation and full sentences in comments.
 - Use `for i in 1:N` instead of `for i = 1:N`.
 - Never use trailing periods for floating point, e.g., `0.`.

### Writing a JuMP Model

These style points address the convention for how users should write JuMP
models.

 - Never use `m` for the name of a model.
 - If the model fits on one screen, short variable names are okay (e.g.,
`@variable(model, x)`).
 - If the model is to large to fit on one screen, its construction should be
 decomposed into functions with long variable names.

#### Inside JuMP model macros

 - Never use `1.0` when `1` is okay.
 - Use the order `coefficient * variable` instead of `variable * coefficient`.
 - Use parenthesis in macros
 - Prefer blocks of `@variables` over lots of `@variable` lines. (The same goes
for the other plural versions of the macros.)


Design principles
-----------------

TODO: How to structure and test large JuMP models, libraries that use JuMP.

For how to write a solver, see MOI.

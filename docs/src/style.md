Style guide and design principles
=================================

This document outlines JuMP's style guide. It also contains a set or principles
for users writing JuMP models.

Style guide
-----------

!!! info
    For developers: if you modify JuMP's source code, please fix the style of
    the surrounding code while you are there.

### JuMP Source Code

 - 80 character line limit. (Except for long URLs in comments.)
 - 4 spaces for indentation.
 - Use `lowercase_with_underscore` for variable and function names.
 - Types and Modules should be in `CamelCase`.
 - Use `UPPERCASE` for constants.
 - No bangs (!) in function names. Use docstrings to indicate what is modified
and variable lifetime issues.
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
 - Never more than one consecutive blank line inside a function.
 - Always use return from a function unless it is a one-line function.
 - Do not split a one-line function over lines. Only if less than 80 characters.
 - Type arguments to user-facing functions so that users do not get internal
errors when passing wrong types.
 - Imports: use `import` instead of `using`. Importing specific things from a
 module is fine, stdlib imports first, alphabetical order after that.

### Doc strings

 - ` ```julia ` not ` ``` ` or 4 spaces indent for pretty printing
 - Every function needs a doc string except if it is small and obvious
 - 4 spaces, signature, empty line and then the doc
 - Always use complete English sentences with proper punctuation
 - Do not terminate lists with punctuation (e.g., as in this style guide)

### Testing guidelines

 - Write tests for all PR’s
 - And documentation for all new features
 - Omit "test" or "tests" in the name of the testset
 - Test one thing per bottom-level testset

Design principles
-----------------

We now outline some principles that we think users should adopt when using JuMP.

TODO: How to structure and test large JuMP models, libraries that use JuMP.
For how to write a solver, see MOI.

### Writing a JuMP Model

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

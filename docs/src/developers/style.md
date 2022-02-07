# Style guide and design principles

## Style guide

This section describes the coding style rules that apply to JuMP code and that
we recommend for JuMP models and surrounding Julia code. The motivations for
a style guide include:

- conveying best practices for writing readable and maintainable code
- reducing the amount of time spent on
  [bike-shedding](https://en.wikipedia.org/wiki/Law_of_triviality) by
  establishing basic naming and formatting conventions
- lowering the barrier for new contributors by codifying the existing practices
  (e.g., you can be more confident your code will pass review if you follow the style guide)

In some cases, the JuMP style guide diverges from the
[Julia style guide](https://docs.julialang.org/en/v1.0.0/manual/style-guide/).
All such cases will be explicitly noted and justified.

The JuMP style guide adopts many recommendations from the
[Google style guides](https://github.com/google/styleguide).

!!! info
    The style guide is always a work in progress, and not all JuMP code
    follows the rules. When modifying JuMP, please fix the style violations
    of the surrounding code (i.e., leave the code tidier than when you
    started). If large changes are needed, consider separating them into
    another PR.

### JuliaFormatter

JuMP uses [JuliaFormatter.jl](https://github.com/domluna/JuliaFormatter.jl) as
an autoformatting tool.

We use the options contained in [`.JuliaFormatter.toml`](https://github.com/jump-dev/JuMP.jl/blob/master/.JuliaFormatter.toml).

To format code, `cd` to the JuMP directory, then run:
```julia
] add JuliaFormatter@0.22.2
using JuliaFormatter
format("src")
format("test")
```

!!! info
    A continuous integration check verifies that all PRs made to JuMP have
    passed the formatter.

The following sections outline extra style guide points that are not fixed
automatically by JuliaFormatter.

### Whitespace

For conciseness, never use more than one blank line within a function, and never
begin a function with a blank line.

Bad:
```julia
function foo(x)
    y = 2 * x


    return y
end

function foo(x)

    y = 2 * x
    return y
end
```

### Juxtaposed multiplication

Only use juxtaposed multiplication when the right-hand side is a symbol.

Good:
```julia
2x  # Acceptable if there are space constraints.
2 * x  # This is preferred if space is not an issue.
2 * (x + 1)
```

Bad:
```julia
2(x + 1)
```

### Empty vectors

For a type `T`, `T[]` and `Vector{T}()` are equivalent ways to create an
empty vector with element type `T`. Prefer `T[]` because it is more concise.

### Comments

For non-native speakers and for general clarity, comments in code must be proper
English sentences with appropriate punctuation.

Good:
```julia
# This is a comment demonstrating a good comment.
```

Bad:
```julia
# a bad comment
```

### JuMP macro syntax

For consistency, always use parentheses.

Good:
```julia
@variable(model, x >= 0)
```

Bad:
```julia
@variable model x >= 0
```

For consistency, always use `constant * variable` as opposed to
`variable * constant`. This makes it easier to read models in
ambiguous cases like `a * x`.

Good:
```julia
a = 4
@constraint(model, 3 * x <= 1)
@constraint(model, a * x <= 1)
```

Bad:
```julia
a = 4
@constraint(model, x * 3 <= 1)
@constraint(model, x * a <= 1)
```

In order to reduce boilerplate code, prefer the plural form of macros over lots
of repeated calls to singular forms.

Good:
```julia
@variables(model, begin
    x >= 0
    y >= 1
    z <= 2
end)
```

Bad:
```julia
@variable(model, x >= 0)
@variable(model, y >= 1)
@variable(model, z <= 2)
```

An exception is made for calls with many keyword arguments, since these need to
be enclosed in parentheses in order to parse properly.

Acceptable:
```julia
@variable(model, x >= 0, start = 0.0, base_name = "my_x")
@variable(model, y >= 1, start = 2.0)
@variable(model, z <= 2, start = -1.0)
```

Also acceptable:
```julia
@variables(model, begin
    x >= 0, (start = 0.0, base_name = "my_x")
    y >= 1, (start = 2.0)
    z <= 2, (start = -1.0)
end)
```

While we always use `in` for `for`-loops, it is acceptable to use `=` in the
container declarations of JuMP macros.

Okay:
```julia
@variable(model, x[i=1:3])
```
Also okay:
```julia
@variable(model, x[i in 1:3])
```

### Naming

```julia
module SomeModule end
function some_function end
const SOME_CONSTANT = ...
struct SomeStruct
  some_field::SomeType
end
@enum SomeEnum ENUM_VALUE_A ENUM_VALUE_B
some_local_variable = ...
some_file.jl # Except for ModuleName.jl.
```

### Exported and non-exported names

Begin private module level functions and constants with an underscore. All other
objects in the scope of a module should be exported. (See JuMP.jl for an example
of how to do this.)

Names beginning with an underscore should only be used for distinguishing
between exported (public) and non-exported (private) objects. Therefore, never
begin the name of a local variable with an underscore.

```julia
module MyModule

export public_function, PUBLIC_CONSTANT

function _private_function()
    local_variable = 1
    return
end

function public_function end

const _PRIVATE_CONSTANT = 3.14159
const PUBLIC_CONSTANT = 1.41421

end
```

### Use of underscores within names

The Julia style guide recommends avoiding underscores "when readable", for
example, `haskey`, `isequal`, `remotecall`, and `remotecall_fetch`. This
convention creates the potential for unnecessary bikeshedding and also forces
the user to recall the presence/absence of an underscore, e.g., "was that
argument named `basename` or `base_name`?". For consistency, *always use
underscores* in variable names and function names to separate words.

### Use of `!`

Julia has a convention of appending `!` to a function name if the function
modifies its arguments. We recommend to:

- Omit `!` when the name itself makes it clear that modification is taking
  place, e.g., `add_constraint` and `set_name`. We depart from the Julia style
  guide because `!` does not provide a reader with any additional information
  in this case, and adherence to this convention is not uniform even in base
  Julia itself (consider `Base.println` and `Base.finalize`).
- Use `!` in all other cases. In particular it can be used to distinguish
  between modifying and non-modifying variants of the same function like `scale`
  and `scale!`.

Note that `!` is *not* a self-documenting feature because it is still
ambiguous which arguments are modified when multiple arguments are present.
Be sure to document which arguments are modified in the method's docstring.

See also the Julia style guide recommendations for
[ordering of function arguments](https://docs.julialang.org/en/v1/manual/style-guide/#Write-functions-with-argument-ordering-similar-to-Julia-Base).

### Abbreviations

Abbreviate names to make the code more readable, not to save typing.
Don't arbitrarily delete letters from a word to abbreviate it (e.g., `indx`).
Use abbreviations consistently within a body of code (e.g., do not mix
`con` and `constr`, `idx` and `indx`).

Common abbreviations:

- `num` for `number`
- `con` for `constraint`

### No one-letter variable names

Where possible, avoid one-letter variable names.

Use `model = Model()` instead of `m = Model()`

Exceptions are made for indices in loops.

### User-facing `MethodError`

Specifying argument types for methods is mostly optional in Julia, which means
that it's possible to find out that you are working with unexpected types deep in
the call chain. Avoid this situation or handle it with a helpful error message.
*A user should see a `MethodError` only for methods that they called directly.*

Bad:
```julia
_internal_function(x::Integer) = x + 1
# The user sees a MethodError for _internal_function when calling
# public_function("a string"). This is not very helpful.
public_function(x) = _internal_function(x)
```

Good:
```julia
_internal_function(x::Integer) = x + 1
# The user sees a MethodError for public_function when calling
# public_function("a string"). This is easy to understand.
public_function(x::Integer) = _internal_function(x)
```

If it is hard to provide an error message at the top of the call chain,
then the following pattern is also ok:
```julia
_internal_function(x::Integer) = x + 1
function _internal_function(x)
    error(
        "Internal error. This probably means that you called " *
        "public_function() with the wrong type.",
    )
end
public_function(x) = _internal_function(x)
```

### `@enum` vs. `Symbol`

The `@enum` macro lets you define types with a finite number of values that
are explicitly enumerated (like `enum` in C/C++). `Symbol`s are lightweight
strings that are used to represent identifiers in Julia (for example, `:x`).

`@enum` provides type safety and can have docstrings attached to explain the
possible values. Use `@enum`s when applicable, e.g., for reporting statuses.
Use strings to provide long-form additional information like error messages.

Use of `Symbol` should typically be reserved for identifiers, e.g., for lookup
in the JuMP model (`model[:my_variable]`).

### `using` vs. `import`

`using ModuleName` brings all symbols exported by the module `ModuleName`
into scope, while `import ModuleName` brings only the module itself into scope.
(See the Julia
[manual](https://docs.julialang.org/en/v1/manual/modules/#modules-1)) for
examples and more details.

For the same reason that `from <module> import *` is not recommended in python
([PEP 8](https://www.python.org/dev/peps/pep-0008/#imports)), avoid
`using ModuleName` except in throw-away scripts or at the REPL. The `using`
statement makes it harder to track where symbols come from and exposes the code
to ambiguities when two modules export the same symbol.

Prefer `using ModuleName: x, p` to `import ModuleName.x, ModuleName.p` and
`import MyModule: x, p` because the `import` versions allow method extension
without qualifying with the module name.

Similarly, `using ModuleName: ModuleName` is an acceptable substitute for
`import ModuleName`, because it does not bring all symbols exported by
`ModuleName` into scope. However, we prefer `import ModuleName` for consistency.

## Documentation

This section describes the writing style that should be used when writing
documentation for JuMP (and supporting packages).

We can recommend the documentation style guides by [Divio](https://www.divio.com/blog/documentation/),
[Google](https://developers.google.com/style/), and [Write the Docs](https://www.writethedocs.org/guide/)
as general reading for those writing documentation. This guide delegates a
thorough handling of the topic to those guides and instead elaborates on the
points more specific to Julia and documentation that use [Documenter](https://github.com/JuliaDocs/Documenter.jl).

 - Be concise
 - Use lists instead of long sentences
 - Use numbered lists when describing a sequence, e.g., (1) do X, (2) then Y
 - Use bullet points when the items are not ordered
 - Example code should be covered by doctests
 - When a word is a Julia symbol and not an English word, enclose it with
   backticks. In addition, if it has a docstring in this doc add a link using
   `@ref`. If it is a plural, add the "s" after the closing backtick. For example,
   ```
   [`VariableRef`](@ref)s
   ```
 - Use [`@meta`](https://juliadocs.github.io/Documenter.jl/v0.21/man/syntax/#@meta-block-1)
   blocks for TODOs and other comments that shouldn't be visible to readers.
   For example,
   ````@markdown
   ```@meta
   # TODO: Mention also X, Y, and Z.
   ```
   ````
### Docstrings

- Every exported object needs a docstring
- All examples in docstrings should be [`jldoctests`](https://juliadocs.github.io/Documenter.jl/stable/man/doctests/)
- Always use complete English sentences with proper punctuation
- Do not terminate lists with punctuation (e.g., as in this doc)

Here is an example:
````julia
"""
    signature(args; kwargs...)

Short sentence describing the function.

Optional: add a slightly longer paragraph describing the function.

## Notes

 - List any notes that the user should be aware of

## Examples

```jldoctest
julia> 1 + 1
2
```
"""
````

## Testing

Use a module to encapsulate tests, and structure all tests as functions. This
avoids leaking local variables between tests.

Here is a basic skeleton:
```julia
module TestPkg

using Test

_helper_function() = 2

function test_addition()
    @test 1 + 1 == _helper_function()
end

function runtests()
    for name in names(@__MODULE__; all = true)
        if startswith("$(name)", "test_")
            @testset "$(name)" begin
                getfield(@__MODULE__, name)()
            end
        end
    end
end

end # TestPkg

TestPkg.runtests()
```

Break the tests into multiple files, with one module per file, so that subsets
of the codebase can be tested by calling `include` with the relevant file.

## Design principles

TODO: How to structure and test large JuMP models, libraries that use JuMP.

For how to write a solver, see MOI.

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

### Formatting

Julia unfortunately does not have an autoformatting tool like
[gofmt](https://blog.golang.org/go-fmt-your-code). Until a reliable
autoformatting tool is available, we adopt the following conventions.

#### Whitespace

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

Julia is mostly insensitive to whitespace characters within lines.
For consistency:

- Use spaces between binary operators (with some exceptions, see below)
- Use a single space after commas and semicolons
- Do not use extra spaces for unary operators, parentheses, or braces
- Indent within new blocks (except `module`) using 4 spaces

Good:
```julia
f(x, y) = [3 * dot(x, y); x']
```

Bad:
```julia
f(x,y) = [ 3*dot(x,y) ; x' ]
```

Good:
```julia
module Foo

function f(x)
    return x + 1
end

end # module Foo
```

##### Exceptions

For aesthetic reasons, we make an exception for whitespace surrounding the
exponential operator `^`.

Good:
```julia
f(x) = x^2
```

Bad:
```julia
f(x) = x ^ 2
```

We also make an exception for the `:` operator when it is used to form a range.

Good:
```julia
x = 1:5
```

Bad:
```julia
x = 1 : 5
```

One reason is that it can be confused with Julia's conditional statement:
`cond ? x : y` which requires whitespace around the `:`.

We also make an exception for juxtaposed multiplication (i.e. dropping the `*`
between a numeric literal and an expression) when the right-hand side is a
symbol.

Good:
```julia
2x  # Acceptable if there are space constraints.
2 * x  # This preferred if space is not an issue.
2 * (x + 1)
```

Bad:
```julia
2(x + 1)
```

#### Return statements

To avoid situations in which it is unclear whether the author intended to return
a certain value or not, always use an explicit `return` statement to exit from a
function. If the return from a function is `nothing`, use `return` instead of
`return nothing`.

We make an exception for assignment-form one-line functions (`f(x) = 2x`).

Good:
```julia
foo(x) = 2x  # Acceptable if one line
function foo(x)
    return 2x
end
function foo(x)
    x[1] += 1
    return
end
```

Bad:
```julia
function foo(x)
    2x
end
function foo(x)
    x[1] += 1
    return nothing
end
```

#### Line length

Line lengths are a contentious issue. Our foremost goal is to maximize code
readability. Very long line lengths can be hard to easily comprehend. However,
arbitrarily enforcing a maximum line length (like 80 characters) inevitably
leads to cases in which slightly longer lines (e.g. 81 characters) might be more
readable.

Therefore, aim to keep line lengths under 80 characters by breaking lines
for maximum readability (examples are given in the [Line breaks](@ref) section),
but don't treat this as a hard rule.

We make exceptions for
 - URLs
 - pathnames

#### Line breaks

The "readability" of a line is subjective. In this section we give suggestions
of good and bad style of how to break a line. These suggestions are inspired by
Google's [Python style guide](https://google.github.io/styleguide/pyguide.html).

!!! note
    If you're unsure about how format your code, you can experiment (in Python)
    using [YAPF](https://yapf.now.sh/).

When defining functions, align arguments vertically after the opening
parenthesis, or list all arguments on a new (indented) line.

Good:
```julia
# Arguments to the function are aligned vertically.
function my_very_long_function_name(with_lots_of_long_arguments_1,
                                    and_another_long_one)
    # First line of the function begins here.
end

# Arguments to the function are listed on a new line and indented.
function my_very_long_function_name(
    with_lots_of_long_arguments_1, and_another_long_one)
    # First line of the function begins here.
end
```

Bad:
```julia
# When defining functions, if vertical alignment is not used, then the arguments
# should not begin on the first line.
function my_very_long_function_name(with_lots_of_long_arguments_1,
    and_another_long_one)
    # First line of the function begins here.
end
```

Don't use vertical alignment if all of the arguments are very far to the right.

Bad:
```julia
a_very_long_variable_name = a_long_variable_name_with_arguments(first_argument,
                                                                second_argument)
```

Better:
```julia
a_very_long_variable_name = a_long_variable_name_with_arguments(
    first_argument, second_argument)
```

Don't use vertical alignment if it would be more readable to place all arguments
on a new indented line.

Bad:
```julia
con_index = MOI.add_constraint(backend(owner_model(variable)),
                               MOI.SingleVariable(index(variable)), set)
```

Better:
```julia
con_index = MOI.add_constraint(
    backend(owner_model(variable)), MOI.SingleVariable(index(variable)), set)
```

Don't break lines at an inner-level of function nesting.

Bad:
```julia
con_index = MOI.add_constraint(
    backend(owner_model(variable)), MOI.SingleVariable(
    index(variable)), new_set)
```

Better:
```julia
con_index = MOI.add_constraint(
    backend(owner_model(variable)),
    MOI.SingleVariable(index(variable)), new_set)
```

For readability, don't split a one-line function over multiple lines.

Bad:
```julia
f(x) = 1 + x +
    x^2 + x^3
```

Better:
```julia
f(x) = 1 + x + x^2 + x^3 + x^3
```

### Syntax

Julia sometimes provides equivalent syntax to express the same basic
operation. We discuss these cases below.

#### `for` loops

Julia allows both `for x = 1:N` and `for x in 1:N`. Always prefer to use
`in` over `=`, because `in` generalizes better to other index sets like `for x in eachindex(A)`.

#### Empty vectors

For a type `T`, `T[]` and `Vector{T}()` are equivalent ways to create an
empty vector with element type `T`. Prefer `T[]` because it is more concise.

#### Trailing periods in floating-point constants

Both `1.0` and `1.` create a `Float64` with value `1.0`. Prefer `1.0` over `1.`
because it is more easily distinguished from the integer constant `1`.

Moreover, as recommended by the [Julia style guide](https://docs.julialang.org/en/v1/manual/style-guide/index.html#Avoid-using-floats-for-numeric-literals-in-generic-code-when-possible-1),
never use `1.0` when `1` is okay.

#### Comments

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

#### JuMP macro syntax

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

#### Exported and non-exported names

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

#### Use of underscores within names

The Julia style guide recommends avoiding underscores "when readable", for
example, `haskey`, `isequal`, `remotecall`, and `remotecall_fetch`. This
convention creates the potential for unnecessary bikeshedding and also forces
the user to recall the presence/absence of an underscore, e.g., "was that
argument named `basename` or `base_name`?". For consistency, *always use
underscores* in variable names and function names to separate words.

#### Use of `!`

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
[ordering of function arguments](https://docs.julialang.org/en/v1.0.0/manual/style-guide/#Write-functions-with-argument-ordering-similar-to-Julia's-Base-1).

#### Abbreviations

Abbreviate names to make the code more readable, not to save typing.
Don't arbitrarily delete letters from a word to abbreviate it (e.g., `indx`).
Use abbreviations consistently within a body of code (e.g., do not mix
`con` and `constr`, `idx` and `indx`).

Common abbreviations:

- `num` for `number`
- `con` for `constraint`

TODO: add more

### Miscellaneous

(TODO: Rethink categories.)

#### User-facing `MethodError`

Specifying argument types for methods is mostly optional in Julia, which means
that it's possible to find out that you are working with unexpected types deep in
the call chain. Avoid this situation or handle it with a helpful error message.
*A user should see a `MethodError` only for methods that they called directly.*

Bad:
```julia
internal_function(x::Integer) = x + 1
# The user sees a MethodError for internal_function when calling
# public_function("a string"). This is not very helpful.
public_function(x) = internal_function(x)
```

Good:
```julia
internal_function(x::Integer) = x + 1
# The user sees a MethodError for public_function when calling
# public_function("a string"). This is easy to understand.
public_function(x::Integer) = internal_function(x)
```

If it is hard to provide an error message at the top of the call chain,
then the following pattern is also ok:
```julia
internal_function(x::Integer) = x + 1
function internal_function(x)
    error("Internal error. This probably means that you called " *
          "public_function() with the wrong type.")
end
public_function(x) = internal_function(x)
```

#### `@enum` vs. `Symbol`

The `@enum` macro lets you define types with a finite number of values that
are explicitly enumerated (like `enum` in C/C++). `Symbol`s are lightweight
strings that are used to represent identifiers in Julia (for example, `:x`).

`@enum` provides type safety and can have docstrings attached to explain the
possible values. Use `@enum`s when applicable, e.g., for reporting statuses.
Use strings to provide long-form additional information like error messages.

Use of `Symbol` should typically be reserved for identifiers, e.g., for lookup
in the JuMP model (`model[:my_variable]`).

#### `using` vs. `import`

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

## Documentation

This section describes the writing style that should be used when writing
documentation for JuMP (and supporting packages).

We can recommend the documentation style guides by [Divio](https://www.divio.com/blog/documentation/),
[Google](https://developers.google.com/style/), and [Write the Docs](https://www.writethedocs.org/guide/)
as general reading for those writing documentation. This guide delegates a
thorough handling of the topic to those guides and instead elaborates on the
more points more specific to Julia and documentation that uses
[Documenter](https://github.com/JuliaDocs/Documenter.jl).

 - Be concise.
 - Use lists instead of long sentences.
 - Use numbered lists when describing a sequence, e.g., (1) do X, (2) then Y.
 - Use bullet points when the items are not ordered.
 - Example code should be covered by doctests. (But it's [unclear what to do](https://github.com/JuliaOpt/JuMP.jl/issues/1175)
   if the code depends on a solver.)
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


## Design principles

TODO: How to structure and test large JuMP models, libraries that use JuMP.

For how to write a solver, see MOI.

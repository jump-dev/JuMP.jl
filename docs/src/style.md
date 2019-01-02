Style guide and design principles
=================================

Style guide
-----------

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


#### TODO: Line breaks

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

Both `1.0` and `1.` create a `Float64` with value `1.0`. Prefer `1.0` over
`1.` because it is more easily distinguished from the integer constant `1`.

### Naming

```julia
module SomeModule end
function some_function end
const SOME_CONSTANT = ...
struct SomeStruct end
@enum SomeEnum ENUM_VALUE_A ENUM_VALUE_B
some_local_variable = ...
some_file.jl # Except for ModuleName.jl.
```

#### Use of underscores

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



Design principles
-----------------

TODO: How to structure and test large JuMP models, libraries that use JuMP.

For how to write a solver, see MOI.

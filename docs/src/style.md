Style guide and design principles
=================================

Style guide
-----------

This section describes the coding style rules that apply for JuMP code and that
we recommend for JuMP models and surrounding Julia code. The motivations for
a style guide include:

- conveying best practices for writing readable and maintainable code
- reducing the amount of time spent on
  [bike-shedding](https://en.wikipedia.org/wiki/Law_of_triviality) by
  establishing basic naming and formatting conventions
- lowering the barrier for new contributors by codifying the existing practices
  (e.g., you can be more confident your code will pass review if you follow the style guide)

In some cases, the JuMP style guide diverges from the [Julia style guide](https://docs.julialang.org/en/v0.6.4/manual/style-guide/). All such cases will be explicitly noted and justified.

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

- Use spaces between binary operators
- Use a single space after commas and semicolons
- Do not use extra spaces for unary operators, parentheses, or braces
- Intent within new blocks (except `module`) using 4 spaces

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

#### TODO: Line breaks

### Syntax

Julia sometimes provides equivalent syntax to express the same basic
operation. We discuss these cases below.

#### `for` loops

Julia allows both `for x = 1:N` and `for x in 1:N`. Always prefer to use
`in` over `=`, because `in` generalizes better to other index sets like `for x in eachindex(A)`.

### Empty vectors

For a type `T`, `T[]` and `Vector{T}()` are equivalent ways to create an
empty vector with element type `T`. Prefer `T[]` because it is more concise.

### Trailing periods in floating-point constants

Both `1.0` and `1.` create a `Float64` with value `1.0`. Prefer `1.0` over
`1.` because it is more easily distinguished from the integer constant `1`.



Design principles
-----------------

TODO: How to structure and test large JuMP models, libraries that use JuMP.

For how to write a solver, see MOI.

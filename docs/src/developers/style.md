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
  (for example, you can be more confident your code will pass review if you follow the
  style guide)

In some cases, the JuMP style guide diverges from the
[Julia style guide](https://docs.julialang.org/en/v1.0.0/manual/style-guide/).
All such cases will be explicitly noted and justified.

The JuMP style guide adopts many recommendations from the
[Google style guides](https://github.com/google/styleguide).

!!! info
    The style guide is always a work in progress, and not all JuMP code
    follows the rules. When modifying JuMP, please fix the style violations
    of the surrounding code (that is, leave the code tidier than when you
    started). If large changes are needed, consider separating them into
    another PR.

### JuliaFormatter

JuMP uses [JuliaFormatter.jl](https://github.com/domluna/JuliaFormatter.jl) as
an auto-formatting tool.

We use the options contained in [`.JuliaFormatter.toml`](https://github.com/jump-dev/JuMP.jl/blob/master/.JuliaFormatter.toml).

To format code, `cd` to the JuMP directory, then run:
```julia
] add JuliaFormatter@1
using JuliaFormatter
format("docs")
format("src")
format("test")
```

!!! info
    A continuous integration check verifies that all PRs made to JuMP have
    passed the formatter.

The following sections outline extra style guide points that are not fixed
automatically by JuliaFormatter.

### Abstract types and composition

Specifying types for method arguments is mostly optional in Julia. The benefit
of abstract method arguments is that it enables functions and types from one
package to be used with functions and types from another package via multiple
dispatch.

However, abstractly typed methods have two main drawbacks:

 1. It's possible to find out that you are working with unexpected types deep
    in the call chain, potentially leading to hard-to-diagnose [`MethodError`s](https://docs.julialang.org/en/v1/manual/methods/#Defining-Methods).
 2. Untyped function arguments can lead to correctness problems if the user's
    choice of input type does not satisfy the assumptions made by the author of
    the function.

As a motivating example, consider the following function:
```jldoctest my_sum
julia> function my_sum(x)
           y = 0.0
           for i in 1:length(x)
               y += x[i]
           end
           return y
       end
my_sum (generic function with 1 method)
```
This function contains a number of implicit assumptions about the type of `x`:
 * `x` supports 1-based `getindex` and implements `length`
 * The element type of `x` supports addition with `0.0`, and then with the
   result of `x + 0.0`.

!!! info
    As a motivating example for the second point, [`VariableRef`](@ref) plus
    `Float64` produces an [`AffExpr`](@ref). Do not assume that `+(::A, ::B)`
    produces an instance of the type `A` or `B`.

`my_sum` works as expected if the user passes in `Vector{Float64}`:
```jldoctest my_sum
julia> my_sum([1.0, 2.0, 3.0])
6.0
```
but it doesn't respect input types, for example returning a `Float64` if the
user passes `Vector{Int}`:
```jldoctest my_sum
julia> my_sum([1, 2, 3])
6.0
```
but it throws a `MethodError` if the user passes `String`:
```jldoctest my_sum
julia> my_sum("abc")
ERROR: MethodError: no method matching +(::Float64, ::Char)
[...]
```
This particular `MethodError` is hard to debug, particularly for new users,
because it mentions `+`, `Float64`, and `Char`, none of which were called or
passed by the user.

#### Dealing with `MethodError`s

This section diverges from the [Julia style guide](https://docs.julialang.org/en/v1.6/manual/style-guide/#Avoid-writing-overly-specific-types),
as well as other common guides like [SciML](https://github.com/SciML/SciMLStyle#generic-code-is-preferred-unless-code-is-known-to-be-specific).
The following suggestions are intended to provide a friendlier experience for
novice Julia programmers, at the cost of limiting the power and flexibility of
advanced Julia programmers.

Code should follow the `MethodError` principle:

!!! info "The MethodError principle"
    A user should see a `MethodError` only for methods that they called
    directly.


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
        "`public_function()`s with the wrong type.",
    )
end
public_function(x) = _internal_function(x)
```

#### Dealing with correctness

Dealing with correctness is harder, because Julia has no way of formally
specifying interfaces that abstract types must implement. Instead, here are two
options that you can use when writing and interacting with generic code:

**Option 1: use concrete types and let users extend new methods.**

In this option, explicitly restrict input arguments to concrete types that
are tested and have been validated for correctness. For example:
```jldoctest my_sum
julia> function my_sum_option_1(x::Vector{Float64})
           y = 0.0
           for i in 1:length(x)
               y += x[i]
           end
           return y
       end
my_sum_option_1 (generic function with 1 method)

julia> my_sum_option_1([1.0, 2.0, 3.0])
6.0
```

Using concrete types satisfies the `MethodError` principle:
```jldoctest my_sum
julia> my_sum_option_1("abc")
ERROR: MethodError: no method matching my_sum_option_1(::String)
```
and it allows other types to be supported in future by defining new methods:
```jldoctest my_sum
julia> function my_sum_option_1(x::Array{T,N}) where {T<:Number,N}
           y = zero(T)
           for i in eachindex(x)
               y += x[i]
           end
           return y
       end
my_sum_option_1 (generic function with 2 methods)
```
Importantly, these methods do not have to be defined in the original package.

!!! info
    Some usage of abstract types is okay. For example, in `my_sum_option_1`, we
    allowed the element type, `T`, to be a subtype of `Number`. This is fairly
    safe, but it still has an implicit assumption that `T` supports `zero(T)`
    and `+(::T, ::T)`.

**Option 2: program defensively, and validate all assumptions.**

An alternative is to program defensively, and to rigorously document and
validate all assumptions that the code makes. In particular:

 1. All assumptions on abstract types that aren't guaranteed by the definition
    of the abstract type (for example, optional methods without a fallback)
    should be documented.
 2. If practical, the assumptions should be checked in code, and informative
    error messages should be provided to the user if the assumptions are not
    met. In general, these checks may be expensive, so you should prefer to do
    this once, at the highest level of the call-chain.
 3. Tests should cover for a range of corner cases and argument types.

For example:
```jldoctest my_sum
"""
    test_my_sum_defensive_assumptions(x::AbstractArray{T}) where {T}

Test the assumptions made by `my_sum_defensive`.
"""
function test_my_sum_defensive_assumptions(x::AbstractArray{T}) where {T}
    try
        # Some types may not define zero.
        @assert zero(T) isa T
        # Check iteration supported
        @assert iterate(x) isa Union{Nothing,Tuple{T,Int}}
        # Check that + is defined
        @assert +(zero(T), zero(T)) isa Any
    catch err
        error(
            "Unable to call my_sum_defensive(::$(typeof(x))) because " *
            "it failed an internal assumption",
        )
    end
    return
end

"""
    my_sum_defensive(x::AbstractArray{T}) where {T}

Return the sum of the elements in the abstract array `x`.

## Assumptions

This function makes the following assumptions:

 * That `zero(T)` is defined
 * That `x` supports the iteration interface
 * That  `+(::T, ::T)` is defined
"""
function my_sum_defensive(x::AbstractArray{T}) where {T}
    test_my_sum_defensive_assumptions(x)
    y = zero(T)
    for xi in x
        y += xi
    end
    return y
end

# output

my_sum_defensive
```

This function works on `Vector{Float64}`:
```jldoctest my_sum
julia> my_sum_defensive([1.0, 2.0, 3.0])
6.0
```
as well as `Matrix{Rational{Int}}`:
```jldoctest my_sum
julia> my_sum_defensive([(1//2) + (4//3)im; (6//5) + (7//11)im])
17//10 + 65//33*im
```
and it throws an error when the assumptions aren't met:
```jldoctest my_sum
julia> my_sum_defensive(['a', 'b', 'c'])
ERROR: Unable to call my_sum_defensive(::Vector{Char}) because it failed an internal assumption
[...]
```

As an alternative, you may choose not to call `test_my_sum_defensive_assumptions`
within `my_sum_defensive`, and instead ask users of `my_sum_defensive` to call
it in their tests.

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

The Julia style guide recommends avoiding underscores "when readable," for
example, `haskey`, `isequal`, `remotecall`, and `remotecall_fetch`. This
convention creates the potential for unnecessary bikeshedding and also forces
the user to recall the presence/absence of an underscore, for example, "was that
argument named `basename` or `base_name`?". For consistency, *always use
underscores* in variable names and function names to separate words.

### Use of `!`

Julia has a convention of appending `!` to a function name if the function
modifies its arguments. We recommend to:

- Omit `!` when the name itself makes it clear that modification is taking
  place, for example, `add_constraint` and `set_name`. We depart from the Julia style
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
Don't arbitrarily delete letters from a word to abbreviate it (for example, `indx`).
Use abbreviations consistently within a body of code (for example, do not mix
`con` and `constr`, `idx` and `indx`).

Common abbreviations:

- `num` for `number`
- `con` for `constraint`

### No one-letter variable names

Where possible, avoid one-letter variable names.

Use `model = Model()` instead of `m = Model()`

Exceptions are made for indices in loops.

### `@enum` vs. `Symbol`

The `@enum` macro lets you define types with a finite number of values that
are explicitly enumerated (like `enum` in C/C++). `Symbol`s are lightweight
strings that are used to represent identifiers in Julia (for example, `:x`).

`@enum` provides type safety and can have docstrings attached to explain the
possible values. Use `@enum`s when applicable, for example, for reporting statuses.
Use strings to provide long-form additional information like error messages.

Use of `Symbol` should typically be reserved for identifiers, for example, for lookup
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
 - Use numbered lists when describing a sequence, for example, (1) do X, (2) then Y
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
- Do not terminate lists with punctuation (for example, as in this doc)

Here is an example:
````julia
"""
    signature(args; kwargs...)

Short sentence describing the function.

Optional: add a slightly longer paragraph describing the function.

## Notes

 - List any notes that the user should be aware of

## Example

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

function runtests()
    for name in names(@__MODULE__; all = true)
        if startswith("$(name)", "test_")
            @testset "$(name)" begin
                getfield(@__MODULE__, name)()
            end
        end
    end
    return
end

_helper_function() = 2

function test_addition()
    @test 1 + 1 == _helper_function()
    return
end

end # module TestPkg

TestPkg.runtests()
```

Break the tests into multiple files, with one module per file, so that subsets
of the codebase can be tested by calling `include` with the relevant file.

# Copyright (c) 2019 Arpit Bhatia and contributors                               #src
#                                                                                #src
# Permission is hereby granted, free of charge, to any person obtaining a copy   #src
# of this software and associated documentation files (the "Software"), to deal  #src
# in the Software without restriction, including without limitation the rights   #src
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell      #src
# copies of the Software, and to permit persons to whom the Software is          #src
# furnished to do so, subject to the following conditions:                       #src
#                                                                                #src
# The above copyright notice and this permission notice shall be included in all #src
# copies or substantial portions of the Software.                                #src
#                                                                                #src
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR     #src
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,       #src
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE    #src
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER         #src
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,  #src
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE  #src
# SOFTWARE.                                                                      #src

# # Getting started with Julia

# Because JuMP is embedded in Julia, knowing some basic Julia is important
# before you start learning JuMP.

# !!! tip
#     This tutorial is designed to provide a minimalist crash course in the
#     basics of Julia. You can find resources that provide a more comprehensive
#     introduction to Julia [here](https://julialang.org/learning/).

# ## Installing Julia

# To install Julia, [download the latest stable release](https://julialang.org/downloads/),
# then follow the [platform specific install instructions](https://julialang.org/downloads/platform/).

# !!! tip
#     Unless you know otherwise, you probably want the 64-bit version.

# Next, you need an IDE to develop in. VS Code is a popular choice, so follow
# [these install instructions](https://www.julia-vscode.org/docs/stable/gettingstarted/).

# Julia can also be used with [Jupyter notebooks](https://github.com/JuliaLang/IJulia.jl)
# or the reactive notebooks of [Pluto.jl](https://github.com/fonsp/Pluto.jl).

# ## The Julia REPL

# The main way of interacting with Julia is via its REPL (Read Evaluate Print
# Loop). To access the REPL, start the Julia executable to arrive at the
# `julia>` prompt, and then start coding:

1 + 1

# As your programs become larger, write a script as a text file, and then run
# that file using:
# ```julia
# julia> include("path/to/file.jl")
# ```

# !!! warning
#     Because of Julia's startup latency, running scripts from the command line
#     like the following is slow:
#     ```
#     $ julia path/to/file.jl
#     ```
#     Use the REPL or a notebook instead, and read [The "time-to-first-solve" issue](@ref)
#     for more information.

# ### Code blocks in this documentation

# In this documentation you'll see a mix of code examples with and without the
# `julia>`.

# The Julia prompt is mostly used to demonstrate short code snippets, and the
# output is exactly what you will see if run from the REPL.

# Blocks without the `julia>` can be copy-pasted into the REPL, but they are
# used because they enable richer output like plots or LaTeX to be displayed in
# the online and [PDF](https://jump.dev/JuMP.jl/stable/JuMP.pdf) versions of the
# documentation. If you run them from the REPL you may see different output.

# ## Where to get help

#  * Read the documentation
#    * JuMP [https://jump.dev/JuMP.jl/stable/](https://jump.dev/JuMP.jl/stable/)
#    * Julia [https://docs.julialang.org/en/v1/](https://docs.julialang.org/en/v1/)
#  * Ask (or browse) the Julia community forum: [https://discourse.julialang.org](https://discourse.julialang.org)
#    * If the question is JuMP-related, ask in the [Optimization (Mathematical)](https://discourse.julialang.org/c/domain/opt/13)
#      section, or tag your question with "jump"

# To access the built-in help at the REPL, type `?` to enter help-mode, followed
# by the name of the function to lookup:
# ```julia
# help?> print
# search: print println printstyled sprint isprint prevind parentindices precision escape_string
#
#   print([io::IO], xs...)
#
#   Write to io (or to the default output stream stdout if io is not given) a canonical
#   (un-decorated) text representation. The representation used by print includes minimal formatting
#   and tries to avoid Julia-specific details.
#
#   print falls back to calling show, so most types should just define show. Define print if your
#   type has a separate "plain" representation. For example, show displays strings with quotes, and
#   print displays strings without quotes.
#
#   string returns the output of print as a string.
#
#   Examples
#   ≡≡≡≡≡≡≡≡≡≡
#
#   julia> print("Hello World!")
#   Hello World!
#   julia> io = IOBuffer();
#
#   julia> print(io, "Hello", ' ', :World!)
#
#   julia> String(take!(io))
#   "Hello World!"
# ```

# ## Numbers and arithmetic

# Since we want to solve optimization problems, we're going to be using a lot of
# math. Luckily, Julia is great for math, with all the usual operators:

1 + 1
1 - 2
2 * 2
4 / 5
3^2

# Did you notice how Julia didn't print `.0` after some of the numbers? Julia is
# a dynamic language, which means you never have to explicitly declare the type
# of a variable. However, in the background, Julia is giving each variable a
# type. Check the type of something using the `typeof` function:

typeof(1)
typeof(1.0)

# Here `1` is an `Int64`, which is an integer with 64 bits of precision, and
# `1.0` is a `Float64`, which is a floating point number with 64-bits of
# precision.

# !!! tip
#     If you aren't familiar with floating point numbers, make sure to read
#     the [Floating point numbers](@ref) section.

# We create complex numbers using `im`:

x = 2 + 1im
real(x)
imag(x)
typeof(x)
x * (1 - 2im)

# !!! info
#     The curly brackets surround what we call the _parameters_ of a type. You
#     can read `Complex{Int64}`  as "a complex number, where the real and
#     imaginary parts are represented by `Int64`." If we call
#     `typeof(1.0 + 2.0im)` it will be `Complex{Float64}`, which a complex
#     number with the parts represented by `Float64`.

# There are also some cool things like an irrational representation of π.

π

# !!! tip
#     To make π (and most other Greek letters), type `\pi` and then press
#     `[TAB]`.

typeof(π)

# However, if we do math with irrational numbers, they get converted to
# `Float64`:

typeof(2π / 3)

# ### Floating point numbers

# !!! warning
#     If you aren't familiar with floating point numbers, make sure to read this
#     section carefully.

# A `Float64` is a [floating point](https://en.wikipedia.org/wiki/Floating-point_arithmetic)
# approximation of a real number using 64-bits of information.

# Because it is an approximation, things we know hold true in mathematics don't
# hold true in a computer. For example:

0.1 * 3 == 0.3

# A more complicated example is:

sin(2π / 3) == √3 / 2

# !!! tip
#     Get `√` by typing `\sqrt` then press `[TAB]`.

# Let's see what the differences are:

0.1 * 3 - 0.3
sin(2π / 3) - √3 / 2

# They are small, but not zero.

# One way of explaining this difference is to consider how we would write
# `1 / 3` and `2 / 3` using only four digits after the decimal point. We would
# write `1 / 3` as `0.3333`, and `2 / 3` as `0.6667`. So, despite the fact that
# `2 * (1 / 3) == 2 / 3`, `2 * 0.3333 == 0.6666 != 0.6667`.

# Let's try that again using ≈ (`\approx + [TAB]`) instead of `==`:

0.1 * 3 ≈ 0.3
sin(2π / 3) ≈ √3 / 2

# `≈` is a clever way of calling the `isapprox` function:

isapprox(sin(2π / 3), √3 / 2; atol = 1e-8)

# !!! warning
#     Floating point is the reason solvers use tolerances when they solve
#     optimization models. A common mistake you're likely to make is checking
#     whether a binary variable is 0 using `value(z) == 0`. Always remember to
#     use something like `isapprox` when comparing floating point numbers.

# Note that `isapprox` will always return `false` if one of the number being
# compared is `0` and `atol` is zero (its default value).

1e-300 ≈ 0.0

# so always set a nonzero value of `atol` if one of the arguments can be zero.

isapprox(1e-9, 0.0; atol = 1e-8)

# !!! tip
#     Gurobi has a [good series of articles](https://www.gurobi.com/documentation/9.0/refman/num_grb_guidelines_for_num.html)
#     on the implications of floating point in optimization if you want to read
#     more.

# If you aren't careful, floating point arithmetic can throw up all manner of
# issues. For example:

1 + 1e-16 == 1

# It even turns out that floating point numbers aren't associative:

(1 + 1e-16) - 1e-16 == 1 + (1e-16 - 1e-16)

# It's important to note that this issue isn't Julia-specific. It happens in
# every programming language (try it out in Python).

# ## Vectors, matrices and arrays

# Similar to MATLAB, Julia has native support for vectors, matrices and tensors;
# all of which are represented by arrays of different dimensions. Vectors are
# constructed by comma-separated elements surrounded by square brackets:

b = [5, 6]

# Matrices can be constructed with spaces separating the columns, and semicolons
# separating the rows:

A = [1.0 2.0; 3.0 4.0]

# We can do linear algebra:

x = A \ b

# !!! info
#     Here is floating point at work again; `x` is approximately `[-4, 4.5]`.

A * x
A * x ≈ b

# Note that when multiplying vectors and matrices, dimensions matter. For
# example, you can't multiply a vector by a vector:

try                         #hide
    b * b
catch err                   #hide
    showerror(stderr, err)  #hide
end                         #hide

# But multiplying transposes works:

b' * b
b * b'

# ## Other common types

# ### Comments

# Although not technically a type, code comments begin with the `#` character:

1 + 1  # This is a comment

# Multiline comments begin with `#=` and end with `=#`:
# ```julia
# #=
# Here is a
# multiline comment
# =#
# ```

# Comments can even be nested inside expressions. This is sometimes helpful when
# documenting inputs to functions:

isapprox(
    sin(π),
    0.0;
    #= We need an explicit atol here because we are comparing with 0 =#
    atol = 0.001,
)

# ### Strings

# Double quotes are used for strings:

typeof("This is Julia")

# Unicode is fine in strings:

typeof("π is about 3.1415")

# Use `println` to print a string:

println("Hello, World!")

# Use `$()` to interpolate values into a string:

x = 123
println("The value of x is: $(x)")

# Use triple-quotes for multiline strings:

s = """
Here is
a
multiline string
"""

println(s)

# ### Symbols

# Julia `Symbol`s are a data structure from the compiler that represent Julia
# identifiers (that is, variable names).

println("The value of x is: $(eval(:x))")

# !!! warning
#     We used `eval` here to demonstrate how Julia links `Symbol`s to variables.
#     However, avoid calling `eval` in your code. It is usually a sign that your
#     code is doing something that could be more easily achieved a different
#     way. The [Community Forum](https://jump.dev/forum) is a good place to ask
#     for advice on alternative approaches.

typeof(:x)

# You can think of a `Symbol` as a `String` that takes up less memory, and that
# can't be modified.

# Convert between `String` and `Symbol` using their constructors:

String(:abc)
Symbol("abc")

# !!! tip
#     `Symbol`s are often (ab)used to stand in for a `String` or an `Enum`, when
#     one of the latter is likely a better choice. The JuMP [Style guide](@ref)
#     recommends reserving `Symbol`s for identifiers. See [@enum vs. Symbol](@ref)
#     for more.

# ### Tuples

# Julia makes extensive use of a simple data structure called Tuples. Tuples are
# immutable collections of values. For example:

t = ("hello", 1.2, :foo)
typeof(t)

# Tuples can be accessed by index, similar to arrays:

t[2]

# And they can be "unpacked" like so:

a, b, c = t
b

# The values can also be given names, which is a convenient way of making
# light-weight data structures.

t = (word = "hello", num = 1.2, sym = :foo)

# Values can be accessed using dot syntax:

t.word

# ## Dictionaries

# Similar to Python, Julia has native support for dictionaries. Dictionaries
# provide a very generic way of mapping keys to values.  For example, a map of
# integers to strings:

d1 = Dict(1 => "A", 2 => "B", 4 => "D")

# !!! info
#     Type-stuff again: `Dict{Int64,String}` is a dictionary with `Int64` keys
#     and `String` values.

# Looking up a value uses the bracket syntax:

d1[2]

# Dictionaries support non-integer keys and can mix data types:

Dict("A" => 1, "B" => 2.5, "D" => 2 - 3im)

# !!! info
#     Julia types form a hierarchy. Here the value type of the dictionary is
#     `Number`, which is a generalization of `Int64`, `Float64`, and
#     `Complex{Int}`. Leaf nodes in this hierarchy are called "concrete" types,
#     and all others are called "Abstract." In general, having variables with
#     abstract types like `Number` can lead to slower code, so you should try to
#     make sure every element in a dictionary or vector is the same type. For
#     example, in this case we could represent every element as a
#     `Complex{Float64}`:

Dict("A" => 1.0 + 0.0im, "B" => 2.5 + 0.0im, "D" => 2.0 - 3.0im)

# Dictionaries can be nested:

d2 = Dict("A" => 1, "B" => 2, "D" => Dict(:foo => 3, :bar => 4))
d2["B"]
d2["D"][:foo]

# ## Structs

# You can define custom datastructures with `struct`:

struct MyStruct
    x::Int
    y::String
    z::Dict{Int,Int}
end

a = MyStruct(1, "a", Dict(2 => 3))
a.x

# By default, these are not mutable

try                         #hide
    a.x = 2
catch err                   #hide
    showerror(stderr, err)  #hide
end                         #hide

# However, you can declare a `mutable struct` which is mutable:

mutable struct MyStructMutable
    x::Int
    y::String
    z::Dict{Int,Int}
end

a = MyStructMutable(1, "a", Dict(2 => 3))
a.x
a.x = 2
a

# ## Loops

# Julia has native support for for-each style loops with the syntax
# `for <value> in <collection> end`:

for i in 1:5
    println(i)
end

# !!! info
#     Ranges are constructed as `start:stop`, or `start:step:stop`.

for i in 1.2:1.1:5.6
    println(i)
end

# This for-each loop also works with dictionaries:

for (key, value) in Dict("A" => 1, "B" => 2.5, "D" => 2 - 3im)
    println("$(key): $(value)")
end

# Note that in contrast to vector languages like MATLAB and R, loops do not
# result in a significant performance degradation in Julia.

# ## Control flow

# Julia control flow is similar to MATLAB, using the keywords
# `if-elseif-else-end`, and the logical operators `||` and `&&` for **or** and
# **and** respectively:

for i in 0:5:15
    if i < 5
        println("$(i) is less than 5")
    elseif i < 10
        println("$(i) is less than 10")
    else
        if i == 10
            println("the value is 10")
        else
            println("$(i) is bigger than 10")
        end
    end
end

# ## Comprehensions

# Similar to languages like Haskell and Python, Julia supports the use of simple
# loops in the construction of arrays and dictionaries, called comprehensions.
#
# A list of increasing integers:

[i for i in 1:5]

# Matrices can be built by including multiple indices:

[i * j for i in 1:5, j in 5:10]

# Conditional statements can be used to filter out some values:

[i for i in 1:10 if i % 2 == 1]

# A similar syntax can be used for building dictionaries:

Dict("$(i)" => i for i in 1:10 if i % 2 == 1)

# ## Functions

# A simple function is defined as follows:

function print_hello()
    return println("hello")
end
print_hello()

# Arguments can be added to a function:

function print_it(x)
    return println(x)
end
print_it("hello")
print_it(1.234)
print_it(:my_id)

# Optional keyword arguments are also possible:

function print_it(x; prefix = "value:")
    return println("$(prefix) $(x)")
end
print_it(1.234)
print_it(1.234; prefix = "val:")

# The keyword `return` is used to specify the return values of a function:

function mult(x; y = 2.0)
    return x * y
end

mult(4.0)
mult(4.0; y = 5.0)

# ### Anonymous functions

# The syntax `input -> output` creates an anonymous function. These are most
# useful when passed to other functions. For example:

f = x -> x^2
f(2)
map(x -> x^2, 1:4)

# ### Type parameters

# We can constrain the inputs to a function using type parameters, which are
# `::` followed by the type of the input we want. For example:

function foo(x::Int)
    return x^2
end

function foo(x::Float64)
    return exp(x)
end

function foo(x::Number)
    return x + 1
end

foo(2)
foo(2.0)
foo(1 + 1im)

# But what happens if we call `foo` with something we haven't defined it for?

try                         #hide
    foo([1, 2, 3])
catch err                   #hide
    showerror(stderr, err)  #hide
end                         #hide

# A `MethodError` means that you passed a
# function something that didn't match the type that it was expecting. In this
# case, the error message says that it doesn't know how to handle an
# `Vector{Int64}`, but it does know how to handle `Float64`, `Int64`, and
# `Number`.
#
# !!! tip
#     Read the "Closest candidates" part of the error message carefully to get a
#     hint as to what was expected.

# ### Broadcasting

# In the example above, we didn't define what to do if `f` was passed a
# `Vector`. Luckily, Julia provides a convenient syntax for mapping `f`
# element-wise over arrays. Just add a `.` between the name of the function and
# the opening `(`. This works for _any_ function, including functions with
# multiple arguments. For example:

f.([1, 2, 3])

# !!! tip
#     Get a `MethodError` when calling a function that takes a `Vector`,
#     `Matrix`, or `Array`? Try broadcasting.

# ## Mutable vs immutable objects

# Some types in Julia are *mutable*, which means you can change the values
# inside them. A good example is an array. You can modify the contents of an
# array without having to make a new array.

# In contrast, types like `Float64` are *immutable*. You cannot modify the
# contents of a `Float64`.

# This is something to be aware of when passing types into functions. For
# example:

function mutability_example(mutable_type::Vector{Int}, immutable_type::Int)
    mutable_type[1] += 1
    immutable_type += 1
    return
end

mutable_type = [1, 2, 3]
immutable_type = 1

mutability_example(mutable_type, immutable_type)

println("mutable_type: $(mutable_type)")
println("immutable_type: $(immutable_type)")

# Because `Vector{Int}` is a mutable type, modifying the variable inside the
# function changed the value outside of the function. In contrast, the change
# to `immutable_type` didn't modify the value outside the function.

# You can check mutability with the `isimmutable` function:

isimmutable([1, 2, 3])
isimmutable(1)

# ## The package manager

# ### Installing packages

# No matter how wonderful Julia's base language is, at some point you will want
# to use an extension package.  Some of these are built-in, for example random
# number generation is available in the `Random` package in the standard
# library. These packages are loaded with the commands `using` and `import`.

using Random  # The equivalent of Python's `from Random import *`
import Random  # The equivalent of Python's `import Random`

Random.seed!(33)

[rand() for i in 1:10]

# The Package Manager is used to install packages that are not part of Julia's
# standard library.

# For example the following can be used to install JuMP,
# ```julia
# using Pkg
# Pkg.add("JuMP")
# ```

# For a complete list of registered Julia packages see the package listing at
# [JuliaHub](https://juliahub.com).

# From time to you may wish to use a Julia package that is not registered.  In
# this case a git repository URL can be used to install the package.
# ```julia
# using Pkg
# Pkg.add("https://github.com/user-name/MyPackage.jl.git")
# ```

# ### Package environments

# By default, `Pkg.add` will add packages to Julia's global environment.
# However, Julia also has built-in support for virtual environments.

# Activate a virtual environment with:
# ```julia
# import Pkg; Pkg.activate("/path/to/environment")
# ```

# You can see what packages are installed in the current environment with
# `Pkg.status()`.

# !!! tip
#     We _strongly_ recommend you create a Pkg environment for each project
#     that you create in Julia, and add only the packages that you need, instead
#     of adding lots of packages to the global environment. The [Pkg manager documentation](https://julialang.github.io/Pkg.jl/v1/environments/)
#     has more information on this topic.

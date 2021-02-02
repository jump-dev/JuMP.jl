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

# **Originally Contributed by**: Juan Pablo Vielma

# Since JuMP is embedded in Julia, knowing some basic Julia is important
# for learning JuMP. This tutorial is designed to provide a minimalist
# crash course in the basics of Julia. You can find resources that provide
# a more comprehensive introduction to Julia [here](https://julialang.org/learning/).

# ## How to print

# In Julia, we usually use println() to print

println("Hello, World!")

# ## Basic data types

# Integers

typeof(1 + -2)

# Floating point numbers

typeof(1.2 - 2.3)

# There are also some cool things like an irrational representation of π. To
# make π (and most other greek letters), type `\pi` and then press `[TAB]`.

π

#-

typeof(π)

# Julia has native support for complex numbers

typeof(2 + 3im)

# Double quotes are used for strings

typeof("This is Julia")

# Unicode is fine in strings

typeof("π is about 3.1415")

# Julia symbols provide a way to make human readable unique identifiers.

:my_id
typeof(:my_id)

# ## Arithmetic and Equality Testing

# Julia is great for math

1 + 1

# Even math involving complex numbers

(2 + 1im) * (1 - 2im)

# We can also write things like the following using √ (\sqrt)

sin(2π/3) == √3/2

# Wait. What???

sin(2π/3) - √3/2

# Let's try again using ≈ (\approx).

sin(2π/3) ≈ √3/2

# Note that this time we used ≈ instead of ==. That is because computers don't
# use real numbers. They use a discrete representation called floating point. If
# you aren't careful, this can throw up all manner of issues. For example:

1 + 1e-16 == 1

# It even turns out that floating point numbers aren't associative!

(1 + 1e-16) - 1e-16 == 1 + (1e-16 - 1e-16)

# ## Vectors, matrices and arrays

# Similar to Matlab, Julia has native support for vectors, matrices and tensors;
# all of which are represented by arrays of different dimensions. Vectors are
# constructed by comma-separated elements surrounded by square brackets:

b = [5, 6]

# Matrices can by constructed with spaces separating the columns, and semicolons
# separating the rows:

A = [1 2; 3 4]

# We can do linear algebra:

x = A \ b

#-

A * x

#-

A * x == b

# Note that when multiplying vectors and matrices, dimensions matter. For
# example, you can't multiply a vector by a vector:

try  #hide
b * b
catch err; showerror(stderr, err); end  #hide

# But multiplying transposes works:

b' * b

#-

b * b'

# ## Tuples

# Julia makes extensive use of a simple data structure called Tuples.  Tuples
# are immutable collections of values. For example,

t = ("hello", 1.2, :foo)

#-

typeof(t)

# Tuples can be accessed by index, similar to arrays,

t[2]

# And can be "unpacked" like so,

a, b, c = t
b

# The values can also be given names, which is a convenient way of making
# light-weight data structures.

t = (word = "hello", num = 1.2, sym = :foo)

# Then values can be accessed using a dot syntax,

t.word

# ## Dictionaries

# Similar to Python, Julia has native support for dictionaries.  Dictionaries
# provide a very generic way of mapping keys to values.  For example, a map of
# integers to strings,

d1 = Dict(1 => "A", 2 => "B", 4 => "D")

# Looking up a values uses the bracket syntax,

d1[2]

# Dictionaries support non-integer keys and can mix data types,

Dict("A" => 1, "B" => 2.5, "D" => 2 - 3im)

# Dictionaries can be nested

d2 = Dict("A" => 1, "B" => 2, "D" => Dict(:foo => 3, :bar => 4))

#-

d2["B"]

#-

d2["D"][:foo]

# ## For-Each Loops

# Julia has native support for for-each style loops with the syntax
# `for <value> in <collection> end`.

for i in 1:5
    println(i)
end

#-

for i in [1.2, 2.3, 3.4, 4.5, 5.6]
    println(i)
end

# This for-each loop also works with dictionaries.

for (key, value) in Dict("A" => 1, "B" => 2.5, "D" => 2 - 3im)
    println("$key: $value")
end

# Note that in contrast to vector languages like Matlab and R, loops do not
# result in a significant performance degradation in Julia.

# ## Control Flow

# Julia control flow is similar to Matlab, using the keywords
# `if-elseif-else-end`, and the logical operators `||` and `&&` for *or* and
# *and* respectively.

i = 10
for i in 0:3:15
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
# loops in the construction of arrays and dictionaries, called comprehenions.
#
# A list of increasing integers,

[i for i in 1:5]

# Matrices can be built by including multiple indices,

[i*j for i in 1:5, j in 5:10]

# Conditional statements can be used to filter out some values,

[i for i in 1:10 if i%2 == 1]

# A similar syntax can be used for building dictionaries

Dict("$i" => i for i in 1:10 if i%2 == 1)

# ## Functions

# A simple function is defined as follows:

function print_hello()
    println("hello")
end
print_hello()

# Arguments can be added to a function,

function print_it(x)
    println(x)
end
print_it("hello")
print_it(1.234)
print_it(:my_id)

# Optional keyword arguments are also possible

function print_it(x; prefix="value:")
    println("$(prefix) $x")
end
print_it(1.234)
print_it(1.234, prefix="val:")

# The keyword `return` is used to specify the return values of a function.

function mult(x; y=2.0)
    return x * y
end
mult(4.0)

#-

mult(4.0, y=5.0)

# ## Other notes on types

# Usually, specifing types is not required to use Julia.  However, it can be
# helpful to understand the basics of Julia types for debugging.
# For example this list has a type of `Array{Int64,1}` indicating that it is a
# one dimensional array of integer values.

[1, 5, -2, 7]

# In this example, the decimal values lead to a one dimensional array of
# floating point values, i.e. `Array{Float64,1}`.  Notice that the integer `7`
# is promoted to a `Float64`, because all elements in the array need share a
# common type.

[1.0, 5.2, -2.1, 7]

# ## Mutable vs immutable objects

# Some types in Julia are *mutable*, which means you can change the values
# inside them. A good example is an array. You can modify the contents of an
# array without having to make a new array.

# In contrast, types like `Float64` are *immutable*. You can't modify the
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
# function changed the value outside of the function. In constrast, the change
# to `immutable_type` didn't modify the value outside the function.

# You can check mutability with the `isimmutable` function.

isimmutable([1, 2, 3])

#-

isimmutable(1)

# ## Using Packages and the Package Manager

# No matter how wonderful Julia's base language is, at some point you will want
# to use an extension package.  Some of these are built-in, for example random
# number generation is available in the `Random` package in the standard
# library. These packages are loaded with the commands `using` and `import`.

using Random

Random.seed!(33)

[rand() for i in 1:10]

# The Package Manager is used to install packages that are not part of Julia's
# standard library.

# For example the following can be used to install JuMP,
# ```julia
# using Pkg
# Pkg.add("JuMP")
# ```

# For a complete list of registed Julia packages see the package listing at
# [JuliaHub](https://juliahub.com).

# From time to you may wish to use a Julia package that is not registered.  In
# this case a git repository URL can be used to install the package.
# ```julia
# using Pkg
# Pkg.add("https://github.com/user-name/MyPackage.jl.git")
# ```

# Note that for clarity this example uses the package manager `Pkg`.  Julia
# includes an interactive package manager that can be accessed using `]`.
# [This video](https://youtu.be/76KL8aSz0Sg) gives an overview of using the
# interactive package manager environment.

# The state of installed packages can also be saved in two files: `Project.toml`
# and `Manifest.toml`. If these files are stored in the directory `/tmp/jump`,
# the state of the packages can be recovered by running

# ```julia
# import Pkg
# Pkg.activate("/tmp/jump")
# Pkg.instantiate()
# ```

# ## HELP!

# Julia includes a help mode that can be accessed using `?`.  Entering any
# object (e.g. function, type, struct, ...) into the help mode will show its
# documentation, if any is available.

# ## Some Common Gotchas

# ### MethodError

# A common error in Julia is `MethodError`, which indicates that the function is
# not defined for the given value.  For example, by default the `ceil` function
# is not defined for complex numbers.  The "closest candidates" list suggest
# some Julia types that the function is defined for.

try  #hide
ceil(1.2 + 2.3im)
catch err; showerror(stderr, err); end  #hide

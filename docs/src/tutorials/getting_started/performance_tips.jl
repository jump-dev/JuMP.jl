# Copyright (c) 2021 Oscar Dowson and contributors                               #src
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

# # Performance tips

# By now you should have read the other "getting started" tutorials. You're
# almost ready to write your own models, but before you do so there are some
# important things to be aware of.

using JuMP
import HiGHS

# ## Read the Julia performance tips

# The first thing to do is read the [Performance tips](https://docs.julialang.org/en/v1/manual/performance-tips/index.html)
# section of the Julia manual. The most important rule is to avoid global
# variables. This is particularly important if you're learning JuMP after using
# a language like MATLAB.

# ## Use macros to build expressions

# Use JuMP's macros (or [`add_to_expression!`](@ref)) to build expressions.
# Avoid constructing expressions outside the macros.

# Constructing an expression outside the macro results in intermediate copies of
# the expression. For example,
# ```julia
# x[1] + x[2] + x[3]
# ```
# is equivalent to
# ```julia
# a = x[1]
# b = a + x[2]
# c = b + x[3]
# ```
# Since we only care about `c`, the `a` and `b` expressions are not needed and
# constructing them slows the program down.

# JuMP's macros rewrite the expressions to operate in-place and avoid these
# extra copies. Because they allocate less memory, they are faster, particularly
# for large expressions.

# Here's an example.

model = Model()
@variable(model, x[1:3])

# Here's what happens if we construct the expression outside the macro:

@allocated x[1] + x[2] + x[3]

# !!! info
#     The `@allocated` measures how many bytes were allocated during the
#     evaluation of an expression. Fewer is better.

# If we use the [`@expression`](@ref) macro, we get many fewer allocations:

@allocated @expression(model, x[1] + x[2] + x[3])

# ## Disable string names

# By default, JuMP creates `String` names for variables and constraints and
# passes these to the solver. The benefit of passing names is that it improves
# the readability of log messages from the solver (for example, "variable x has
# invalid bounds" instead of "variable v1203 has invalid bounds"), but for
# larger models the overhead of passing names can be non-trivial.

# Disable the creation of `String` names by setting `set_string_name = false` in
# the [`@variable`](@ref) and [`@constraint`](@ref) macros, or by calling
# [`set_string_names_on_creation`](@ref) to disable all names for a particular
# model:

model = Model();
set_string_names_on_creation(model, false)
@variable(model, x)
@constraint(model, c, 2x <= 1)

# Note that this doesn't change how symbolic names and bindings are stored:

x
model[:x]
x === model[:x]

# But you can no longer look up the variable by the string name:

variable_by_name(model, "x") === nothing

# !!! info
#     For more information on the difference between string names, symbolic
#     names, and bindings, see
#     [String names, symbolic names, and bindings](@ref variable_names_and_bindings).

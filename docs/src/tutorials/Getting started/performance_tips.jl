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

# This tutorial explains how to write performant JuMP code.

# !!! tip
#     Read the [Performance tips](https://docs.julialang.org/en/v1/manual/performance-tips/index.html)
#     section of the Julia manual. The most important rule is to avoid global
#     variables!

using JuMP  # hide
using GLPK  # hide

# ## The "time-to-first-solve" issue

# Similar to the infamous [time-to-first-plot](https://discourse.julialang.org/t/roadmap-for-a-faster-time-to-first-plot/22956)
# plotting problem, JuMP suffers from time-to-first-solve latency. This latency
# occurs because the first time you call JuMP code in each session, Julia needs
# to compile a lot of code specific to your problem. This issue is actively being
# worked on, but there are a few things you can do to improve things.

# ### Don't call JuMP from the command line

# In other languages, you might be used to a workflow like:
# ```
# $ julia my_script.jl
# ```
# This doesn't work for JuMP, because we have to pay the compilation latency
# every time you run the script. Instead, use one of the [suggested workflows](https://docs.julialang.org/en/v1/manual/workflow-tips/)
# from the Julia documentation.

# ### Use PackageCompiler

# As a final option, consider using [PackageCompiler.jl](https://julialang.github.io/PackageCompiler.jl/dev/)
# to [create a custom sysimage](https://julialang.github.io/PackageCompiler.jl/dev/examples/plots/).

# This is a good option if you have finished prototyping a model, and you now
# want to call it frequently from the command line without paying the
# compilation price.

# ## Use macros to build expressions

# ### What

# Use JuMP's macros (or [`add_to_expression!`](@ref)) to build expressions.
# Avoid constructing expressions outside the macros.

# ### Why

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
# constructing them slows the program down!

# JuMP's macros rewrite the expressions to operate in-place and avoid these
# extra copies. Because they allocate less memory, they are faster, particularly
# for large expressions.

# ### Example

model = Model()
@variable(model, x[1:3])

# Here's what happens if we construct the expression outside the macro:

@allocated x[1] + x[2] + x[3]

# !!! info
#     The `@allocated` measures how many bytes were allocated during the
#     evaluation of an expression. Fewer is better.

# If we use the [`@expression`](@ref) macro, we get many fewer allocations:

@allocated @expression(model, x[1] + x[2] + x[3])

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

# The Julia manual has an excellent section on [Performance tips](https://docs.julialang.org/en/v1/manual/performance-tips/index.html).
# The purpose of this tutorial is to highlight a number of performance issues
# that are specific to JuMP.

# ## Required packages

# This tutorial uses the following packages:

using JuMP
import HiGHS

# ## Use macros to build expressions

# Use JuMP's macros to build expressions.

# Constructing an expression outside the macro results in intermediate copies of
# the expression. For example,
# ```julia
# result = x[1] + x[2] + x[3]
# ```
# is equivalent to
# ```julia
# tmp = x[1] + x[2]
# result = tmp + x[3]
# ```
# Since we only care about `result`, the `tmp` expression is not needed, and
# constructing it slows the program down.

# JuMP's macros rewrite the expressions to operate in-place and avoid temporary
# expressions. Because they allocate less memory, they are faster, particularly
# for large expressions.

# Here's an example.

model = Model()
@variable(model, x[1:3])
x[1] + x[2] + x[3]                                                 #hide
@assert @allocated(x[1] + x[2] + x[3]) > 1000                      #hide
@expression(model, x[1] + x[2] + x[3])                             #hide
@assert @allocated(@expression(model, x[1] + x[2] + x[3])) < 1000  #hide
@allocated x[1] + x[2] + x[3]
@allocated @expression(model, x[1] + x[2] + x[3])

# ## Use `add_to_expression!` to build summations

# If you don't want to use the expression macros, use
# [`add_to_expression!`](@ref) to build summations. For example, instead of:

expr = zero(AffExpr)
for i in 1:3
    global expr  #hide
    expr += x[i]
end
expr

# do

expr = zero(AffExpr)
for i in 1:3
    add_to_expression!(expr, x[i])
end
expr

# The former is equivalent to:

expr0 = zero(AffExpr)
expr1 = expr0 + x[1]
expr2 = expr1 + x[2]
expr = expr2 + x[3]

# which allocates four unique [`AffExpr`](@ref) objects. The latter efficiently
# updates `expr` in-place so that only one [`AffExpr`](@ref) object is
# allocated.

# The function [`add_to_expression!`](@ref) also supports terms like `y += a * x`
# where `a` is a constant. For example, instead of:

expr = zero(AffExpr)
for i in 1:3
    global expr  #hide
    expr += i * x[i]
end
expr

# do

expr = zero(AffExpr)
for i in 1:3
    add_to_expression!(expr, i, x[i])
end
expr

# Don't do this, because `i * x[i]` will allocate a new [`AffExpr`](@ref) in
# each iteration:

expr = zero(AffExpr)
for i in 1:3
    add_to_expression!(expr, i * x[i])
end
expr

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

# ## Sparsity

# JuMP executes the code that you write without trying to be clever about
# multiplication by `0.0`. For example, consider this model:

d = 100_000
a = zeros(d);
a[3] = 1.0
model = Model();
@variable(model, x[1:d]);
complicated_expression(x, i) = x[i]^2
@expression(model, sum(a[i] * complicated_expression(x, i) for i in 1:d))

# Although the final expression consists of  a single element, the sum is over
# 100,000 elements, all but one of which are then multiplied by `0.0`. The
# `@expression` line is equivalent to:

expr = zero(QuadExpr)
for i in 1:d
    tmp = complicated_expression(x, i)
    global expr = add_to_expression!(expr, a[i], tmp)
end

# Notice how we compute `complicated_expression` in every iteration, even though
# most results will be discarded. You can improve the performance of model
# construction by pre-computing the set of non-zero indices:

indices = [i for i in 1:d if !iszero(a[i])]
@expression(model, sum(a[i] * complicated_expression(x, i) for i in indices))

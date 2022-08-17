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

# # Debugging

# Dealing with bugs is an unavoidable part of coding optimization models in
# JuMP. This includes bugs related to general Julia code such as syntax errors,
# method errors, typos, and off-by-one indexing errors, but it also includes
# optimization-specific bugs related to the formulation and solution of your
# model.

# This tutorial explains some common sources of bugs and modeling issues that
# you might encounter when writing models in JuMP, and it suggests a variety of
# strategies to deal with them.

# !!! tip
#     This tutorial is more advanced than the other "Getting started" tutorials.
#     It's in the "Getting started" section to give you an early preview of how
#     to test and debug JuMP models. However, if you are new to JuMP, you may
#     want to briefly skim the tutorial, and come back to it once you have
#     written a few JuMP models.

using JuMP
import HiGHS

# ## General rules for debugging

#

# Before all else, read the [Debugging chapter](https://benlauwens.github.io/ThinkJulia.jl/latest/book.html#chap21)
# in the book [ThinkJulia.jl](https://benlauwens.github.io/ThinkJulia.jl/latest/book.html).

#

# * Simplify the problem

# ## Debugging an infeasible model

# ## Debugging an unbounded model

# A model is unbounded if there is no limit on how good the objective value can
# get. In general, an unbounded model means that you have an error in your
# modeling, because all physical systems have limits. (You cannot make an
# infinite amount of profit.)

# A simple example of an unbounded model is:

model = Model(HiGHS.Optimizer)
set_silent(model)
@variable(model, x >= 0)
@objective(model, Max, 2x + 1)

# because we can increase `x` without limit, and the objective value `2x + 1`
# gets better as `x` increases.

# JuMP doesn't have an `UNBOUNDED` termination status. Instead, unbounded models
# will return `DUAL_INFEASIBLE`:

optimize!(model)
termination_status(model)

# Common sources of unboundedness are:
#
#  * Using `Max` instead of `Min`
#  * Omitting variable bounds, such as `0 <= x <= 1`
#  * Using `+` instead of `-` in a term of the objective function.

# Strategies to debug sources of unboundedness are:

#  * Double check whether you intended `Min` or `Max` in the [`@objective`](@ref)
#    line.
#  * Print the objective function with `print(objective_function(model))` and
#    verify that the value and sign of each coefficient is as you expect.
#  * Add large bounds to all variables that are free or have one-sided bounds,
#    then re-solve the problem. Because all variables are now bounded, the
#    problem will have a finite optimal solution. Look at the value of each
#    variable in the optimal solution to see if it is at one of the new bounds.
#    If it is, you either need to specify a better bound for that variable, or
#    there might be a mistake in the objective function associated with that
#    variable (for example, a `+` instead of a `-`).

# If there are too many variables to add bounds to, or there are too many terms
# to examine by hand, another strategy is to create a new variable with a large
# upper bound (if maximizing, lower bound if minimizing) and a constraint that
# the variable must be less-than or equal to the expression of the objective
# function. For example:

model = Model(HiGHS.Optimizer)
set_silent(model)
@variable(model, x >= 0)
## @objective(model, Max, 2x + 1)
@variable(model, objective <= 10_000)
@constraint(model, objective <= 2x + 1)
@objective(model, Max, objective)

# This new model has a finite optimal solution, so we can solve it and then look
# for variables with large positive or negative values in the optimal solution.

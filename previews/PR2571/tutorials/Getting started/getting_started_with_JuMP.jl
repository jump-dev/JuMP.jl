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

# # Getting started with JuMP

# This tutorial is aimed at providing a quick introduction to writing JuMP code.
# If you're new to Julia, you should start with [Getting started with Julia](@ref).

# ## What is JuMP?

# JuMP ("Julia for Mathematical Programming") is an open-source modeling
# language that is embedded in Julia. It allows users to users formulate various
# classes of optimization problems (linear, mixed-integer, quadratic, conic
# quadratic, semidefinite, and nonlinear) with easy-to-read code. These problems
# can then be solved using state-of-the-art open-source and commercial solvers.

# JuMP also makes advanced optimization techniques easily accessible from a
# high-level language.

# ## Installation

# JuMP is a package for Julia. From Julia, JuMP is installed by using the
# built-in package manager.

# ```julia
# import Pkg
# Pkg.add("JuMP")
# ```

# You also need to include a Julia package which provides an appropriate solver.
# One such solver is `GLPK.Optimizer`, which is provided by the
# [GLPK.jl package](https://github.com/JuliaOpt/GLPK.jl).
# ```julia
# import Pkg
# Pkg.add("GLPK")
# ```
# See [Installation Guide](@ref) for a list of other solvers you can use.

# ## An example

# Let's try to solve the following linear programming problem by using JuMP and
# GLPK. We will first look at the complete code to solve the problem and then go
# through it step by step.

# ```math
# \begin{aligned}
# & \min & 12x + 20y \\
# & \;\;\text{s.t.} & 6x + 8y \geq 100 \\
# & & 7x + 12y \geq 120 \\
# & & x \geq 0 \\
# & & y \in [0, 3] \\
# \end{aligned}
# ```

using JuMP
using GLPK
model = Model(GLPK.Optimizer)
@variable(model, x >= 0)
@variable(model, 0 <= y <= 3)
@objective(model, Min, 12x + 20y)
@constraint(model, c1, 6x + 8y >= 100)
@constraint(model, c2, 7x + 12y >= 120)
print(model)
optimize!(model)
@show termination_status(model)
@show primal_status(model)
@show dual_status(model)
@show objective_value(model)
@show value(x)
@show value(y)
@show shadow_price(c1)
@show shadow_price(c2)
nothing #hide

# ## Step-by-step

# Once JuMP is installed, to use JuMP in your programs, we just need to write:

using JuMP

# We also need to include a Julia package which provides an appropriate solver.
# We want to use `GLPK.Optimizer` here which is provided by the `GLPK.jl`
# package.

using GLPK

# A model object is a container for variables, constraints, solver options, etc.
# Models are created with the [`Model`](@ref) function. The model can be created
# with an optimizer attached with default arguments by calling the constructor
# with the optimizer type, as follows:

model = Model(GLPK.Optimizer)

# Variables are modeled using [`@variable`](@ref):

@variable(model, x >= 0)

# They can have lower and upper bounds.

@variable(model, 0 <= y <= 30)

# The objective is set using [`@objective`](@ref):

@objective(model, Min, 12x + 20y)

# Constraints are modeled using [`@constraint`](@ref). Here `c1` and `c2` are
# the names of our constraint.

@constraint(model, c1, 6x + 8y >= 100)

#-

@constraint(model, c2, 7x + 12y >= 120)

#- Call `print` to display the model:

print(model)

# To solve the optimization problem, call the [`optimize!`] function.

optimize!(model)

# !!! info
#     The `!` after optimize is just part of the name. It's nothing special.
#     Julia has a convention that functions which mutate their arguments should
#     end in `!`. A common example is `push!`.

# Now let's see what information we can query about the solution.

# [`termination_status`](@ref) tells us why the solver stopped:

termination_status(model)

# In this case, the solver found an optimal solution. We should also check
# [`primal_status`](@ref) to see if the solver found a primal feasible point:

primal_status(model)

# and [`dual_status`](@ref) to see if the solver found a dual feasible point:

dual_status(model)

# Now we know that our solver found an optimal solution, and has a primal and a
# dual solution to query.

# Query the objective value using [`objective_value`](@ref):

objective_value(model)

# The primal solution using [`value`](@ref):

value(x)

#-

value(y)

# and the dual solution using [`shadow_price`](@ref):

shadow_price(c1)

#-

shadow_price(c2)

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

# This tutorial is aimed at providing a quick introduction to writing and
# solving optimization models with JuMP.

# If you're new to Julia, start by reading [Getting started with Julia](@ref).

# ## What is JuMP?

# JuMP ("Julia for Mathematical Programming") is an open-source modeling
# language that is embedded in Julia. It allows users to formulate various
# classes of optimization problems (linear, mixed-integer, quadratic, conic
# quadratic, semidefinite, and nonlinear) with easy-to-read code. These problems
# can then be solved using state-of-the-art open-source and commercial solvers.

# JuMP also makes advanced optimization techniques easily accessible from a
# high-level language.

# ## What is a solver?

# A solver is a software package that incorporates algorithms for finding
# solutions to one or more classes of problem.

# For example, HiGHS is a solver for linear programming (LP) and mixed integer
# programming (MIP) problems. It incorporates algorithms such as the simplex
# method and the interior-point method.

# The [Supported-solvers](@ref) table lists the open-source and commercial
# solvers that JuMP currently supports.

# ## What is MathOptInterface?

# Each solver has its own concepts and data structures for representing
# optimization models and obtaining results.

# [MathOptInterface](https://github.com/jump-dev/MathOptInterface.jl) (MOI) is
# an abstraction layer that JuMP uses to convert from the problem written in
# JuMP to the solver-specific data structures for each solver.

# MOI can be used directly, or through a higher-level modeling interface like
# JuMP.

# Because JuMP is built on top of MOI, you'll often see the `MathOptInterface.`
# prefix displayed when JuMP types are printed. However, you'll only need to
# understand and interact with MOI to accomplish advanced tasks such as creating
# [solver-independent callbacks](@ref callbacks_manual).

# ## Installation

# JuMP is a package for Julia. From Julia, JuMP is installed by using the
# built-in package manager.

# ```julia
# import Pkg
# Pkg.add("JuMP")
# ```

# You also need to include a Julia package which provides an appropriate solver.
# One such solver is `HiGHS.Optimizer`, which is provided by the
# [HiGHS.jl package](https://github.com/jump-dev/HiGHS.jl).
# ```julia
# import Pkg
# Pkg.add("HiGHS")
# ```
# See [Installation Guide](@ref) for a list of other solvers you can use.

# ## An example

# Let's solve the following linear programming problem using JuMP and HiGHS.
# We will first look at the complete code to solve the problem and then go
# through it step by step.

# Here's the problem:
# ```math
# \begin{aligned}
# & \min & 12x + 20y \\
# & \;\;\text{s.t.} & 6x + 8y \geq 100 \\
# & & 7x + 12y \geq 120 \\
# & & x \geq 0 \\
# & & y \in [0, 3] \\
# \end{aligned}
# ```

# And here's the code to solve this problem:

using JuMP
using HiGHS
model = Model(HiGHS.Optimizer)
@variable(model, x >= 0)
@variable(model, 0 <= y <= 3)
@objective(model, Min, 12x + 20y)
@constraint(model, c1, 6x + 8y >= 100)
@constraint(model, c2, 7x + 12y >= 120)
print(model)
optimize!(model)
termination_status(model)
primal_status(model)
dual_status(model)
objective_value(model)
value(x)
value(y)
shadow_price(c1)
shadow_price(c2)

# ## Step-by-step

# Once JuMP is installed, to use JuMP in your programs write:

using JuMP

# We also need to include a Julia package which provides an appropriate solver.
# We want to use `HiGHS.Optimizer` here which is provided by the `HiGHS.jl`
# package:

using HiGHS

# JuMP builds problems incrementally in a `Model` object. Create a model by
# passing an optimizer to the [`Model`](@ref) function:

model = Model(HiGHS.Optimizer)

# Variables are modeled using [`@variable`](@ref):

@variable(model, x >= 0)

# !!! info
#     The macro creates a new Julia object, `x`, in the current scope. We could
#     have made this more explicit by writing:
#     ```julia
#     x = @variable(model, x >= 0)
#     ```
#     but the macro does this automatically for us to save writing `x` twice.

# Variables can have lower and upper bounds:

@variable(model, 0 <= y <= 30)

# The objective is set using [`@objective`](@ref):

@objective(model, Min, 12x + 20y)

# Constraints are modeled using [`@constraint`](@ref). Here, `c1` and `c2` are
# the names of our constraint:

@constraint(model, c1, 6x + 8y >= 100)
@constraint(model, c2, 7x + 12y >= 120)

# Call `print` to display the model:

print(model)

# To solve the optimization problem, call the [`optimize!`](@ref) function:

optimize!(model)

# !!! info
#     The `!` after optimize is part of the name. It's nothing special.
#     Julia has a convention that functions which mutate their arguments should
#     end in `!`. A common example is `push!`.

# Now let's see what information we can query about the solution.

# [`termination_status`](@ref) tells us why the solver stopped:

termination_status(model)

# In this case, the solver found an optimal solution.

# Check [`primal_status`](@ref) to see if the solver found a primal feasible
# point:

primal_status(model)

# and [`dual_status`](@ref) to see if the solver found a dual feasible point:

dual_status(model)

# Now we know that our solver found an optimal solution, and that it has a
# primal and a dual solution to query.

# Query the objective value using [`objective_value`](@ref):

objective_value(model)

# the primal solution using [`value`](@ref):

value(x)
value(y)

# and the dual solution using [`shadow_price`](@ref):

shadow_price(c1)
shadow_price(c2)

# That's it for our simple model. In the rest of this tutorial, we expand on
# some of the basic JuMP operations.

# ## Model basics

# Create a model by passing an optimizer:

model = Model(HiGHS.Optimizer)

# Alternatively, call [`set_optimizer`](@ref) at any point before calling
# [`optimize!`](@ref):

model = Model()
set_optimizer(model, HiGHS.Optimizer)

# For some solvers, you can also use [`direct_model`](@ref), which offers a more
# efficient connection to the underlying solver:

model = direct_model(HiGHS.Optimizer())

# !!! warning
#     Some solvers do not support [`direct_model`](@ref)!

# ### Solver Options

# Pass options to solvers with [`optimizer_with_attributes`](@ref):

model =
    Model(optimizer_with_attributes(HiGHS.Optimizer, "output_flag" => false))

# !!! note
#     These options are solver-specific. To find out the various options
#     available, see the GitHub README of the individual solver packages. The
#     link to each solver's GitHub page is in the [Supported solvers](@ref)
#     table.

# You can also pass options with [`set_attribute`](@ref):

model = Model(HiGHS.Optimizer)
set_attribute(model, "output_flag", false)

# ## Solution basics

# We saw above how to use [`termination_status`](@ref) and
# [`primal_status`](@ref) to understand the solution returned by the solver.

# However, only query solution attributes like [`value`](@ref) and
# [`objective_value`](@ref) if there is an available solution. Here's a
# recommended way to check:

function solve_infeasible()
    model = Model(HiGHS.Optimizer)
    @variable(model, 0 <= x <= 1)
    @variable(model, 0 <= y <= 1)
    @constraint(model, x + y >= 3)
    @objective(model, Max, x + 2y)
    optimize!(model)
    if termination_status(model) != OPTIMAL
        @warn("The model was not solved correctly.")
        return
    end
    return value(x), value(y)
end

solve_infeasible()

# ## Variable basics

# Let's create a new empty model to explain some of the variable syntax:

model = Model()

# ### Variable bounds

# All of the variables we have created till now have had a bound. We can also
# create a free variable.

@variable(model, free_x)

# While creating a variable, instead of using the <= and >= syntax, we can also
# use the `lower_bound` and `upper_bound` keyword arguments.

@variable(model, keyword_x, lower_bound = 1, upper_bound = 2)

# We can query whether a variable has a bound using the `has_lower_bound` and
# `has_upper_bound` functions. The values of the bound can be obtained using the
# `lower_bound` and `upper_bound` functions.

has_upper_bound(keyword_x)
upper_bound(keyword_x)

# Note querying the value of a bound that does not exist will result in an error.

try                         #hide
    lower_bound(free_x)
catch err                   #hide
    showerror(stderr, err)  #hide
end                         #hide

# ### [Containers](@id tutorial_variable_container)

# We have already seen how to add a single variable to a model using the
# [`@variable`](@ref) macro. Now let's look at ways to add multiple variables to
# a model.

# JuMP provides data structures for adding collections of variables to a model.
# These data structures are referred to as _containers_ and are of three types:
# `Arrays`, `DenseAxisArrays`, and `SparseAxisArrays`.

# #### Arrays

# JuMP arrays are created when you have integer indices that start at `1`:

@variable(model, a[1:2, 1:2])

# Index elements in `a` as follows:

a[1, 1]

#-

a[2, :]

# Create an n-dimensional variable $x \in {R}^n$ with bounds $l \le x \le u$
# ($l, u \in {R}^n$) as follows:

n = 10
l = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10];
u = [10, 11, 12, 13, 14, 15, 16, 17, 18, 19];

@variable(model, l[i] <= x[i = 1:n] <= u[i])

# We can also create variable bounds that depend upon the indices:

@variable(model, y[i = 1:2, j = 1:2] >= 2i + j)

# #### DenseAxisArrays

# `DenseAxisArrays` are used when the indices are not one-based integer ranges.
# The syntax is similar except with an arbitrary vector as an index as opposed
# to a one-based range:

@variable(model, z[i = 2:3, j = 1:2:3] >= 0)

# Indices do not have to be integers. They can be any Julia type:

@variable(model, w[1:5, ["red", "blue"]] <= 1)

# Index elements in a `DenseAxisArray` as follows:

z[2, 1]

#-

w[2:3, ["red", "blue"]]

# See [Forcing the container type](@ref variable_forcing) for more details.

# #### SparseAxisArrays

# `SparseAxisArrays` are created when the indices do not form a Cartesian product.
# For example, this applies when indices have a dependence upon previous indices
# (called triangular indexing):

@variable(model, u[i = 1:2, j = i:3])

# We can also conditionally create variables by adding a comparison check that
# depends upon the named indices and is separated from the indices by a
# semi-colon `;`:

@variable(model, v[i = 1:9; mod(i, 3) == 0])

# Index elements in a `DenseAxisArray` as follows:

u[1, 2]

#-

v[[3, 6]]

# ### Integrality

# JuMP can create binary and integer variables. Binary variables are constrained
# to the set  ``\{0, 1\}``, and integer variables are constrained to the set
# ``\mathbb{Z}``.

# #### Integer variables

# Create an integer variable by passing  `Int`:

@variable(model, integer_x, Int)

# or setting the `integer` keyword to `true`:

@variable(model, integer_z, integer = true)

# #### Binary variables

# Create a binary variable by passing `Bin`:

@variable(model, binary_x, Bin)

# or setting the `binary` keyword to `true`:

@variable(model, binary_z, binary = true)

# ## Constraint basics

# We'll need a new model to explain some of the constraint  basics:

model = Model()
@variable(model, x)
@variable(model, y)
@variable(model, z[1:10]);

# ### [Containers](@id tutorial_constraint_container)

# Just as we had containers for variables, JuMP also provides `Arrays`,
# `DenseAxisArrays`, and `SparseAxisArrays` for storing collections of
# constraints. Examples for each container type are given below.

# #### Arrays

# Create an `Array` of constraints:

@constraint(model, [i = 1:3], i * x <= i + 1)

# #### DenseAxisArrays

# Create an `DenseAxisArray` of constraints:

@constraint(model, [i = 1:2, j = 2:3], i * x <= j + 1)

# #### SparseAxisArrays

# Create an `SparseAxisArray` of constraints:

@constraint(model, [i = 1:2, j = 1:2; i != j], i * x <= j + 1)

# ### Constraints in a loop

# We can add constraints using regular Julia loops:

for i in 1:3
    @constraint(model, 6x + 4y >= 5i)
end

# or use for each loops inside the `@constraint` macro:

@constraint(model, [i in 1:3], 6x + 4y >= 5i)

# We can also create constraints such as $\sum _{i = 1}^{10} z_i \leq 1$:

@constraint(model, sum(z[i] for i in 1:10) <= 1)

# ## Objective functions

# Set an objective function with [`@objective`](@ref):

model = Model(HiGHS.Optimizer)
@variable(model, x >= 0)
@variable(model, y >= 0)
@objective(model, Min, 2x + y)

# Create a maximization objective using `Max`:

@objective(model, Max, 2x + y)

# !!! tip
#     Calling [`@objective`](@ref) multiple times will over-write the previous
#     objective. This can be useful when you want to solve the same problem with
#     different objectives.

# ## Vectorized syntax

# We can also add constraints and an objective to JuMP using vectorized linear
# algebra. We'll illustrate this by solving an LP in standard form that is,

# ```math
# \begin{aligned}
# & \min & c^T x \\
# & \;\;\text{s.t.} & A x = b \\
# & & x \ge 0
# \end{aligned}
# ```

vector_model = Model(HiGHS.Optimizer)
A = [1 1 9 5; 3 5 0 8; 2 0 6 13]
b = [7, 3, 5]
c = [1, 3, 5, 2]
@variable(vector_model, x[1:4] >= 0)
@constraint(vector_model, A * x .== b)
@objective(vector_model, Min, c' * x)
optimize!(vector_model)
objective_value(vector_model)

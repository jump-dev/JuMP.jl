# Copyright (c) 2022 Oscar Dowson and contributors                               #src
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

# # Nested optimization problems

# In this tutorial we explain how to write a *nested* optimization problem,
# where an *outer* problem uses the results from the optimization of an *inner*
# subproblem.
#
# To model the problem, we define a user-defined function to handle the decomposition
# of the inner problem inside the outer one. Finally, we show how to improve
# the performance by using a cache that avoids resolving the inner problem.
#
# For a simpler example of writing a user-defined function,
# see the [User-defined Hessians](@ref) tutorial.

# This tutorial uses the following packages:

using JuMP
import Ipopt

# In the rest of this tutorial, our goal is to solve the bilevel
# optimization problem:

# ```math
# \begin{array}{r l}
# \min\limits_{x,z} & x_1^2 + x_2^2 + z \\
# s.t.            & \begin{array}{r l}
#                       z = \max\limits_{y} & x_1^2 y_1 + x_2^2 y_2  - x_1 y_1^4 - 2 x_2 y_2^4 \\
#                       s.t.                  & (y_1 - 10)^2 + (y_2 - 10)^2 \le 25
#                   \end{array} \\
#                 & x \ge 0.
# \end{array}
# ```

# This bilevel optimization problem is composed of two nested optimization
# problems. An _upper_ level, involving variables ``x``, and a _lower_ level,
# involving variables ``y``. From the perspective of the lower-level problem,
# the values of ``x`` are fixed parameters, and so the model optimizes ``y``
# given those fixed parameters. Simultaneously, the upper-level problem
# optimizes ``x`` and ``z`` given the response of ``y``.

# ## Decomposition

# There are a few ways to solve this problem, but we are going to use a
# nonlinear decomposition method. The first step is to write a function to
# compute the lower-level problem:

# ```math
# \begin{array}{r l}
#   V(x_1, x_2) = \max\limits_{y} & x_1^2 y_1 + x_2^2 y_2  - x_1 y_1^4 - 2 x_2 y_2^4 \\
#                            s.t. & (y_1 - 10)^2 + (y_2 - 10)^2 \le 25
# \end{array}
# ```

function solve_lower_level(x...)
    model = Model(Ipopt.Optimizer)
    set_silent(model)
    @variable(model, y[1:2])
    @NLobjective(
        model,
        Max,
        x[1]^2 * y[1] + x[2]^2 * y[2] - x[1] * y[1]^4 - 2 * x[2] * y[2]^4,
    )
    @constraint(model, (y[1] - 10)^2 + (y[2] - 10)^2 <= 25)
    optimize!(model)
    @assert termination_status(model) == LOCALLY_SOLVED
    return objective_value(model), value.(y)
end

# The next function takes a value of ``x`` and returns the optimal lower-level
# objective-value and the optimal response ``y``. The reason why we need both
# the objective and the optimal ``y`` will be made clear shortly, but for now
# let us define:

function V(x...)
    f, _ = solve_lower_level(x...)
    return f
end

# Then, we can substitute ``V`` into our full problem to create:

# ```math
# \begin{array}{r l}
# \min\limits_{x} & x_1^2 + x_2^2 + V(x_1, x_2) \\
# s.t.            & x \ge 0.
# \end{array}
# ```

# This looks like a nonlinear optimization problem with a user-defined function
# ``V``! However, because ``V`` solves an optimization problem internally, we
# can't use automatic differentiation to compute the first and second
# derivatives. Instead, we can use JuMP's ability to pass callback functions
# for the gradient and Hessian instead.

# First up, we need to define the gradient of ``V`` with respect to ``x``. In
# general, this may be difficult to compute, but because ``x`` appears only in
# the objective, we can just differentiate the objective function with respect
# to ``x``, giving:

function ∇V(g::AbstractVector, x...)
    _, y = solve_lower_level(x...)
    g[1] = 2 * x[1] * y[1] - y[1]^4
    g[2] = 2 * x[2] * y[2] - 2 * y[2]^4
    return
end

# Second, we need to define the Hessian of ``V`` with respect to ``x``. This is
# a symmetric matrix, but in our example only the diagonal elements are
# non-zero:

function ∇²V(H::AbstractMatrix, x...)
    _, y = solve_lower_level(x...)
    H[1, 1] = 2 * y[1]
    H[2, 2] = 2 * y[2]
    return
end

# !!! info
#     Providing an explicit Hessian function is optional
#     if first derivatives are already available.

# We now have enough to define our bilevel optimization problem:

model = Model(Ipopt.Optimizer)
@variable(model, x[1:2] >= 0)
register(model, :V, 2, V, ∇V, ∇²V)
@NLobjective(model, Min, x[1]^2 + x[2]^2 + V(x[1], x[2]))
optimize!(model)
solution_summary(model)

# The optimal objective value is:

objective_value(model)

# and the optimal upper-level decision variables ``x`` are:

value.(x)

# To compute the optimal lower-level decision variables, we need to call
# `solve_lower_level` with the optimal upper-level decision variables:

_, y = solve_lower_level(value.(x)...)
y

# ## Improving performance

# Our solution approach works, but it has a performance problem: every time
# we need to compute the value, gradient, or Hessian of ``V``, we have to
# re-solve the lower-level optimization problem. This is wasteful, because we
# will often call the gradient and Hessian at the same point, and so solving the
# problem twice with the same input repeats work unnecessarily.

# We can work around this by using a cache:

mutable struct Cache
    x::Any
    f::Float64
    y::Vector{Float64}
end

# with a function to update the cache if needed:

function _update_if_needed(cache::Cache, x...)
    if cache.x !== x
        cache.f, cache.y = solve_lower_level(x...)
        cache.x = x
    end
    return
end

# Then, we define cached versions of out three functions which call
# `_updated_if_needed` and return values from the cache.

function cached_f(cache::Cache, x...)
    _update_if_needed(cache, x...)
    return cache.f
end

function cached_∇f(cache::Cache, g::AbstractVector, x...)
    _update_if_needed(cache, x...)
    g[1] = 2 * x[1] * cache.y[1] - cache.y[1]^4
    g[2] = 2 * x[2] * cache.y[2] - 2 * cache.y[2]^4
    return
end

function cached_∇²f(cache::Cache, H::AbstractMatrix, x...)
    _update_if_needed(cache, x...)
    H[1, 1] = 2 * cache.y[1]
    H[2, 2] = 2 * cache.y[2]
    return
end

# Now we're ready to setup and solve the upper level optimization problem:

model = Model(Ipopt.Optimizer)
@variable(model, x[1:2] >= 0)
cache = Cache(Float64[], NaN, Float64[])
register(
    model,
    :V,
    2,
    (x...) -> cached_f(cache, x...),
    (g, x...) -> cached_∇f(cache, g, x...),
    (H, x...) -> cached_∇²f(cache, H, x...),
)
@NLobjective(model, Min, x[1]^2 + x[2]^2 + V(x[1], x[2]))
optimize!(model)
solution_summary(model)

# an we can check we get the same objective value:

objective_value(model)

# and upper-level decision variable ``x``:

value.(x)

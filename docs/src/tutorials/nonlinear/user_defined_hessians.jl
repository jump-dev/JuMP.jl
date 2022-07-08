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

# # User-defined hessians

# In this tutorial, we explain how to write a user-defined function with an
# explicit hessian.

# This tutorial uses the following packages:

using JuMP
import Ipopt

# Our goal for this tutorial is to solve the bilevel optimization problem:

# ```math
# \begin{array}{r l}
# \min\limits_{x} & x_1^2 + x_2^2 + z \\
# s.t.            & \begin{array}{r l}
#                       z \ge \max\limits_{y} & x_1^2 y_1 + x_2^2 y_2  - x_1 y_1^4 - 2 x_2 y_2^4 \\
#                       s.t.                  & (y_1 - 10)^2 + (y_2 - 10)^2 \le 25
#                   \end{array} \\
#                 & x \ge 0.
# \end{array}
# ```

# This bilevel optimization problem is composed of two nested optimization
# problems. An _upper_ level, involving variables ``x``, and a _lower_ level,
# involving variables ``y``. From the perspective of the lower-level, the
# values of ``x`` are fixed parameters, and so the model optimizes ``y`` given
# those fixed parameters. Simultaneously, the upper level is optimizing ``x``
# given the response of ``yy``.

# There are a few ways to solve this problem, but we are going to use a
# nonlinear decomposition method. The first step is to write a function to
# compute:

# ```math
# \begin{array}{r l}
#   V(x_1, x_z) = \max\limits_{y} & x_1^2 y_1 + x_2^2 y_2  - x_1 y_1^4 - 2 x_2 y_2^4 \\
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

# This function takes a guess of ``x`` and returns the optimal lower-level
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
# for the gradient and hessian instead.

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

# Second, we need to define the hessian of ``V`` with respect to ``x``. This is
# a symmetric matrix, but in our example only the diagonal elements are
# non-zero:

function ∇²V(H::AbstractMatrix, x...)
    _, y = solve_lower_level(x...)
    H[1, 1] = 2 * y[1]
    H[2, 2] = 2 * y[2]
    return
end

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

# This solution approach worked, but it has a performance problem: every time
# we needed to compute the value, gradient, or hessian of ``V``, we had to
# re-solve the lower-level optimization problem! This is wasteful, because we
# will often call the gradient and hessian at the same point, and so solving the
# problem twice with the same input repeats work unnecessarily.

# We can work around this by using memoization:

function memoized_solve_lower_level()
    last_x, f, y = nothing, NaN, [NaN, NaN]
    function _update_if_needed(x...)
        if last_x != x
            f, y = solve_lower_level(x...)
            last_x = x
        end
        return
    end
    function memoized_f(x...)
        _update_if_needed(x...)
        return f
    end
    function memoized_∇f(g::AbstractVector, x...)
        _update_if_needed(x...)
        g[1] = 2 * x[1] * y[1] - y[1]^4
        g[2] = 2 * x[2] * y[2] - 2 * y[2]^4
        return
    end
    function memoized_∇²f(H::AbstractMatrix, x...)
        _update_if_needed(x...)
        H[1, 1] = 2 * y[1]
        H[2, 2] = 2 * y[2]
        return
    end
    return memoized_f, memoized_∇f, memoized_∇²f
end

f, ∇f, ∇²f = memoized_solve_lower_level()

# The function above is a little confusing, but it returns three new functions
# `f`, `∇f`, and `∇²f`, each of which call `_update_if_needed(x...)`. This
# function only updates the cached values of `f` and `y` if the input `x` is
# different to what is last saw.

model = Model(Ipopt.Optimizer)
@variable(model, x[1:2] >= 0)
register(model, :V, 2, f, ∇f, ∇²f)
@NLobjective(model, Min, x[1]^2 + x[2]^2 + V(x[1], x[2]))
optimize!(model)
solution_summary(model)

# an we can check we get the same objective value:

objective_value(model)

# and upper-level decision variable ``x``:

value.(x)

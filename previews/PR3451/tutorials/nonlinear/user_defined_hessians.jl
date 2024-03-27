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

# # User-defined Hessians

# In this tutorial, we explain how to write a user-defined function (see
# [User-defined functions](@ref)) with a Hessian matrix explicitly provided by
# the user.
#
# For a more advanced example, see [Nested optimization problems](@ref).

# This tutorial uses the following packages:

using JuMP
import Ipopt

# ## Rosenbrock example

# As a simple example, we consider the Rosenbrock function:

rosenbrock(x...) = (1 - x[1])^2 + 100 * (x[2] - x[1]^2)^2

# which has the gradient vector:

function ∇rosenbrock(g::AbstractVector, x...)
    g[1] = 400 * x[1]^3 - 400 * x[1] * x[2] + 2 * x[1] - 2
    g[2] = 200 * (x[2] - x[1]^2)
    return
end

# and the Hessian matrix:

function ∇²rosenbrock(H::AbstractMatrix, x...)
    H[1, 1] = 1200 * x[1]^2 - 400 * x[2] + 2
    ## H[1, 2] = -400 * x[1] <-- not needed because Hessian is symmetric
    H[2, 1] = -400 * x[1]
    H[2, 2] = 200.0
    return
end

# You may assume the Hessian matrix `H` is initialized with zeros, and
# because it is symmetric you need only to fill in the non-zero of the
# lower-triangular terms.

# The matrix type passed in as `H` depends on the automatic differentiation
# system, so make sure the first argument to the Hessian function supports an
# `AbstractMatrix` (it may be something other than `Matrix{Float64}`). However,
# you may assume only that `H` supports `size(H)` and `setindex!`.

# Now that we have the function, its gradient, and its Hessian, we can construct
# a JuMP model, register the function, and use it in a macro:

model = Model(Ipopt.Optimizer)
@variable(model, x[1:2])
@register(model, f_rosenbrock, 2, rosenbrock, ∇rosenbrock, ∇²rosenbrock)
@objective(model, Min, f_rosenbrock(x[1], x[2]))
optimize!(model)
solution_summary(model; verbose = true)

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

# # Tips and Tricks

# **Originally Contributed by**: Arpit Bhatia

# This tutorial is aimed at providing a simplistic introduction to conic
# programming using JuMP.

# It uses the following packages:

using JuMP
import SCS
import LinearAlgebra

# !!! info
#     This tutorial uses sets from [MathOptInterface](@ref moi_documentation).
#     By default, JuMP exports the `MOI` symbol as an alias for the
#     MathOptInterface.jl package. We recommend making this more explicit in
#     your code by adding the following lines:
#     ```julia
#     import MathOptInterface
#     const MOI = MathOptInterface
#     ```

import Random      # hide
Random.seed!(1234) # hide
nothing            # hide

# !!! tip
#     A good resource for learning more about functions which can be modeled
#     using cones is the [MOSEK Modeling Cookbook](https://docs.mosek.com/modeling-cookbook/index.html).

# ## What is a cone?

# A subset $C$ of a vector space $V$ is a cone if $\forall x \in C$ and positive
# scalars $\lambda > 0$, the product $\lambda x \in C$.

# A cone $C$ is a convex cone if $\lambda x + (1 - \lambda) y \in C$, for any
# $\lambda \in [0, 1]$, and any $x, y \in C$.

# ## What is a conic program?

# Conic programming problems are convex optimization problems in which a convex
# function is minimized over the intersection of an affine subspace and a convex
# cone. An example of a conic-form minimization problems, in the primal form is:

# ```math
# \begin{aligned}
# & \min_{x \in \mathbb{R}^n} & a_0^T x + b_0 \\
# & \;\;\text{s.t.} & A_i x + b_i & \in \mathcal{C}_i & i = 1 \ldots m
# \end{aligned}
# ```

# The corresponding dual problem is:

# ```math
# \begin{aligned}
# & \max_{y_1, \ldots, y_m} & -\sum_{i=1}^m b_i^T y_i + b_0 \\
# & \;\;\text{s.t.} & a_0 - \sum_{i=1}^m A_i^T y_i & = 0 \\
# & & y_i & \in \mathcal{C}_i^* & i = 1 \ldots m
# \end{aligned}
# ```
#
# where each $\mathcal{C}_i$ is a closed convex cone and $\mathcal{C}_i^*$ is
# its dual cone.

# ## Second-Order Cone

# The Second-Order Cone (or Lorentz Cone) of dimension $n$ is of the form:
# ```math
# Q^n = \{ (t, x) \in \mathbb{R}^n : t \ge ||x||_2 \}
# ```

# ### Example

# Minimize the L2 norm of a vector $x$.

model = Model()
@variable(model, x[1:3])
@variable(model, norm_x)
@constraint(model, [norm_x; x] in SecondOrderCone())
@objective(model, Min, norm_x)

# ## Rotated Second-Order Cone

# A Second-Order Cone rotated by $\pi/4$ in the $(x_1,x_2)$ plane is called a
# Rotated Second-Order Cone. It is of the form:
# ```math
# Q_r^n = \{ (t,u,x) \in \mathbb{R}^n : 2tu \ge ||x||_2^2, t,u \ge 0 \}
# ```

# ### Example

# Given a set of predictors $x$, and observations $y$, find the parameter
# $\theta$ that minimizes the sum of squares loss between $y_i$ and
# $\theta x_i$.

x = [1.0, 2.0, 3.0, 4.0]
y = [0.45, 1.04, 1.51, 1.97]
model = Model()
@variable(model, θ)
@variable(model, loss)
@constraint(model, [loss; 0.5; θ .* x .- y] in RotatedSecondOrderCone())
@objective(model, Min, loss)

# ## Exponential Cone

# An Exponential Cone is a set of the form:

# ```math
# K_{exp} = \{ (x,y,z) \in \mathbb{R}^3 : y \exp (x/y) \le z, y > 0 \}
# ```

model = Model()
@variable(model, x[1:3] >= 0)
@constraint(model, x in MOI.ExponentialCone())
@objective(model, Min, x[3])

# ### Example: Entropy Maximization

# The entropy maximization problem consists of maximizing the entropy function,
# $H(x) = -x\log{x}$ subject to linear inequality constraints.

# ```math
# \begin{aligned}
# & \max & - \sum_{i=1}^n x_i \log x_i \\
# & \;\;\text{s.t.} & \mathbf{1}' x = 1 \\
# & & Ax \leq b
# \end{aligned}
# ```

# We can model this problem using an exponential cone by using the following
# transformation:

# ```math
# t\leq -x\log{x} \iff t\leq x\log(1/x)  \iff (t, x, 1) \in K_{exp}
# ```

# Thus, our problem becomes,

# ```math
# \begin{aligned}
# & \max & 1^Tt \\
# & \;\;\text{s.t.} & Ax \leq b \\
# & & 1^T x = 1 \\
# & & (t_i, x_i, 1) \in K_{exp} && \forall i = 1 \ldots n \\
# \end{aligned}
# ```

n = 15
m = 10
A = randn(m, n)
b = rand(m, 1)

model = Model(SCS.Optimizer)
set_silent(model)
@variable(model, t[1:n])
@variable(model, x[1:n])
@objective(model, Max, sum(t))
@constraint(model, sum(x) == 1)
@constraint(model, A * x .<= b)
@constraint(model, con[i = 1:n], [t[i], x[i], 1] in MOI.ExponentialCone())
optimize!(model)

#-

objective_value(model)

# ### Positive Semidefinite Cone

# The set of positive semidefinite matrices (PSD) of dimension $n$ form a cone
# in $\mathbb{R}^n$. We write this set mathematically as:

# ```math
# \mathcal{S}_{+}^n = \{ X \in \mathcal{S}^n \mid z^T X z \geq 0, \: \forall z\in \mathbb{R}^n \}.
# ```

# A PSD cone is represented in JuMP using the MOI sets
# `PositiveSemidefiniteConeTriangle` (for upper triangle of a PSD matrix) and
# `PositiveSemidefiniteConeSquare` (for a complete PSD matrix). However, it is
# preferable to use the `PSDCone` shortcut as illustrated below.

# #### Example: largest eigenvalue of a symmetric matrix

# Suppose $A$ has eigenvalues $\lambda_{1} \geq \lambda_{2} \ldots \geq \lambda_{n}$.
# Then the matrix $t I-A$ has eigenvalues $t-\lambda_{1}, t-\lambda_{2}, \ldots, t-\lambda_{n}$.
# Note that $t I-A$ is PSD exactly when all these eigenvalues are non-negative,
# and this happens for values $t \geq \lambda_{1}$. Thus, we can model the
# problem of finding the largest eigenvalue of a symmetric matrix as:

# ```math
# \begin{aligned}
# \lambda_{1} = \min t \\
# \text { s.t. } t I-A \succeq 0
# \end{aligned}
# ```

A = [3 2 4; 2 0 2; 4 2 3]
I = Matrix{Float64}(LinearAlgebra.I, 3, 3)
model = Model(SCS.Optimizer)
set_silent(model)
@variable(model, t)
@objective(model, Min, t)
@constraint(model, t .* I - A in PSDCone())

optimize!(model)

#-

objective_value(model)

# ## Other Cones and Functions

# For other cones supported by JuMP, check out the
# [MathOptInterface Manual](https://jump.dev/MathOptInterface.jl/stable).

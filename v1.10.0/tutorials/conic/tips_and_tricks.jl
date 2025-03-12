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

# # [Tips and Tricks](@id conic_tips_and_tricks)

# **This tutorial was originally contributed by Arpit Bhatia.**

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
#     import MathOptInterface as MOI
#     ```

import Random      # hide
Random.seed!(1234) # hide
nothing            # hide

# !!! tip
#     A good resource for learning more about functions which can be modeled
#     using cones is the [MOSEK Modeling Cookbook](https://docs.mosek.com/modeling-cookbook/index.html).

# ## Background theory

# A subset $C$ of a vector space $V$ is a cone if $\forall x \in C$ and positive
# scalars $\lambda > 0$, the product $\lambda x \in C$.

# A cone $C$ is a convex cone if $\lambda x + (1 - \lambda) y \in C$, for any
# $\lambda \in [0, 1]$, and any $x, y \in C$.

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

# The [`SecondOrderCone`](@ref) (or Lorentz Cone) of dimension $n$ is a cone of
# the form:
# ```math
# K_{soc} = \{ (t, x) \in \mathbb{R}^n : t \ge ||x||_2 \}
# ```

# It is most commonly used to represent the L2-norm of the vector $x$:

model = Model(SCS.Optimizer)
set_silent(model)
@variable(model, x[1:3])
@variable(model, t)
@constraint(model, sum(x) == 1)
@constraint(model, [t; x] in SecondOrderCone())
@objective(model, Min, t)
optimize!(model)
value(t), value.(x)

# ## Rotated Second-Order Cone

# A Second-Order Cone rotated by $\pi/4$ in the $(x_1,x_2)$ plane is called a
# [`RotatedSecondOrderCone`](@ref). It is a cone of the form:
# ```math
# K_{rsoc} = \{ (t,u,x) \in \mathbb{R}^n : 2tu \ge ||x||_2^2, t,u \ge 0 \}
# ```

# When `u = 0.5`, it represents the sum of squares of a vector $x$:

data = [1.0, 2.0, 3.0, 4.0]
target = [0.45, 1.04, 1.51, 1.97]
model = Model(SCS.Optimizer)
set_silent(model)
@variable(model, θ)
@variable(model, t)
@variable(model, u == 0.5)
@expression(model, residuals, θ * data .- target)
@constraint(model, [t; u; residuals] in RotatedSecondOrderCone())
@objective(model, Min, t)
optimize!(model)
value(θ), value(t)

# ## Exponential Cone

# The [`MOI.ExponentialCone`](@ref) is a set of the form:

# ```math
# K_{exp} = \{ (x,y,z) \in \mathbb{R}^3 : y \exp (x/y) \le z, y > 0 \}
# ```

# It can be used to model problems involving `log` and `exp`. For example, the
# entropy maximization problem consists of maximizing the entropy function,
# $H(x) = -x\log{x}$ subject to linear inequality constraints.

# ```math
# \begin{aligned}
# & \max & - \sum_{i=1}^n x_i \log x_i \\
# & \;\;\text{s.t.} & \mathbf{1}^\top x = 1 \\
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

m, n = 10, 15
A, b = randn(m, n), rand(m, 1)
model = Model(SCS.Optimizer)
set_silent(model)
@variable(model, t[1:n])
@variable(model, x[1:n])
@objective(model, Max, sum(t))
@constraint(model, sum(x) == 1)
@constraint(model, A * x .<= b)
@constraint(model, [i = 1:n], [t[i], x[i], 1] in MOI.ExponentialCone())
optimize!(model)
objective_value(model)

# The [`MOI.ExponentialCone`](@ref) has a dual, the [`MOI.DualExponentialCone`](@ref),
# that offers an alternative formulation that can be more efficient for some
# formulations.

# There is also the [`MOI.RelativeEntropyCone`](@ref) for explicitly encoding
# the relative entropy function

model = Model(SCS.Optimizer)
set_silent(model)
@variable(model, t)
@variable(model, x[1:n])
@objective(model, Max, -t)
@constraint(model, sum(x) == 1)
@constraint(model, A * x .<= b)
@constraint(model, [t; ones(n); x] in MOI.RelativeEntropyCone(2n + 1))
optimize!(model)
objective_value(model)

# ## PowerCone

# The [`MOI.PowerCone`](@ref) is a three-dimensional set parameterized by a
# scalar value `α`. It has the form:
# ```math
# K_{p} = \{ (x,y,z) \in \mathbb{R}^3 : x^{\alpha} y^{1-\alpha} \ge |z|, x \ge 0, y \ge 0 \}
# ```

# The power cone permits a number of reformulations. For example, when ``p > 1``,
# we can model ``t \ge x^p`` using the power cone ``(t, 1, x)`` with
# ``\alpha = 1 / p``. Thus, to model ``t \ge x^3`` with ``x \ge 0``

model = Model(SCS.Optimizer)
set_silent(model)
@variable(model, t)
@variable(model, x >= 1.5)
@constraint(model, [t, 1, x] in MOI.PowerCone(1 / 3))
@objective(model, Min, t)
optimize!(model)
value(t), value(x)

# The [`MOI.PowerCone`](@ref) has a dual, the [`MOI.DualPowerCone`](@ref),
# that offers an alternative formulation that can be more efficient for some
# formulations.

# ## P-Norm

# The p-norm ``||x||_p = \left(\sum\limits_{i} |x_i|^p\right)^{\frac{1}{p}}``
# can be modeled using [`MOI.PowerCone`](@ref)s. See the [Mosek Modeling Cookbook](https://docs.mosek.com/modeling-cookbook/powo.html#p-norm-cones)
# for the derivation.

function p_norm(x::Vector, p)
    N = length(x)
    model = Model(SCS.Optimizer)
    set_silent(model)
    @variable(model, r[1:N])
    @variable(model, t)
    @constraint(model, [i = 1:N], [r[i], t, x[i]] in MOI.PowerCone(1 / p))
    @constraint(model, sum(r) == t)
    @objective(model, Min, t)
    optimize!(model)
    return value(t)
end

x = rand(5);
LinearAlgebra.norm(x, 4), p_norm(x, 4)

# ## Positive Semidefinite Cone

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
objective_value(model)

# ## GeometricMeanCone

# The [`MOI.GeometricMeanCone`](@ref) is a cone of the form:
# ```math
# K_{geo} = \{ (t, x) \in \mathbb{R}^n : x \ge 0, t \le \sqrt[n-1]{x_1 x_2 \cdots x_{n-1}} \}
# ```

model = Model(SCS.Optimizer)
set_silent(model)
@variable(model, x[1:4])
@variable(model, t)
@constraint(model, sum(x) == 1)
@constraint(model, [t; x] in MOI.GeometricMeanCone(5))
optimize!(model)
value(t), value.(x)

# ## Other Cones and Functions

# For other cones supported by JuMP, check out the
# [MathOptInterface Manual](https://jump.dev/MathOptInterface.jl/stable).

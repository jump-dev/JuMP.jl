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

# # Conic Programming

# **Originally Contributed by**: Arpit Bhatia

# This tutorial is aimed at providing a simplistic introduction to conic
# programming using JuMP.

# It uses the following packages:

using JuMP
import SCS
import LinearAlgebra
import Random

Random.seed!(1234)

# ## What is a Cone?

# A subset $C$ of a vector space $V$ is a cone if $\forall x \in C$ and positive
# scalars $\alpha$, the product $\alpha x \in C$. A cone C is a convex cone if
# $\alpha x + \beta y \in C$, for any positive scalars $\alpha, \beta$, and any
# $x, y \in C$.

# ## Conic Programming

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

# ## Some of the Types of Cones Supported by JuMP

# ### Second-Order Cone

# The Second-Order Cone (or Lorentz Cone) of dimension $n$ is of the form:
#
# ```math
# Q^n = \{ (t,x) \in \mathbb{R}^\mbox{n} : t \ge ||x||_2 \}
# ```

# A Second-Order Cone rotated by $\pi/4$ in the $(x_1,x_2)$ plane is called a
# Rotated Second-Order Cone. It is of the form:

# ```math
# Q_r^n = \{ (t,u,x) \in \mathbb{R}^\mbox{n} : 2tu \ge ||x||_2^2, t,u \ge 0 \}
# ```

# These cones are represented in JuMP using the MOI sets [`SecondOrderCone`](@ref)
# and [`RotatedSecondOrderCone`](@ref0).

# #### Example: Euclidean Projection on a hyperplane

# For a given point $u_{0}$ and a set $K$, we refer to any point $u \in K$ which
# is closest to $u_{0}$ as a projection of $u_{0}$ on $K$.
#
# The projection of a point $u_{0}$ on a hyperplane $K = \{u | p' \cdot u = q\}$
# is given by:

# ```math
# \begin{aligned}
# & \min & ||u - u_{0}|| \\
# & \;\;\text{s.t.} & p' \cdot u = q  \\
# \end{aligned}
# ```

u0 = rand(10)
p = rand(10)
q = rand();

# We can model the above problem as the following conic program:

# ```math
# \begin{aligned}
# & \min & t \\
# & \;\;\text{s.t.} & p' \cdot u = q \\
# & & (t, u - u_{0}) \in Q^{n+1}
# \end{aligned}
# ```

# On comparing this with the primal form of a conic problem we saw above,

# ```math
# \begin{aligned}
# & x = (t , u) &\\
# & a_0 = e_1 &\\
# & b_0 = 0 &\\
# & A_1 = (0, p) &\\
# & b_1 = -q &\\
# & C_1 = \mathbb{R}_- &\\
# & A_2 = 1 &\\
# & b_2 = -(0, u_0) &\\
# & C_2 = Q^{n+1} &
# \end{aligned}
# ```

# Thus, we can obtain the dual problem as:

# ```math
# \begin{aligned}
# & \max & y_1 + (0, u_0)^T y_2 \\
# & \;\;\text{s.t.} & e_1 - (0,p)^T y_1 - y_2 = 0 \\
# & & y_1 \in \mathbb{R}_- \\
# & & y_2 \in Q^{n+1}
# \end{aligned}
# ```

model = Model(SCS.Optimizer)
set_silent(model)
@variable(model, u[1:10])
@variable(model, t)
@objective(model, Min, t)
@constraint(model, [t, (u - u0)...] in SecondOrderCone())
@constraint(model, u' * p == q)
optimize!(model)

#-

objective_value(model)

#-

value.(u)

#-

e1 = [1, zeros(10)...]
dual_model = Model(SCS.Optimizer)
set_silent(model)
@variable(dual_model, y1 <= 0)
@variable(dual_model, y2[1:11])
@objective(dual_model, Max, q * y1 + vcat(0, u0)' * y2)
@constraint(dual_model, e1 - [0, p...] .* y1 - y2 .== 0)
@constraint(dual_model, y2 in SecondOrderCone())
optimize!(dual_model)

#+

objective_value(dual_model)


# We can also have an equivalent formulation using a Rotated Second-Order Cone:

# ```math
# \begin{aligned}
# & \min & t \\
# & \;\;\text{s.t.} & p' \cdot u = q \\
# & & (t, 1/2, u - u_{0})\in Q_r^{n+2}
# \end{aligned}
# ```

model = Model(SCS.Optimizer)
set_silent(model)
@variable(model, u[1:10])
@variable(model, t)
@objective(model, Min, t)
@constraint(model, [t, 0.5, (u - u0)...] in RotatedSecondOrderCone())
@constraint(model, u' * p == q)
optimize!(model)

#+

value.(u)

# The difference here is that the objective in the case of the Second-Order Cone
# is $||u - u_{0}||_2$, while in the case of a Rotated Second-Order Cone is
# $||u - u_{0}||_2^2$. However, the value of x is the same for both.

# ### Exponential Cone

# An Exponential Cone is a set of the form:

# ```math
# K_{exp} = \{ (x,y,z) \in \mathbb{R}^3 : y \exp (x/y) \le z, y > 0 \}
# ```

# It is represented in JuMP using the MOI set `MOI.ExponentialCone`.

# #### Example: Entropy Maximization

# As the name suggests, the entropy maximization problem consists of maximizing
# the entropy function, $H(x) = -x\log{x}$ subject to linear inequality
# constraints.

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

n = 15;
m = 10;
A = randn(m, n);
b = rand(m, 1);

model = Model(SCS.Optimizer)
set_silent(model)
@variable(model, t[1:n])
@variable(model, x[1:n])
@objective(model, Max, sum(t))
@constraint(model, sum(x) == 1)
@constraint(model, A * x .<= b )
# Cannot use the exponential cone directly in JuMP, hence we use MOI to specify
# the set.
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
# \lambda_{1} = \max t \\
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
# [MathOptInterface Manual](http://jump.dev/MathOptInterface.jl/stable).

# A good resource for learning more about functions which can be modelled using
# cones is the [MOSEK Modeling Cookbook](https://docs.mosek.com/modeling-cookbook/index.html).

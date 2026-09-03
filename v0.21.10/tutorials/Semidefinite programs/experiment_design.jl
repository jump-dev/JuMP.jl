# Copyright (c) 2019 Arpit Bhatia, Chris Coey, and contributors                  #src
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

# # Experiment Design

# **Originally Contributed by**: Arpit Bhatia, Chris Coey

# This tutorial covers experiment design examples (D-optimal, A-optimal, and
# E-optimal) from section 7.5 of the book [Convex Optimization](https://doi.org/10.1017/CBO9780511804441)
# by Boyd and Vandenberghe.

# The tutorial uses the following packages
using JuMP
import SCS
import LinearAlgebra
import Random

# We set a seed so the random numbers are repeatable:
Random.seed!(1234)

# ## Relaxed Experiment Design Problem

# The basic experiment design problem is as follows.

# Given the menu of possible choices for experiments, $v_{1}, \ldots, v_{p}$,
# and the total number $m$ of experiments to be carried out, choose the numbers
# of each type of experiment, i.e., $m_{1}, \ldots, m_{p}$ to make the error
# covariance $E$ small (in  some sense).

# The variables $m_{1}, \ldots, m_{p}$ must, of course, be integers and sum to
# $m$ the given total number of experiments. This leads to the optimization
# problem:

# ```math
# \begin{aligned}
# \min\left(\mathrm{w.r.t.} \mathbf{S}_{+}^{n}\right) & E=\left(\sum_{j=1}^{p} m_{j} v_{j} v_{j}^{T}\right)^{-1} \\
# \text{subject to} & m_{i} \geq 0 \\
# & \sum\limits_{i=1}^p m_{i} = m \\
# & m_{i} \in \mathbb{Z},\quad i=1,\ldots,p
# \end{aligned}
# ```

# The basic experiment design problem can be a hard combinatorial problem when
# $m$, the total number of experiments, is comparable to $n$, since in this case
# the $m_{i}$ are all small integers.

# In the case when $m$ is large compared to $n$, however, a good approximate
# solution can be found by ignoring, or relaxing, the constraint that the
# $m_{i}$ are integers.

# Let $\lambda_{i}=m_{i} / m,$ which is the fraction of the total number of
# experiments for which $a_{j}=v_{i},$ or the relative frequency of experiment
# $i$. We can express the error covariance in terms of $\lambda_{i}$ as:

# ```math
# E=\frac{1}{m}\left(\sum_{i=1}^{p} \lambda_{i} v_{i} v_{i}^{T}\right)^{-1}
# ```

# The vector $\lambda \in \mathbf{R}^{p}$ satisfies
# $\lambda \succeq 0, \mathbf{1}^{T} \lambda=1,$ and also, each
# $\lambda_{i}$ is an integer multiple of $1 / m$. By ignoring this last
# constraint, we arrive at the problem:

# ```math
# \begin{aligned}
# \min\left(\mathrm{w.r.t.} \mathbf{S}_{+}^{n}\right) & E=(1 / m)\left(\sum_{i=1}^{p} \lambda_{i} v_{i} v_{i}^{T}\right)^{-1} \\
# \text{subject to:}\quad & \lambda \succeq 0 \\
# & \mathbf{1}^{T} \lambda=1
# \end{aligned}
# ```

# Several scalarizations have been proposed for the experiment design problem,
# which is a vector optimization problem over the positive semidefinite cone.

q = 4 # dimension of estimate space
p = 8 # number of experimental vectors
nmax = 3 # upper bound on lambda
n = 12

V = randn(q, p)

eye = Matrix{Float64}(LinearAlgebra.I, q, q);

# ## A-optimal design

# In A-optimal experiment design, we minimize tr $E$, the trace of the
# covariance matrix. This objective is simply the mean of the norm of the error
# squared:

# ```math
# \mathbf{E}\|e\|_{2}^{2}=\mathbf{E} \operatorname{tr}\left(e e^{T}\right)=\operatorname{tr} E
# ```

# The A-optimal experiment design problem in SDP form is

# ```math
# \begin{aligned}
# \min & \mathbf{1}^{T} u \\
# \text{subject to} & \left[\begin{aligned}{\sum_{i=1}^{p} \lambda_{i} v_{i} v_{i}^{T}} & {e_{k}} \\ {e_{k}^{T}} & {u_{k}}\end{aligned}\right] \succeq 0, \quad k=1, \ldots, n \\
# & \lambda \succeq 0 \\
# &  \mathbf{1}^{T} \lambda=1
# \end{aligned}
# ```

aOpt = Model(SCS.Optimizer)
set_silent(aOpt)
@variable(aOpt, np[1:p], lower_bound = 0, upper_bound = nmax)
@variable(aOpt, u[1:q], lower_bound = 0)
@constraint(aOpt, sum(np) <= n)
for i in 1:q
    matrix = [
        V*LinearAlgebra.diagm(0 => np ./ n)*V' eye[:, i]
        eye[i, :]' u[i]
    ]
    @SDconstraint(aOpt, matrix >= 0)
end
@objective(aOpt, Min, sum(u))
optimize!(aOpt)
objective_value(aOpt)

#-

value.(np)

# ## E-optimal design

# In $E$ -optimal design, we minimize the norm of the error covariance matrix,
# i.e. the maximum eigenvalue of $E$.

# Since the diameter (twice the longest semi-axis) of the confidence ellipsoid
# $\mathcal{E}$ is proportional to $\|E\|_{2}^{1 / 2}$, minimizing $\|E\|_{2}$
# can be interpreted geometrically as minimizing the diameter of the confidence
# ellipsoid.

# E-optimal design can also be interpreted as minimizing the maximum variance of
# $q^{T} e$, over all $q$ with $\|q\|_{2}=1$. The E-optimal experiment design
# problem in SDP form is:

# ```math
# \begin{aligned}
# \min & t \\
# \text{subject to} & \sum_{i=1}^{p} \lambda_{i} v_{i} v_{i}^{T} \succeq t I \\
# & \lambda \succeq 0 \\
# & \mathbf{1}^{T} \lambda=1
# \end{aligned}
# ```

eOpt = Model(SCS.Optimizer)
set_silent(eOpt)
@variable(eOpt, 0 <= np[1:p] <= nmax)
@variable(eOpt, t)
@SDconstraint(
    eOpt,
    V * LinearAlgebra.diagm(0 => np ./ n) * V' - (t .* eye) >= 0
)
@constraint(eOpt, sum(np) <= n)
@objective(eOpt, Max, t)
optimize!(eOpt)
objective_value(eOpt)
#-
value.(np)

# ## D-optimal design

# The most widely used scalarization is called $D$ -optimal design,
# in which we minimize the determinant of the error covariance matrix $E$.
# This corresponds to designing the experiment to minimize the volume of the
# resulting confidence ellipsoid (for a fixed confidence level). Ignoring the
# constant factor $1 / m$ in $E$, and taking the logarithm of the objective, we
# can pose this problem as convex optimization problem:

# ```math
# \begin{aligned}
# \min & \log \operatorname{det}\left(\sum_{i=1}^{p} \lambda_{i} v_{i} v_{i}^{T}\right)^{-1} \\
# \text{subject to} & \lambda \succeq 0 \\
# & \mathbf{1}^{T} \lambda=1
# \end{aligned}
# ```

dOpt = Model(SCS.Optimizer)
set_silent(dOpt)
@variable(dOpt, np[1:p], lower_bound = 0, upper_bound = nmax)
@variable(dOpt, t)
@objective(dOpt, Max, t)
@constraint(dOpt, sum(np) <= n)
E = V * LinearAlgebra.diagm(0 => np ./ n) * V'
@constraint(
    dOpt,
    [t, 1, (E[i, j] for i in 1:q for j in 1:i)...] in MOI.LogDetConeTriangle(q)
)
optimize!(dOpt)
objective_value(dOpt)
#-
value.(np)

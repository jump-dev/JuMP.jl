# MIT License                                                                    #src
#                                                                                #src
# Copyright (c) 2020 François Pacaud                                             #src
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

# # Logistic Regression

# **Originally Contributed by**: François Pacaud

# This tutorial shows how to solve a logistic regression problem
# with JuMP. Logistic regression is a well known method in machine learning,
# useful when we want to classify binary variables with the help of
# a given set of features. To this goal,
# we find the optimal combination of features maximizing
# the (log)-likelihood onto a training set. From a modern optimization glance,
# the resulting problem is convex and differentiable. On a modern optimization
# glance, it is even conic representable.
#
# ## Formulating the logistic regression problem
#
# Suppose we have a set of training data-point $i = 1, \cdots, n$, where
# for each $i$ we have a vector of features $x_i \in \mathbb{R}^p$ and a
# categorical observation $y_i \in \{-1, 1\}$.
#
# The log-likelihood is given by
#
# ```math
# l(\theta) = \sum_{i=1}^n \log(\dfrac{1}{1 + \exp(-y_i \theta^\top x_i)})
# ```
#
# and the optimal $\theta$ minimizes the logistic loss function:
#
# ```math
# \min_{\theta}\; \sum_{i=1}^n \log(1 + \exp(-y_i \theta^\top x_i)).
# ```
#
# Most of the time, instead of solving directly the previous optimization problem, we
# prefer to add a regularization term:
#
# ```math
# \min_{\theta}\; \sum_{i=1}^n \log(1 + \exp(-y_i \theta^\top x_i)) + \lambda \| \theta \|
# ```
#
# with $\lambda \in \mathbb{R}_+$ a penalty and $\|.\|$ a norm function. By adding
# such a regularization term, we avoid overfitting on the training set and usually
# achieve a greater score in cross-validation.

# ## Reformulation as a conic optimization problem
#
# By introducing auxiliary variables $t_1, \cdots, t_n$ and $r$,
# the optimization problem is equivalent to
#
# ```math
# \begin{aligned}
# \min_{t, r, \theta} \;& \sum_{i=1}^n t_i + \lambda r \\
# \text{subject to } & \quad t_i \geq \log(1 + \exp(- y_i \theta^\top x_i)) \\
#                    & \quad r \geq \|\theta\|
# \end{aligned}
# ```
#
# Now, the trick is to reformulate the constraints $t_i \geq \log(1 + \exp(- y_i \theta^\top x_i))$
# with the help of the *exponential cone*
#
# ```math
# K_{exp} = \{ (x, y, z) \in \mathbb{R}^3 : \; y \exp(x / y) \leq z \} .
# ```
#
# Indeed, by passing to the exponential, we
# see that for all $i=1, \cdots, n$, the constraint $t_i \geq \log(1 + \exp(- y_i \theta^\top x_i))$
# is equivalent to
#
# ```math
# \exp(-t_i) + \exp(u_i - t_i) \leq 1
# ```
#
# with $u_i = -y_i \theta^\top x_i$. Then, by adding two auxiliary variables
# $z_{i1}$ and $z_{i2}$ such that $z_{i1} \geq \exp(u_i-t_i)$ and $z_{i2} \geq \exp(-t_i)$, we get
# the equivalent formulation
#
# ```math
# \left\{
# \begin{aligned}
# (u_i -t_i , 1, z_{i1}) & \in  K_{exp}  \\
# (-t_i , 1, z_{i2}) & \in  K_{exp}  \\
# z_{i1} + z_{i2} & \leq  1
# \end{aligned}
# \right.
# ```
#
# In this setting, the conic version of the logistic regression problems writes out
#
# ```math
# \begin{aligned}
# \min_{t, z, r, \theta}&  \; \sum_{i=1}^n t_i + \lambda r \\
# \text{subject to } & \quad  (u_i -t_i , 1, z_{i1})  \in  K_{exp}  \\
#                    & \quad  (-t_i , 1, z_{i2})  \in  K_{exp}  \\
#                    & \quad  z_{i1} + z_{i2}  \leq  1 \\
#                    & \quad u_i = -y_i x_i^\top \theta \\
#                    & \quad r \geq \|\theta\|
# \end{aligned}
# ```
#
# and thus encompasses $3n + p + 1$ variables and $3n + 1$ constraints ($u_i = -y_i \theta^\top x_i$
# is only a virtual constraint used to clarify the notation).
# Thus, if $n \gg 1$, we get a large number of variables and constraints.

# ## Fitting logistic regression with a conic solver
#
# It is now time to pass to the implementation. We choose SCS as a conic solver.
using JuMP
import Random
import SCS

Random.seed!(2713);

# We start by implementing a function to generate a fake dataset, and where
# we could tune the correlation between the feature variables. The function
# is a direct transcription of the one used in [this blog post](http://fa.bianp.net/blog/2013/numerical-optimizers-for-logistic-regression/).
function generate_dataset(n_samples = 100, n_features = 10; shift = 0.0)
    X = randn(n_samples, n_features)
    w = randn(n_features)
    y = sign.(X * w)
    X .+= 0.8 * randn(n_samples, n_features) # add noise
    X .+= shift # shift the points in the feature space
    X = hcat(X, ones(n_samples, 1))
    return X, y
end

# We write a `softplus` function to formulate each constraint
# $t \geq \log(1 + \exp(u))$ with two exponential cones.
function softplus(model, t, u)
    z = @variable(model, [1:2], lower_bound = 0.0)
    @constraint(model, sum(z) <= 1.0)
    @constraint(model, [u - t, 1, z[1]] in MOI.ExponentialCone())
    @constraint(model, [-t, 1, z[2]] in MOI.ExponentialCone())
end

# ### $\ell_2$ regularized logistic regression
# Then, with the help of the `softplus` function, we could write our
# optimization model. In the $\ell_2$ regularization case, the constraint
# $r \geq \|\theta\|_2$ rewrites as a second order cone constraint.
function build_logit_model(X, y, λ)
    n, p = size(X)
    model = Model()
    @variable(model, θ[1:p])
    @variable(model, t[1:n])
    for i in 1:n
        u = -(X[i, :]' * θ) * y[i]
        softplus(model, t[i], u)
    end
    ## Add ℓ2 regularization
    @variable(model, 0.0 <= reg)
    @constraint(model, [reg; θ] in MOI.SecondOrderCone(p + 1))
    ## Define objective
    @objective(model, Min, sum(t) + λ * reg)
    return model
end

# We generate the dataset.
#
# !!! warning
#     Be careful here, for large n and p SCS could fail to converge!
#
n, p = 200, 10
X, y = generate_dataset(n, p, shift = 10.0);

## We could now solve the logistic regression problem
λ = 10.0
model = build_logit_model(X, y, λ)
set_optimizer(model, SCS.Optimizer)
set_silent(model)
JuMP.optimize!(model)

#-

θ♯ = JuMP.value.(model[:θ])

# It appears that the speed of convergence is not that impacted by the correlation
# of the dataset, nor by the penalty $\lambda$.

# ### $\ell_1$ regularized logistic regression
#
# We now formulate the logistic problem with a $\ell_1$ regularization term.
# The $\ell_1$ regularization ensures sparsity in the optimal
# solution of the resulting optimization problem. Luckily, the $\ell_1$ norm
# is implemented as a set in `MathOptInterface`. Thus, we could easily formulate
# the sparse logistic regression problem with the help of a `MOI.NormOneCone`
# set.
function build_sparse_logit_model(X, y, λ)
    n, p = size(X)
    model = Model()
    @variable(model, θ[1:p])
    @variable(model, t[1:n])
    for i in 1:n
        u = -(X[i, :]' * θ) * y[i]
        softplus(model, t[i], u)
    end
    ## Add ℓ1 regularization
    @variable(model, 0.0 <= reg)
    @constraint(model, [reg; θ] in MOI.NormOneCone(p + 1))
    ## Define objective
    @objective(model, Min, sum(t) + λ * reg)
    return model
end

## Auxiliary function to count non-null components:
count_nonzero(v::Vector; tol = 1e-6) = sum(abs.(v) .>= tol)

## We solve the sparse logistic regression problem on the same dataset as before.
λ = 10.0
sparse_model = build_sparse_logit_model(X, y, λ)
set_optimizer(sparse_model, SCS.Optimizer)
set_silent(sparse_model)
JuMP.optimize!(sparse_model)

#-

θ♯ = JuMP.value.(sparse_model[:θ])
println(
    "Number of non-zero components: ",
    count_nonzero(θ♯),
    " (out of ",
    p,
    " features)",
)

# ### Extensions
# A direct extension would be to consider the sparse logistic regression with
# *hard* thresholding, which, on contrary to the *soft* version using a $\ell_1$ regularization,
# adds an explicit cardinality constraint in its formulation:
#
# ```math
# \begin{aligned}
# \min_{\theta} & \; \sum_{i=1}^n \log(1 + \exp(-y_i \theta^\top x_i)) + \lambda \| \theta \|_2^2 \\
# \text{subject to } & \quad \| \theta \|_0 <= k
# \end{aligned}
# ```
#
# where $k$ is the maximum number of non-zero components in the vector $\theta$,
# and $\|.\|_0$ is the $\ell_0$ pseudo-norm:
#
# ```math
# \| x\|_0 = \#\{i : \; x_i \neq 0\}
# ```
#
# The cardinality constraint $\|\theta\|_0 \leq k$ could be reformulated with
# binary variables. Thus the hard sparse regression problem could be solved
# by any solver supporting mixed integer conic problems.
#

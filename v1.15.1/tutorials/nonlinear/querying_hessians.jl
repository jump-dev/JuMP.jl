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

# # Computing Hessians

# The purpose of this tutorial is to demonstrate how to compute the Hessian of
# the Lagrangian of a nonlinear program.

# !!! warning
#     This is an advanced tutorial that interacts with the low-level nonlinear
#     interface of MathOptInterface.
#
#     By default, JuMP exports the `MOI` symbol as an alias for the
#     MathOptInterface.jl package. We recommend making this more explicit in
#     your code by adding the following lines:
#     ```julia
#     import MathOptInterface as MOI
#     ```

# Given a nonlinear program:
# ```math
# \begin{align}
# & \min_{x \in \mathbb{R}^n} & f(x) \\
# & \;\;\text{s.t.} & l \le g_i(x) \le u
# \end{align}
# ```
# the Hessian of the Lagrangian is computed as:
# ```math
# H(x, \sigma, \mu) = \sigma\nabla^2 f(x) + \sum_{i=1}^m \mu_i \nabla^2 g_i(x)
# ```
# where ``x`` is a primal point, ``\sigma`` is a scalar (typically ``1``), and
# ``\mu`` is a vector of weights corresponding to the Lagrangian dual of the
# constraints.

# This tutorial uses the following packages:

using JuMP
import Ipopt
import LinearAlgebra
import Random
import SparseArrays

# ## The basic model

# To demonstrate how to interact with the lower-level nonlinear interface, we
# need an example model. The exact model isn't important; we use the model from
# [The Rosenbrock function](@ref) tutorial, with some additional constraints to
# demonstrate various features of the lower-level interface.

model = Model(Ipopt.Optimizer)
set_silent(model)
@variable(model, x[i = 1:2], start = -i)
@constraint(model, g_1, x[1]^2 <= 1)
@constraint(model, g_2, (x[1] + x[2])^2 <= 2)
@objective(model, Min, (1 - x[1])^2 + 100 * (x[2] - x[1]^2)^2)
optimize!(model)

# ## The analytic solution

# With a little work, it is possible to analytically derive the correct hessian:

function analytic_hessian(x, σ, μ)
    g_1_H = [2.0 0.0; 0.0 0.0]
    g_2_H = [2.0 2.0; 2.0 2.0]
    f_H = zeros(2, 2)
    f_H[1, 1] = 2.0 + 1200.0 * x[1]^2 - 400.0 * x[2]
    f_H[1, 2] = f_H[2, 1] = -400.0 * x[1]
    f_H[2, 2] = 200.0
    return σ * f_H + μ' * [g_1_H, g_2_H]
end

# Here are various points:

analytic_hessian([1, 1], 0, [0, 0])

#-

analytic_hessian([1, 1], 0, [1, 0])

#-

analytic_hessian([1, 1], 0, [0, 1])

#-

analytic_hessian([1, 1], 1, [0, 0])

# ## Create a nonlinear model

# JuMP delegates automatic differentiation to the `MOI.Nonlinear` submodule.
# Therefore, to compute the Hessian of the Lagrangian, we need to create a
# [`MOI.Nonlinear.Model`](@ref) object:

rows = Any[]
nlp = MOI.Nonlinear.Model()
for (F, S) in list_of_constraint_types(model)
    if F <: VariableRef
        continue  # Skip variable bounds
    end
    for ci in all_constraints(model, F, S)
        push!(rows, ci)
        object = constraint_object(ci)
        MOI.Nonlinear.add_constraint(nlp, object.func, object.set)
    end
end
MOI.Nonlinear.set_objective(nlp, objective_function(model))
nlp

# It is important that we save the constraint indices in a vector `rows`, so
# that we know the order of the constraints in the nonlinear model.

# Next, we need to convert our model into an [`MOI.Nonlinear.Evaluator`](@ref),
# specifying an automatic differentiation backend. In this case, we use
# [`MOI.Nonlinear.SparseReverseMode`](@ref):

evaluator = MOI.Nonlinear.Evaluator(
    nlp,
    MOI.Nonlinear.SparseReverseMode(),
    index.(all_variables(model)),
)

# Before computing anything with the evaluator, we need to initialize it.
# Use [`MOI.features_available`](@ref) to see what we can query:

MOI.features_available(evaluator)

# Consult the MOI documentation for specifics, but to obtain the Hessian matrix,
# we need to initialize `:Hess`:

MOI.initialize(evaluator, [:Hess])

# MOI represents the Hessian as a sparse matrix. Get the sparsity pattern as
# follows:

hessian_sparsity = MOI.hessian_lagrangian_structure(evaluator)

# The sparsity pattern has a few properties of interest:
# * Each element `(i, j)` indicates a structural non-zero in row `i` and column
#   `j`
# * The list may contain duplicates, in which case we should add the values
#   together
# * The list does not need to be sorted
# * The list may contain any mix of lower- or upper-triangular indices
# This format matches Julia's sparse-triplet form of a SparseArray, so we can
# convert from the sparse Hessian representation to a Julia SparseArray as
# follows:

I = [i for (i, _) in hessian_sparsity]
J = [j for (_, j) in hessian_sparsity]
V = zeros(length(hessian_sparsity))
n = num_variables(model)
H = SparseArrays.sparse(I, J, V, n, n)

# Of course, knowing where the zeros are isn't very interesting. We really want
# to compute the value of the Hessian matrix at a point.

MOI.eval_hessian_lagrangian(evaluator, V, ones(n), 1.0, ones(length(rows)))
H = SparseArrays.sparse(I, J, V, n, n)

# In practice, we often want to compute the value of the hessian at the optimal
# solution.

# First, we compute the primal solution. To do so, we need a vector of the
# variables in the order that they were passed to the solver:

x = all_variables(model)

# Here `x[1]` is the variable that corresponds to column 1, and so on. Here's
# the optimal primal solution:

x_optimal = value.(x)

# Next, we need the optimal dual solution associated with the nonlinear
# constraints (this is where it is important to record the order of the
# constraints as we added them to `nlp`):

y_optimal = dual.(rows)

# Now we can compute the Hessian at the optimal primal-dual point:

MOI.eval_hessian_lagrangian(evaluator, V, x_optimal, 1.0, y_optimal)
H = SparseArrays.sparse(I, J, V, n, n)

# However, this Hessian isn't quite right because it isn't symmetric. We can fix
# this by filling in the appropriate off-diagonal terms:

function fill_off_diagonal(H)
    ret = H + H'
    row_vals = SparseArrays.rowvals(ret)
    non_zeros = SparseArrays.nonzeros(ret)
    for col in 1:size(ret, 2)
        for i in SparseArrays.nzrange(ret, col)
            if col == row_vals[i]
                non_zeros[i] /= 2
            end
        end
    end
    return ret
end

fill_off_diagonal(H)

# Putting everything together:

function compute_optimal_hessian(model::Model)
    rows = Any[]
    nlp = MOI.Nonlinear.Model()
    for (F, S) in list_of_constraint_types(model)
        for ci in all_constraints(model, F, S)
            push!(rows, ci)
            object = constraint_object(ci)
            MOI.Nonlinear.add_constraint(nlp, object.func, object.set)
        end
    end
    MOI.Nonlinear.set_objective(nlp, objective_function(model))
    x = all_variables(model)
    backend = MOI.Nonlinear.SparseReverseMode()
    evaluator = MOI.Nonlinear.Evaluator(nlp, backend, index.(x))
    MOI.initialize(evaluator, [:Hess])
    hessian_sparsity = MOI.hessian_lagrangian_structure(evaluator)
    I = [i for (i, _) in hessian_sparsity]
    J = [j for (_, j) in hessian_sparsity]
    V = zeros(length(hessian_sparsity))
    MOI.eval_hessian_lagrangian(evaluator, V, value.(x), 1.0, dual.(rows))
    H = SparseArrays.sparse(I, J, V, length(x), length(x))
    return Matrix(fill_off_diagonal(H))
end

H_star = compute_optimal_hessian(model)

# If we compare our solution against the analytical solution:

analytic_hessian(value.(x), 1.0, dual.([g_1, g_2]))

# If we look at the eigenvalues of the Hessian:

LinearAlgebra.eigvals(H_star)

# we see that they are all positive. Therefore, the Hessian is positive
# definite, and so the solution found by Ipopt is a local minimizer.

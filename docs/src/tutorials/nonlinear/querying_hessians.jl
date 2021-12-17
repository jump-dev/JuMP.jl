# Copyright (c) 2021 Oscar Dowson and contributors                               #src
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

using JuMP
import Ipopt
import Random
import SparseArrays

n = 1_000
Random.seed!(1234)
data = randn(n)
model = Model(Ipopt.Optimizer)
set_silent(model)
@variable(model, μ, start = 0.0)
@variable(model, σ >= 0.0, start = 1.0)
quad_ref = @constraint(model, μ^2 <= 0.01)
nlp_ref = @NLconstraint(model, (σ + μ)^2 >= 0)
@NLobjective(
    model,
    Max,
    n / 2 * log(1 / (2 * π * σ^2)) -
    sum((data[i] - μ)^2 for i in 1:n) / (2 * σ^2)
)
optimize!(model)

# To access

d = NLPEvaluator(model)

# Before computing anything with the NLPEvaluator, we need to initialize it.
# Use [`MOI.features_available`](@ref) to see what we can query:

MOI.features_available(d)

# Consult the MOI documentation for specifics. But to obtain the Hessian matrix,
# we need to initialize `:Hess`:

MOI.initialize(d, [:Hess])

# MOI represents the Hessian as a sparse matrix. Get the sparsity pattern as
# follows:

hessian_sparsity = MOI.hessian_lagrangian_structure(d)

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

num_g = num_nl_constraints(model)
MOI.eval_hessian_lagrangian(d, V, ones(n), 1.0, ones(num_g))
H = SparseArrays.sparse(I, J, V, n, n)

# In practice, we often want to compute the value of the hessian at the optimal
# solution.

# First, we compute the primal solution. To do so, we need a vector of the
# variables in the order that they were passed to the solver:

x = all_variables(model)

# Here `x[1]` is the variable that corresponds to column 1, and so on. Here's
# the optimal primal solution:

x_optimal = value.(x)

# The next step is a litte trickier.

nlp_cons = all_nl_constraints(model)
y_optimal = dual.(nlp_cons)

# Now we can compute the Hessian at the optimal primal-dual point:

MOI.eval_hessian_lagrangian(d, V, x_optimal, 1.0, y_optimal)
H = SparseArrays.sparse(I, J, V, n, n)

# However, this Hessian only accounts for the objective and constraints entered
# using [`@NLobjective`](@ref) and [`@NLconstraint`](@ref). If we want to take
# quadratic objectives and constraints written using [`@objective`](@ref) or
# [`@constraint`](@ref) into account, we'll need to handle them separately.

# !!! tip
#     If you don't want to do this, you can replace calls to [`@objective`](@ref)
#     and [`@constraint`](@ref) with [`@NLobjective`](@ref) and
#     [`@NLconstraint`](@ref).

# ## Hessians from QuadExpr functions

# To compute the hessian from a quadratic expression, let's see how JuMP
# represents a quadratic constraint:

f = constraint_object(quad_ref).func

# `f`` is a quadratic expression of the form:
# ```
# f(x) = Σqᵢⱼ * xᵢ * xⱼ + Σaᵢ xᵢ + c
# ```
# So `∇²f(x)` is the matrix formed by `[qᵢⱼ]ᵢⱼ`.

variables_to_column = Dict(xi => i for (i, xi) in enumerate(x))

function add_to_hessian(H, f::QuadExpr, μ)
    for (vars, coef) in f.terms
        i = variables_to_column[vars.a]
        j = variables_to_column[vars.b]
        H[i, j] += μ * coef
    end
    return
end

# Then we iterate over all constraints in the model and add their Hessian
# components:

for (F, S) in list_of_constraint_types(model)
    if F <: QuadExpr
        for cref in all_constraints(model, F, S)
            f = constraint_object(cref).func
            add_to_hessian(H, f, dual(cref))
        end
    end
end

H

# Finally, we need to take into account the objective function:

f_obj = objective_function(model)
if f_obj isa QuadExpr
    add_to_hessian(H, f_obj, 1.0)
end

H

# Putting everything together:

"""
    optimal_hessian(model::Model)

Return the Hessian matrix of a nonlinear model, evaluated at the optimal
solution.
"""
function optimal_hessian(model)
    d = NLPEvaluator(model)
    MOI.initialize(d, [:Hess])
    hessian_sparsity = MOI.hessian_lagrangian_structure(d)
    I = [i for (i, _) in hessian_sparsity]
    J = [j for (_, j) in hessian_sparsity]
    V = zeros(length(hessian_sparsity))
    n = num_variables(model)
    H = SparseArrays.sparse(I, J, V, n, n)
    x = all_variables(model)
    x_optimal = value.(x)
    nlp_cons = all_nl_constraints(model)
    y_optimal = dual.(nlp_cons)
    MOI.eval_hessian_lagrangian(d, V, x_optimal, 1.0, y_optimal)
    H = SparseArrays.sparse(I, J, V, n, n)
    variables_to_column = Dict(xi => i for (i, xi) in enumerate(x))
    function add_to_hessian(H, f::QuadExpr, μ)
        for (vars, coef) in f.terms
            i = variables_to_column[vars.a]
            j = variables_to_column[vars.b]
            H[i, j] += μ * coef
        end
        return
    end
    for (F, S) in list_of_constraint_types(model)
        if F <: QuadExpr
            for cref in all_constraints(model, F, S)
                f = constraint_object(cref).func
                add_to_hessian(H, f, dual(cref))
            end
        end
    end
    f_obj = objective_function(model)
    if f_obj isa QuadExpr
        add_to_hessian(H, f_obj, 1.0)
    end
    return H
end

optimal_hessian(model)

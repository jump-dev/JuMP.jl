# Copyright 2017, Iain Dunning, Joey Huchette, Miles Lubin, and contributors    #src
# This Source Code Form is subject to the terms of the Mozilla Public License   #src
# v.2.0. If a copy of the MPL was not distributed with this file, You can       #src
# obtain one at https://mozilla.org/MPL/2.0/.                                   #src

# # Quantum state discrimination

# This tutorial solves the problem of [quantum state discrimination](https://en.wikipedia.org/wiki/Quantum_state_discrimination).

# The purpose of this tutorial to demonstrate how to solve problems involving
# complex-valued decision variables and the [`HermitianPSDCone`(@ref). See
# [Complex number support](@ref) for more details.

# ## Required packages

# This tutorial makes use of the following packages:

using JuMP
import LinearAlgebra
import SCS

# ## formulation

# A `d`-dimensional quantum state, ``\rho``, can be defined by a complex-valued
# Hermitian matrix with a trace of `1`. Assume we have `N` `d`-dimensional
# quantum states, ``\{\rho_i}_{i=1}^n``, each of which is equally likely.

# The goal of the Quantum state discrimination problem is to choose a set of
# positive-operator-valued-measures (POVMs), ``E_i`` such that if we observe
# ``E_i`` then the most probable state that we are in is ``\rho_i``.

# Each POVM ``E_i`` is a complex-valued Hermitian matrix, and there is a
# requirement that ``\sum\limits E_i = \mathbf{I}``.

# To choose the set of POVMs, we want to maximize the probability that we guess
# the quantum state corrrectly. This can be formulated as the following
# optimization problem:

# ```math
# \begin{aligned}
# \max\limits_{E} \;\; & \\mathbb{E}_i[ tr(\rho_i \times E_i)] \\
# \text{s.t.}     \;\; & \sum\limits_i E_i = \mathbf{I} \\
#                      & E_i \succeq 0 \forall i.
# ```

# ## Data

# To setup our problem, we need `N` `d-`dimensional quantum states. To keep the
# problem simple, we use `N = 2` and `d = 2`.

N, d = 2, 2

# We then generated `N` random `d`-dimensional quantum states:

function random_state(d)
    x = randn(ComplexF64, (d, d))
    y = x * x'
    return LinearAlgebra.Hermitian(y / LinearAlgebra.tr(y))
end

ρ = [random_state(d) for i in 1:N]

# ## JuMP formulation

# To model the problem in JuMP, we need a solver that supports positive
# semidefinite matrices:

model = Model(SCS.Optimizer)
set_silent(model)

# Then, we construct our set of `E` variables:

E = [@variable(model, [1:d, 1:d] in HermitianPSDCone()) for i in 1:N]

# Here we have created a vector of matrices. This is different to other modeling
# languages such as YALMIP, which allow you to create a multi-dimensional array
# in which 2-dimensional slices of the array are Hermitian matrices.

# We also need to enforce the constraint that
# ``\sum\limits_i E_i = \mathbf{I}``:

@constraint(model, sum(E) .== LinearAlgebra.I)

# This constraint is a complex-valued equality constraint. In the solver, it
# will be decomposed onto two types of equality constraints: one to enforce
# equality of the real components, and one to enforce equality of the imaginary
# components.

# Our objective is to maximize the expected probability of guessing correctly:

@objective(
    model,
    Max,
    sum(real(LinearAlgebra.tr(ρ[i] * E[i])) for i in 1:N) / N,
)

# Now we optimize:

optimize!(model)
solution_summary(model)

# The POVMs are:

solution = [value.(e) for e in E]

# ## Alternative formulation

# The formulation above includes `n` Hermitian matrices, and a set of linear
# equality constraints. We can simplify the problem by replacing `E[n]` with
# ``I - \sum E_i``, where ``I`` is the identity matrix. This results in:

model = Model(SCS.Optimizer)
set_silent(model)
E = [@variable(model, [1:d, 1:d] in HermitianPSDCone()) for i in 1:N-1]
E_n = LinearAlgebra.Hermitian(LinearAlgebra.I - sum(E))
@constraint(model, E_n in HermitianPSDCone())
push!(E, E_n)

# The objective can also be simplified, by observing that it is equivalent to:

@objective(model, Max, real(LinearAlgebra.dot(ρ, E)) / N)

# Then we can check that we get the same solution:

optimize!(model)
solution_summary(model)

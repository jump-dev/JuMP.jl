# Copyright 2017, Iain Dunning, Joey Huchette, Miles Lubin, and contributors    #src
# This Source Code Form is subject to the terms of the Mozilla Public License   #src
# v.2.0. If a copy of the MPL was not distributed with this file, You can       #src
# obtain one at https://mozilla.org/MPL/2.0/.                                   #src

# # Quantum state discrimination

# This tutorial solves the problem of [quantum state discrimination](https://en.wikipedia.org/wiki/Quantum_state_discrimination).#
# The purpose is to demonstrate how you can solve problems involving
# complex-valued decision variables and the [`HermitianPSDCone`(@ref). See
# [Complex number support](@ref) for more details.

# ## Required packages

# This tutorial makes use of the following packages:

using JuMP
import LinearAlgebra
import SCS

# ## Data

function random_state(d)
    x = randn(ComplexF64, (d, d))
    y = x * x'
    return LinearAlgebra.Hermitian(round.(y / LinearAlgebra.tr(y); digits = 3))
end

N, d = 2, 2

states = [random_state(d) for i in 1:N]

# ## JuMP formulation

model = Model(SCS.Optimizer)
set_silent(model)
E = [@variable(model, [1:d, 1:d] in HermitianPSDCone()) for i in 1:N]
@constraint(model, sum(E) .== LinearAlgebra.I)
@objective(model, Max, real(LinearAlgebra.dot(states, E)) / N)
optimize!(model)
objective_value(model)
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
@objective(model, Max, real(LinearAlgebra.dot(states, E)) / N)
optimize!(model)
objective_value(model)

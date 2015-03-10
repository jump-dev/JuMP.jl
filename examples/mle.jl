#  Copyright 2015, Iain Dunning, Joey Huchette, Miles Lubin, and contributors
#  This Source Code Form is subject to the terms of the Mozilla Public
#  License, v. 2.0. If a copy of the MPL was not distributed with this
#  file, You can obtain one at http://mozilla.org/MPL/2.0/.
using JuMP

# Use nonlinear optimization to compute the maximum likelihood estimate (MLE)
# of the parameters of a normal distribution
# aka the sample mean and variance

n = 1000
data = rand(n)

m = Model()

@defVar(m, μ, start = 0.0)
@defVar(m, σ >= 0.0, start = 1.0)

@setNLObjective(m, Max, (n/2)*log(1/(2π*σ^2))-sum{(data[i]-μ)^2, i=1:n}/(2σ^2))

solve(m)

println("μ = ", getValue(μ))
println("mean(data) = ", mean(data))
println("σ^2 = ", getValue(σ)^2)
println("var(data) = ", var(data))
println("MLE objective: ", getObjectiveValue(m))

# constrained MLE?
@addNLConstraint(m, μ == σ^2)

solve(m)
println("\nWith constraint μ == σ^2:")
println("μ = ", getValue(μ))
println("σ^2 = ", getValue(σ)^2)

println("Constrained MLE objective: ", getObjectiveValue(m))




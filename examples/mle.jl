#  Copyright 2017, Iain Dunning, Joey Huchette, Miles Lubin, and contributors
#  This Source Code Form is subject to the terms of the Mozilla Public
#  License, v. 2.0. If a copy of the MPL was not distributed with this
#  file, You can obtain one at http://mozilla.org/MPL/2.0/.
using JuMP, Ipopt, Statistics

# Use nonlinear optimization to compute the maximum likelihood estimate (MLE)
# of the parameters of a normal distribution
# aka the sample mean and variance

n = 1000
data = randn(n)

m = Model(with_optimizer(Ipopt.Optimizer, print_level=0))

@variable(m, μ, start = 0.0)
@variable(m, σ >= 0.0, start = 1.0)

@NLobjective(m, Max, (n/2)*log(1/(2π*σ^2))-sum((data[i]-μ)^2 for i=1:n)/(2σ^2))

JuMP.optimize!(m)

println("μ = ", JuMP.result_value(μ))
println("mean(data) = ", mean(data))
println("σ^2 = ", JuMP.result_value(σ)^2)
println("var(data) = ", var(data))
println("MLE objective: ", JuMP.objective_value(m))

# constrained MLE?
@NLconstraint(m, μ == σ^2)

JuMP.optimize!(m)
println("\nWith constraint μ == σ^2:")
println("μ = ", JuMP.result_value(μ))
println("σ^2 = ", JuMP.result_value(σ)^2)

println("Constrained MLE objective: ", JuMP.objective_value(m))

#  Copyright 2017, Iain Dunning, Joey Huchette, Miles Lubin, and contributors
#  This Source Code Form is subject to the terms of the Mozilla Public
#  License, v. 2.0. If a copy of the MPL was not distributed with this
#  file, You can obtain one at http://mozilla.org/MPL/2.0/.
#############################################################################
# JuMP
# An algebraic modeling langauge for Julia
# See http://github.com/JuliaOpt/JuMP.jl
#############################################################################
# test/snoop_commands.jl
# Illustrative uses of JuMP for compilation snooping. Reproduced from model.jl.
# Should be referenced in tests in compilation.jl, not included directly in runtests.jl.
#############################################################################
using JuMP, Compat

# If solvers not loaded, load them (i.e running just these tests)
!isdefined(@__MODULE__, :lp_solvers) && include("solvers.jl")

# Create and solve MILP models
for solver in ip_solvers
    modA = Model(solver=solver)
    @variable(modA, x >= 0)
    @variable(modA, y <= 5, Int)
    @variable(modA, 2 <= z <= 4)
    @variable(modA, 0 <= r[i=3:6] <= i)
    @objective(modA, Max, ((x + y)/2.0 + 3.0)/3.0 + z + r[3])
    @constraint(modA, 2 <= x+y)
    @constraint(modA, sum(r[i] for i=3:5) <= (2 - x)/2.0)
    @constraint(modA, 7.0*y <= z + r[6]/1.9)
    solve(modA)
end

# Create and solve LP models (minimization)
for solver in lp_solvers
    modA = Model(solver=solver)
    @variable(modA, x >= 0)
    @variable(modA, y <= 5)
    @variable(modA, 2 <= z <= 4)
    @variable(modA, 0 <= r[i=3:6] <= i)
    @objective(modA, Min, -((x + y)/2.0 + 3.0)/3.0 - z - r[3])
    @constraintref cons[1:3]
    cons[1] = @constraint(modA, x+y >= 2)
    cons[2] = @constraint(modA, sum(r[i] for i=3:5) <= (2 - x)/2.0)
    cons[3] = @constraint(modA, 7.0*y <= z + r[6]/1.9)

    # Solution
    solve(modA)
end

# Create and solve LP models (maxmization)
for solver in lp_solvers
    modA = Model(solver=solver)
    @variable(modA, x >= 0)
    @variable(modA, y <= 5)
    @variable(modA, 2 <= z <= 4)
    @variable(modA, 0 <= r[i=3:6] <= i)
    @objective(modA, Max, ((x + y)/2.0 + 3.0)/3.0 + z + r[3])
    @constraintref cons[1:3]
    cons[1] = @constraint(modA, x+y >= 2)
    cons[2] = @constraint(modA, sum(r[i] for i=3:5) <= (2 - x)/2.0)
    cons[3] = @constraint(modA, 7.0*y <= z + r[6]/1.9)

    # Solution
    solve(modA)
end

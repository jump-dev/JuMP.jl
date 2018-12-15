#############################################################################
#  Copyright 2017, Iain Dunning, Joey Huchette, Miles Lubin, and contributors
#  This Source Code Form is subject to the terms of the Mozilla Public
#  License, v. 2.0. If a copy of the MPL was not distributed with this
#  file, You can obtain one at http://mozilla.org/MPL/2.0/.
#
# Author: Louis Luangkesorn
# Date: Jan 30, 2015
#
# Allocation of passenger cars to trains to minimize cars required or car-miles run.
# Based on
#   Fourer, D.M. Gay and Brian W. Kernighan, A Modeling Language for Mathematical
#   Programming, http://www.ampl.com/REFS/amplmod.ps.gz
#   Appendix D
#############################################################################
using JuMP, Clp, Printf
const MOI = JuMP.MathOptInterface

solver = Clp.Optimizer

ORIG = ["GARY", "CLEV", "PITT"]
DEST = ["FRA", "DET", "LAN", "WIN", "STL", "FRE", "LAF"]

supply = [1400 2600 2900]
demand = [ 900 1200  600 400 1700 1100 1000]

@assert sum(supply) == sum(demand)

cost = [
39   14   11   14   16   82    8;
27    9   12    9   26   95   17;
24   14   17   13   28   99   20
]

m = Model(with_optimizer(solver))

@variable(m, Trans[i in 1:length(ORIG), j in 1:length(DEST)] >= 0)


@objective(m, Min, sum(cost[i,j] * Trans[i,j] for i in 1:length(ORIG), j in 1:length(DEST)))

@constraint(m, [i in 1:length(ORIG)], sum(Trans[i,j] for j in 1:length(DEST)) == supply[i])

@constraint(m, [j in 1:length(DEST)], sum(Trans[i,j] for i in 1:length(ORIG)) == demand[j])

println("Solving original problem...")
JuMP.optimize!(m)

term_status = JuMP.termination_status(m)
primal_status = JuMP.primal_status(m)
is_optimal = term_status == MOI.Optimal

if is_optimal
	@printf("Optimal!\n")
	@printf("Objective value: %d\n", JuMP.objective_value(m))
	@printf("Transpotation:\n")
	for j in 1:length(DEST)
		@printf("\t%s", DEST[j])
	end
	@printf("\n")
	for i in 1:length(ORIG)
		@printf("%s", ORIG[i])
		for j in 1:length(DEST)
			@printf("\t%d", JuMP.value(Trans[i,j]))
		end
		@printf("\n")
	end
else
    @printf("No solution\n")
end

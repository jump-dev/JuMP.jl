#  Copyright 2017, Iain Dunning, Joey Huchette, Miles Lubin, and contributors
#  This Source Code Form is subject to the terms of the Mozilla Public
#  License, v. 2.0. If a copy of the MPL was not distributed with this
#  file, You can obtain one at http://mozilla.org/MPL/2.0/.
#############################################################################
# JuMP
# An algebraic modeling language for Julia
# See http://github.com/JuliaOpt/JuMP.jl
#############################################################################

using JuMP, GLPK, Test
const MOI = JuMP.MathOptInterface

"""
	example_transp()

Allocation of passenger cars to trains to minimize cars required or car-miles
run. Based on:  Fourer, D.M. Gay and Brian W. Kernighan, A Modeling Language for
Mathematical Programming, http://www.ampl.com/REFS/amplmod.ps.gz Appendix D.

Author: Louis Luangkesorn
Date: Jan 30, 2015
"""
function example_transp()
	ORIG = ["GARY", "CLEV", "PITT"]
	DEST = ["FRA", "DET", "LAN", "WIN", "STL", "FRE", "LAF"]

	supply = [1_400, 2_600, 2_900]
	demand = [900, 1_200, 600, 400, 1_700, 1_100, 1_000]

	@assert sum(supply) == sum(demand)

	cost = [
		39   14   11   14   16   82    8;
		27    9   12    9   26   95   17;
		24   14   17   13   28   99   20
	]

	model = Model(GLPK.Optimizer)

	@variable(model, trans[1:length(ORIG), 1:length(DEST)] >= 0)
	@objective(model, Min, sum(cost[i, j] * trans[i, j] for i in 1:length(ORIG), j in 1:length(DEST)))

	@constraint(model, [i in 1:length(ORIG)],
		sum(trans[i, j] for j in 1:length(DEST)) == supply[i])
	@constraint(model, [j in 1:length(DEST)],
		sum(trans[i, j] for i in 1:length(ORIG)) == demand[j])

	JuMP.optimize!(model)

	@test JuMP.termination_status(model) == MOI.OPTIMAL
	@test JuMP.primal_status(model) == MOI.FEASIBLE_POINT
	@test JuMP.objective_value(model) == 196200.0
end

example_transp()

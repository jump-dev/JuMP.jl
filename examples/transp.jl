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
using JuMP, Clp

solver = ClpSolver()

ORIG = ["GARY", "CLEV", "PITT"];
DEST = ["FRA", "DET", "LAN", "WIN", "STL", "FRE", "LAF"];

supply = [1400 2600 2900];
demand = [ 900 1200  600 400 1700 1100 1000];

@assert sum(supply) == sum(demand);

cost = [
39   14   11   14   16   82    8;
27    9   12    9   26   95   17;
24   14   17   13   28   99   20
]

m = Model(solver=solver);

@variable(m, Trans[i=1:length(ORIG), j=1:length(DEST)] >= 0);


@objective(m, Min, sum(cost[i,j] * Trans[i,j] for i=1:length(ORIG), j=1:length(DEST)));

@constraint(m, [i=1:1:length(ORIG)], sum(Trans[i,j] for j=1:length(DEST)) == supply[i]);

@constraint(m, [j = 1:length(DEST)], sum(Trans[i,j] for i=1:length(ORIG)) == demand[j]);

println("Solving original problem...")
status = solve(m);

if status == :Optimal
	@printf("Optimal!\n");
	@printf("Objective value: %d\n", getobjectivevalue(m));
	@printf("Transpotation:\n");
	for j = 1:length(DEST)
		@printf("\t%s", DEST[j]);
	end
	@printf("\n");
	for i = 1:length(ORIG)
		@printf("%s", ORIG[i]);
		for j = 1:length(DEST)
			@printf("\t%d", getvalue(Trans[i,j]));
		end
		@printf("\n");
	end
else
    @printf("No solution\n");
end

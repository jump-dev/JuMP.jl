# # The transportation problem

# Allocation of passenger cars to trains to minimize cars required or car-miles
# run. Based on:
#
# Fourer, D.M. Gay and Brian W. Kernighan, A Modeling Language for Mathematical
# Programming, https://ampl.com/REFS/amplmod.ps.gz Appendix D.
#
# Originally contributed by Louis Luangkesorn, January 30, 2015.

using JuMP
import GLPK
import Test

function example_transp()
	ORIG = ["GARY", "CLEV", "PITT"]
	DEST = ["FRA", "DET", "LAN", "WIN", "STL", "FRE", "LAF"]
	supply = [1_400, 2_600, 2_900]
	demand = [900, 1_200, 600, 400, 1_700, 1_100, 1_000]
	Test.@test sum(supply) == sum(demand)
	cost = [
		39   14   11   14   16   82    8;
		27    9   12    9   26   95   17;
		24   14   17   13   28   99   20
	]
	model = Model(GLPK.Optimizer)
	@variable(model, trans[1:length(ORIG), 1:length(DEST)] >= 0)
	@objective(
		model,
		Min,
		sum(
			cost[i, j] * trans[i, j]
			for i in 1:length(ORIG), j in 1:length(DEST)
		)
	)
	@constraints(model, begin
		[i in 1:length(ORIG)], sum(trans[i, :]) == supply[i]
		[j in 1:length(DEST)], sum(trans[:, j]) == demand[j]
	end)
	optimize!(model)
	Test.@test termination_status(model) == MOI.OPTIMAL
	Test.@test primal_status(model) == MOI.FEASIBLE_POINT
	Test.@test objective_value(model) == 196200.0
	return
end

example_transp()

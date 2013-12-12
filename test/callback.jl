using Gurobi

let
	mod = Model(solver=GurobiSolver(LazyConstraints=1))

	@defVar(mod, 0 <= x <= 2, Int)
	@defVar(mod, 0 <= y <= 2, Int)
	@setObjective(mod, Max, y + 0.5x)
	function corners(cb)
	    x_val = getValue(x)
	    y_val = getValue(y)
    	TOL = 1e-6
	    # Check top right
    	if y_val + x_val > 3 + TOL
        	@addLazyConstraint(cb, y + x <= 3)
    	end
    end
	setlazycallback(mod, corners)
	solve(mod)
	@test_approx_eq_eps getValue(x) 1.0 1e-6
	@test_approx_eq_eps getValue(y) 2.0 1e-6
end


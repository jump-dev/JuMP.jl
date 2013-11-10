m = Model()
@defVar(m, 0 <= x <= 3)
@defVar(m, 1 <= y <= 3)
@setObjective(m, :Max, 1.1x + 1.0y)
maincon = @addConstraint(m, x + y <= 3)
solve(m)
@test_approx_eq_eps getValue(x) 2.0 1e-6
@test_approx_eq_eps getValue(y) 1.0 1e-6

# Test adding a variable
@defVar(m, 0 <= z <= 5, 100.0, [maincon], [1.0])
@test !m.firstsolve
solve(m)
@test_approx_eq_eps getValue(x) 0.0 1e-6
@test_approx_eq_eps getValue(y) 1.0 1e-6
@test_approx_eq_eps getValue(z) 2.0 1e-6
@test_approx_eq_eps getObjectiveValue(m) 201.0 1e-6

# Test changing bounds
setLower(y, 0.0)
setUpper(z, 2.0)
solve(m)
@test_approx_eq_eps getValue(x) 1.0 1e-6
@test_approx_eq_eps getValue(y) 0.0 1e-6
@test_approx_eq_eps getValue(z) 2.0 1e-6
m.firstsolve = true

# Test changing problem type
setUpper(z, 1.5)
m.colCat[3] = JuMP.INTEGER
solve(m)
@test_approx_eq_eps getValue(x) 2.0 1e-6
@test_approx_eq_eps getValue(y) 0.0 1e-6
@test_approx_eq_eps getValue(z) 1.0 1e-6

# Test changing constraint bound
chgConstrRHS(maincon, 2.0)
solve(m)
@test_approx_eq_eps getValue(x) 1.0 1e-6
@test_approx_eq_eps getValue(y) 0.0 1e-6
@test_approx_eq_eps getValue(z) 1.0 1e-6

using JuMP, FactCheck
using Compat

!isdefined(:conic_solvers_with_duals) && include("solvers.jl")

const TOL = 1e-4

facts("[socduals] SOC dual") do
for solver in conic_solvers_with_duals
context("With solver $(typeof(solver))") do
    m = Model(solver=solver)

    @variable(m, x[1:5])
    @objective(m, Max, x[1] + x[2] + x[3] + x[4] + 2*x[5])
    @constraint(m, constrcon1, norm(x[2:5]) <= x[1])
    @constraint(m, constrlin1, x[1] <= 5)
    @constraint(m, constrlin2, x[1]+2*x[2] - x[3] <= 10)
    @constraint(m, constrcon2, norm([2 3;1 1]*x[2:3]-[3;4]) <= x[5] - 2)

    solve(m)
    @fact length(getdual(constrcon1)) --> 5
    @fact length(getdual(constrcon2)) --> 3
    @fact length(m.conicconstrDuals) --> 10
    @fact length(m.linconstrDuals) --> 2

    @fact dot(getdual(constrcon2),[getvalue(x[5]) - 2;[2 3;1 1]*[getvalue(x[2]);getvalue(x[3])]-[3;4]]) --> less_than_or_equal(TOL)

end
end
end

# Test for consistency between LP duals and SOC duals

facts("[socduals] LP dual vs SOC dual (Max)") do
for lp_solver in lp_solvers
context("With LP solver $(typeof(lp_solver))") do

    m1 = Model(solver=lp_solver)
    @variable(m1, x1[1:2] >= 0)
    @objective(m1, Max, x1[1] + 2x1[2])
    @constraint(m1, c11, 3x1[1] + x1[2] <= 4)
    @constraint(m1, c12, x1[1] + 2x1[2] >= 1)
    @constraint(m1, c13, -x1[1] + x1[2] == 0.5)

    solve(m1)

    @fact getdual(c11) --> roughly(0.75, TOL)
    @fact getdual(c12) --> roughly(0.0,TOL)
    @fact getdual(c13) --> roughly(1.25,TOL)

end
end
for  conic_solver in conic_solvers_with_duals
context("With conic solver $(typeof(conic_solver))") do

    m2 = Model(solver=conic_solver)
    @variable(m2, x2[1:2] >= 0)
    @objective(m2, Max, x2[1] + 2x2[2])
    @constraint(m2, c21, 3x2[1] + x2[2] <= 4)
    @constraint(m2, c22, x2[1] + 2x2[2] >= 1)
    @constraint(m2, c23, -x2[1] + x2[2] == 0.5)
    @constraint(m2, norm(x2[1]) <= x2[2])

    solve(m2)

    @fact getdual(c21) --> roughly(0.75, TOL)
    @fact getdual(c22) --> roughly(0.0,TOL)
    @fact getdual(c23) --> roughly(1.25,TOL)

end
end
end

facts("[socduals] LP vs SOC dual (Min)") do
for lp_solver in lp_solvers
context("With LP solver $(typeof(lp_solver))") do

    m1 = Model(solver=lp_solver)
    @variable(m1, x1[1:2] >= 0)
    @objective(m1, Min, -x1[1] - 2x1[2])
    @constraint(m1, c11, 3x1[1] + x1[2] <= 4)
    @constraint(m1, c12, x1[1] + 2x1[2] >= 1)
    @constraint(m1, c13, -x1[1] + x1[2] == 0.5)

    solve(m1)

    @fact getdual(c11) --> roughly(-0.75, TOL)
    @fact getdual(c12) --> roughly(-0.0,TOL)
    @fact getdual(c13) --> roughly(-1.25,TOL)

end
end
for  conic_solver in conic_solvers_with_duals
context("With conic solver $(typeof(conic_solver))") do

    m2 = Model(solver=conic_solver)
    @variable(m2, x2[1:2] >= 0)
    @objective(m2, Min, -x2[1] - 2x2[2])
    @constraint(m2, c21, 3x2[1] + x2[2] <= 4)
    @constraint(m2, c22, x2[1] + 2x2[2] >= 1)
    @constraint(m2, c23, -x2[1] + x2[2] == 0.5)
    @constraint(m2, norm(x2[1]) <= x2[2])

    solve(m2)
    @fact getdual(c21) --> roughly(-0.75, TOL)
    @fact getdual(c22) --> roughly(-0.0,TOL)
    @fact getdual(c23) --> roughly(-1.25,TOL)

end
end
end

facts("[socduals] LP vs SOC reduced costs") do
for lp_solver in lp_solvers
context("With LP solver $(typeof(lp_solver))") do

    m1 = Model(solver=lp_solver)

    @variable(m1, x1 >= 0)
    @variable(m1, y1 >= 0)
    @constraint(m1, x1 + y1 == 1)
    @objective(m1, Max, y1)

    solve(m1)

    @fact dot([getdual(y1), getdual(x1)],[getvalue(y1); getvalue(x1)]) --> less_than_or_equal(TOL)

    @fact getvalue(x1) --> roughly(0.0,TOL)
    @fact getvalue(y1) --> roughly(1.0,TOL)

    m1 = Model(solver=lp_solver)

    @variable(m1, x1 >= 0)
    @constraint(m1, x1 <= 1)
    @objective(m1, Max, x1)

    solve(m1)

    @fact getdual(x1) --> roughly(0.0, TOL)
    @fact getvalue(x1) --> roughly(1.0,TOL)

    RHS = 2
    X_LB = 0.5
    Y_UB = 5

    m1 = Model(solver=lp_solver)

    @variable(m1, x1 >= X_LB)
    @variable(m1, y1 <= Y_UB)
    @constraint(m1, c1,  x1 + y1 == RHS)
    @objective(m1, Max, 0.8*y1 + 0.1*x1)

    solve(m1)

    @fact getvalue(x1) --> roughly(0.5, TOL)
    @fact getvalue(y1) --> roughly(1.5, TOL)

    @fact getdual(x1) --> roughly(-0.7, TOL)
    @fact getdual(y1) --> roughly(0.0, TOL)

    lp_dual_obj = getdual(c1)*RHS + getdual(x1) * X_LB + getdual(y1) * Y_UB
    @fact getobjectivevalue(m1) --> roughly(lp_dual_obj, TOL)

end
end
for  conic_solver in conic_solvers_with_duals
context("With conic solver $(typeof(conic_solver))") do

    m2 = Model(solver=conic_solver)

    @variable(m2, x2 >= 0)
    @variable(m2, y2 >= 0)
    @constraint(m2, x2 + y2 == 1)
    @objective(m2, Max, y2)
    @constraint(m2, norm(x2) <= y2)

    solve(m2)

    @fact dot([getdual(y2), getdual(x2)],[getvalue(y2); getvalue(x2)]) --> less_than_or_equal(TOL)

    @fact getvalue(x2) --> roughly(0.0,TOL)
    @fact getvalue(y2) --> roughly(1.0,TOL)

    m2 = Model(solver=conic_solver)

    @variable(m2, x2 >= 0)
    @constraint(m2, x2 <= 1)
    @objective(m2, Max, x2)
    @constraint(m2, norm(x2) <= x2)

    solve(m2)

    @fact getdual(x2) --> roughly(0.0, TOL)
    @fact getvalue(x2) --> roughly(1.0,TOL)

    RHS = 2
    X_LB = 0.5
    Y_UB = 5

    m2 = Model(solver=conic_solver)

    @variable(m2, x2 >= X_LB)
    @variable(m2, y2 <= Y_UB)
    @constraint(m2, c2, x2 + y2 == RHS)
    @objective(m2, Max, 0.8*y2 + 0.1*x2)
    @constraint(m2, norm(x2) <= y2)

    solve(m2)

    @fact getvalue(x2) --> roughly(0.5, TOL)
    @fact getvalue(y2) --> roughly(1.5, TOL)

    @fact getdual(x2) --> roughly(-0.7, TOL)
    @fact getdual(y2) --> roughly(0.0, TOL)

    conic_dual_obj = getdual(c2)*RHS + getdual(x2) * X_LB + getdual(y2) * Y_UB
    @fact getobjectivevalue(m2) --> roughly(conic_dual_obj, TOL)

end
end
end

# || b-Ax || <= a - c^T x, x\in K
# (a - c^Tx, b-Ax) \in SOC, x\in K
# ((a,b) - (c^T,A)x) \in SOC, x\in K
# if the problem is infeasible, we get a y satisfying -(a,b)^Ty > 0, (c,A^T)y \in K^*, y \in SOC^*

facts("[socduals] SOC infeasibility ray test") do
for  conic_solver in conic_solvers_with_duals
context("With conic solver $(typeof(conic_solver))") do

    m2 = Model(solver=conic_solver)

    @variable(m2, x2 >= 0)
    @constraint(m2, x2 <= 1)
    @objective(m2, Max, x2)
    @constraint(m2, c2, norm(x2) <= x2 - 1)

    status = solve(m2)

    inf_ray = getdual(c2)
    @fact status --> :Infeasible
    @fact inf_ray[1] ≥ abs(inf_ray[2]) - TOL --> true
    @fact -(-inf_ray[1]) --> greater_than_or_equal(TOL)

end
end
end

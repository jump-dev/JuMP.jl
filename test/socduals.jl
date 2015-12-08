using JuMP, FactCheck
using Compat

!isdefined(:conic_solvers_with_duals) && include("solvers.jl")

const TOL = 1e-4

facts("[socduals] SOC dual") do
for solver in conic_solvers_with_duals
context("With solver $(typeof(solver))") do
    m = Model(solver=solver)

    @defVar(m, x[1:5])
    @setObjective(m, Max, x[1] + x[2] + x[3] + x[4] + 2*x[5])
    @addConstraint(m, constrcon1, norm(x[2:5]) <= x[1])
    @addConstraint(m, constrlin1, x[1] <= 5)
    @addConstraint(m, constrlin2, x[1]+2*x[2] - x[3] <= 10)
    @addConstraint(m, constrcon2, norm([2 3;1 1]*x[2:3]-[3;4]) <= x[5] - 2)

    solve(m)
    @fact length(getDual(constrcon1)) --> 5
    @fact length(getDual(constrcon2)) --> 3
    @fact length(m.conicconstrDuals) --> 10
    @fact length(m.linconstrDuals) --> 2

    @fact dot(getDual(constrcon2),[getValue(x[5]) - 2;[2 3;1 1]*[getValue(x[2]);getValue(x[3])]-[3;4]]) --> less_than_or_equal(TOL)

end
end
end

# Test for consistency between LP duals and SOC duals

facts("[socduals] LP dual vs SOC dual (Max)") do
for lp_solver in lp_solvers
context("With LP solver $(typeof(lp_solver))") do

    m1 = Model(solver=lp_solver)
    @defVar(m1, x1[1:2] >= 0)
    @setObjective(m1, Max, x1[1] + 2x1[2])
    @addConstraint(m1, c11, 3x1[1] + x1[2] <= 4)
    @addConstraint(m1, c12, x1[1] + 2x1[2] >= 1)
    @addConstraint(m1, c13, -x1[1] + x1[2] == 0.5)

    solve(m1)

    @fact getDual(c11) --> roughly(0.75, TOL)
    @fact getDual(c12) --> roughly(0.0,TOL)
    @fact getDual(c13) --> roughly(1.25,TOL)

end
end
for  conic_solver in conic_solvers_with_duals
context("With conic solver $(typeof(conic_solver))") do

    m2 = Model(solver=conic_solver)
    @defVar(m2, x2[1:2] >= 0)
    @setObjective(m2, Max, x2[1] + 2x2[2])
    @addConstraint(m2, c21, 3x2[1] + x2[2] <= 4)
    @addConstraint(m2, c22, x2[1] + 2x2[2] >= 1)
    @addConstraint(m2, c23, -x2[1] + x2[2] == 0.5)
    @addConstraint(m2, norm(x2[1]) <= x2[2])

    solve(m2)

    @fact getDual(c21) --> roughly(0.75, TOL)
    @fact getDual(c22) --> roughly(0.0,TOL)
    @fact getDual(c23) --> roughly(1.25,TOL)

end
end
end

facts("[socduals] LP vs SOC dual (Min)") do
for lp_solver in lp_solvers
context("With LP solver $(typeof(lp_solver))") do

    m1 = Model(solver=lp_solver)
    @defVar(m1, x1[1:2] >= 0)
    @setObjective(m1, Min, -x1[1] - 2x1[2])
    @addConstraint(m1, c11, 3x1[1] + x1[2] <= 4)
    @addConstraint(m1, c12, x1[1] + 2x1[2] >= 1)
    @addConstraint(m1, c13, -x1[1] + x1[2] == 0.5)

    solve(m1)

    @fact getDual(c11) --> roughly(-0.75, TOL)
    @fact getDual(c12) --> roughly(-0.0,TOL)
    @fact getDual(c13) --> roughly(-1.25,TOL)

end
end
for  conic_solver in conic_solvers_with_duals
context("With conic solver $(typeof(conic_solver))") do

    m2 = Model(solver=conic_solver)
    @defVar(m2, x2[1:2] >= 0)
    @setObjective(m2, Min, -x2[1] - 2x2[2])
    @addConstraint(m2, c21, 3x2[1] + x2[2] <= 4)
    @addConstraint(m2, c22, x2[1] + 2x2[2] >= 1)
    @addConstraint(m2, c23, -x2[1] + x2[2] == 0.5)
    @addConstraint(m2, norm(x2[1]) <= x2[2])

    solve(m2)
    @fact getDual(c21) --> roughly(-0.75, TOL)
    @fact getDual(c22) --> roughly(-0.0,TOL)
    @fact getDual(c23) --> roughly(-1.25,TOL)

end
end
end

facts("[socduals] LP vs SOC reduced costs") do
for lp_solver in lp_solvers
context("With LP solver $(typeof(lp_solver))") do

    m1 = Model(solver=lp_solver)

    @defVar(m1, x1 >= 0)
    @defVar(m1, y1 >= 0)
    @addConstraint(m1, x1 + y1 == 1)
    @setObjective(m1, Max, y1)

    solve(m1)

    @fact dot([getDual(y1), getDual(x1)],[getValue(y1); getValue(x1)]) --> less_than_or_equal(TOL)

    @fact getValue(x1) --> roughly(0.0,TOL)
    @fact getValue(y1) --> roughly(1.0,TOL)

    m1 = Model(solver=lp_solver)

    @defVar(m1, x1 >= 0)
    @addConstraint(m1, x1 <= 1)
    @setObjective(m1, Max, x1)

    solve(m1)

    @fact getDual(x1) --> roughly(0.0, TOL)
    @fact getValue(x1) --> roughly(1.0,TOL)

    RHS = 2
    X_LB = 0.5
    Y_UB = 5

    m1 = Model(solver=lp_solver)

    @defVar(m1, x1 >= X_LB)
    @defVar(m1, y1 <= Y_UB)
    @addConstraint(m1, c1,  x1 + y1 == RHS)
    @setObjective(m1, Max, 0.8*y1 + 0.1*x1)

    solve(m1)

    @fact getValue(x1) --> roughly(0.5, TOL)
    @fact getValue(y1) --> roughly(1.5, TOL)

    @fact getDual(x1) --> roughly(-0.7, TOL)
    @fact getDual(y1) --> roughly(0.0, TOL)

    lp_dual_obj = getDual(c1)*RHS + getDual(x1) * X_LB + getDual(y1) * Y_UB
    @fact getObjectiveValue(m1) --> roughly(lp_dual_obj, TOL)

end
end
for  conic_solver in conic_solvers_with_duals
context("With conic solver $(typeof(conic_solver))") do

    m2 = Model(solver=conic_solver)

    @defVar(m2, x2 >= 0)
    @defVar(m2, y2 >= 0)
    @addConstraint(m2, x2 + y2 == 1)
    @setObjective(m2, Max, y2)
    @addConstraint(m2, norm(x2) <= y2)

    solve(m2)

    @fact dot([getDual(y2), getDual(x2)],[getValue(y2); getValue(x2)]) --> less_than_or_equal(TOL)

    @fact getValue(x2) --> roughly(0.0,TOL)
    @fact getValue(y2) --> roughly(1.0,TOL)

    m2 = Model(solver=conic_solver)

    @defVar(m2, x2 >= 0)
    @addConstraint(m2, x2 <= 1)
    @setObjective(m2, Max, x2)
    @addConstraint(m2, norm(x2) <= x2)

    solve(m2)

    @fact getDual(x2) --> roughly(0.0, TOL)
    @fact getValue(x2) --> roughly(1.0,TOL)

    RHS = 2
    X_LB = 0.5
    Y_UB = 5

    m2 = Model(solver=conic_solver)

    @defVar(m2, x2 >= X_LB)
    @defVar(m2, y2 <= Y_UB)
    @addConstraint(m2, c2, x2 + y2 == RHS)
    @setObjective(m2, Max, 0.8*y2 + 0.1*x2)
    @addConstraint(m2, norm(x2) <= y2)

    solve(m2)

    @fact getValue(x2) --> roughly(0.5, TOL)
    @fact getValue(y2) --> roughly(1.5, TOL)

    @fact getDual(x2) --> roughly(-0.7, TOL)
    @fact getDual(y2) --> roughly(0.0, TOL)

    conic_dual_obj = getDual(c2)*RHS + getDual(x2) * X_LB + getDual(y2) * Y_UB
    @fact getObjectiveValue(m2) --> roughly(conic_dual_obj, TOL)

end
end
end



facts("[socduals] SOC infeasibility ray test") do
for  conic_solver in conic_solvers_with_duals
context("With conic solver $(typeof(conic_solver))") do

    m2 = Model(solver=conic_solver)

    @defVar(m2, x2 >= 0)
    @addConstraint(m2, x2 <= -1)
    @setObjective(m2, Max, x2)
    @addConstraint(m2, c2, norm(x2) <= x2)

    status = solve(m2)

    inf_ray = getDual(c2)
    @fact status --> :Infeasible
    @fact (-inf_ray[1] - inf_ray[2]) --> less_than_or_equal(-TOL)

end
end
end

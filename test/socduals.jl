using JuMP, FactCheck
using Compat

!isdefined(:conic_solvers_with_duals) && include("solvers.jl")

const TOL = 1e-4

facts("SOC dual vector return test") do
for solver in conic_solvers_with_duals
context("With solver $(typeof(solver))") do
    m = Model(solver=solver)

    @defVar(m, x[1:5])
    @setObjective(m, Max, x[1] + x[2] + x[3] + x[4] + 2*x[5])
    @addConstraint(m, constrcon1, norm(x[2:5]) <= x[1])
    @addConstraint(m, constrlin1, x[1] <= 5)
    @addConstraint(m, constrlin2, x[1]+2*x[2] - x[3] <= 10)
    @addConstraint(m, constrcon2, norm([2 3;1 1]*x[2:3]-[3;4]) <= x[5] - 2)
    #@addConstraint(m, constrcon3, x[5]^2 + x[4]^2 <= 10)
    #@addConstraint(m, constrcon2, norm2{x[i], i = 2:5} <= 10)

    #@show typeof(m)
    solve(m)
    #@show MathProgBase.getconstrduals(m.internalModel)
    @fact length(getDual(constrcon1)) --> 5
    @fact length(getDual(constrcon2)) --> 3
    @fact length(m.conicconstrDuals) --> 10
    @fact length(m.linconstrDuals) --> 2

    @fact (getDual(constrcon2)'*[getValue(x[5]) - 2;[2 3;1 1]*[getValue(x[2]);getValue(x[3])]-[3;4]])[1] --> less_than_or_equal(TOL)

end
end
end

facts("LP dual vs SOC dual test / MAX") do
for _lp_solver in lp_solvers
context("With lp solver $(typeof(_lp_solver))") do 

    m1 = Model(solver=_lp_solver)
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
for  _conic_solver in conic_solvers_with_duals
context("With conic solver $(typeof(_conic_solver))") do

    m2 = Model(solver=_conic_solver)
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

facts("LP vs SOC dual test / MIN") do
for _lp_solver in lp_solvers
context("With lp solver $(typeof(_lp_solver))") do 

    m1 = Model(solver=_lp_solver)
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
for  _conic_solver in conic_solvers_with_duals
context("With conic solver $(typeof(_conic_solver))") do

    m2 = Model(solver=_conic_solver)
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

facts("LP vs SOC reduced costs test") do
for _lp_solver in lp_solvers
context("With lp solver $(typeof(_lp_solver))") do 

    m1 = Model(solver=_lp_solver)
    
    @defVar(m1, x1 >= 0)
    @defVar(m1, y1 >= 0)
    @addConstraint(m1, x1 + y1 == 1)
    @setObjective(m1, Max, y1)

    solve(m1)

    @fact dot([getDual(y1), getDual(x1)],[getValue(y1); getValue(x1)]) --> less_than_or_equal(TOL)

    @fact getValue(x1) --> roughly(0.0,TOL)
    @fact getValue(y1) --> roughly(1.0,TOL)

end
end
for  _conic_solver in conic_solvers_with_duals
context("With conic solver $(typeof(_conic_solver))") do

    m2 = Model(solver=_conic_solver)

    @defVar(m2, x2 >= 0)
    @defVar(m2, y2 >= 0)
    @addConstraint(m2, x2 + y2 == 1)
    @setObjective(m2, Max, y2)
    @addConstraint(m2, norm(x2) <= y2)
    
    solve(m2)

    @fact dot([getDual(y2), getDual(x2)],[getValue(y2); getValue(x2)]) --> less_than_or_equal(TOL)

    @fact getValue(x2) --> roughly(0.0,TOL)
    @fact getValue(y2) --> roughly(1.0,TOL)

end
end
end


facts("LP vs SOC reduced costs test") do
for _lp_solver in lp_solvers
context("With lp solver $(typeof(_lp_solver))") do 

    m1 = Model(solver=_lp_solver)
    
    @defVar(m1, x1 >= 0)
    @addConstraint(m1, x1 <= 1)
    @setObjective(m1, Max, x1)

    solve(m1)

    @fact getDual(x1) --> roughly(0.0, TOL)
    @fact getValue(x1) --> roughly(1.0,TOL)

end
end
for  _conic_solver in conic_solvers_with_duals
context("With conic solver $(typeof(_conic_solver))") do

    m2 = Model(solver=_conic_solver)

    @defVar(m2, x2 >= 0)
    @addConstraint(m2, x2 <= 1)
    @setObjective(m2, Max, x2)
    @addConstraint(m2, norm(x2) <= x2)
    
    solve(m2)

    @fact getDual(x2) --> roughly(0.0, TOL)
    @fact getValue(x2) --> roughly(1.0,TOL)

end
end
end

facts("LP vs SOC reduced costs test") do
for _lp_solver in lp_solvers
context("With lp solver $(typeof(_lp_solver))") do 
for  _conic_solver in conic_solvers_with_duals
context("With conic solver $(typeof(_conic_solver))") do

    RHS = 2
    X_LB = 0.5
    Y_UB = 5

    m1 = Model(solver=_lp_solver)
    
    @defVar(m1, x1 >= X_LB)
    @defVar(m1, y1 <= Y_UB)
    @addConstraint(m1, c1,  x1 + y1 == RHS)
    @setObjective(m1, Max, 0.8*y1 + 0.1*x1)

    solve(m1)

    m2 = Model(solver=_conic_solver)

    @defVar(m2, x2 >= X_LB)
    @defVar(m2, y2 <= Y_UB)
    @addConstraint(m2, c2, x2 + y2 == RHS)
    @setObjective(m2, Max, 0.8*y2 + 0.1*x2)
    @addConstraint(m2, norm(x2) <= y2)
    
    solve(m2)

    @fact getValue(x1) --> roughly(getValue(x2), TOL)
    @fact getValue(y1) --> roughly(getValue(y2), TOL)

    @fact getDual(x1) --> roughly(getDual(x2), TOL)
    @fact getDual(y1) --> roughly(getDual(y2), TOL)

    lp_dual_x = getDual(x1)
    lp_dual_y = getDual(y1)
    @fact getDual(x1) --> roughly(lp_dual_x, TOL)
    @fact getDual(y1) --> roughly(lp_dual_y, TOL)
    lp_dual_obj = getDual(c1)*RHS + lp_dual_x * X_LB + lp_dual_y * Y_UB
    @fact getObjectiveValue(m1) --> roughly(lp_dual_obj, TOL)

    conic_dual_x = m2.conicconstrDuals[2]
    conic_dual_y = m2.conicconstrDuals[3]
    @fact getDual(x2) --> roughly(conic_dual_x, TOL)
    @fact getDual(y2) --> roughly(conic_dual_y, TOL)
    conic_dual_obj = getDual(c2)*RHS + conic_dual_x * X_LB + conic_dual_y * Y_UB
    @fact getObjectiveValue(m1) --> roughly(conic_dual_obj, TOL)
end
end
end
end
end

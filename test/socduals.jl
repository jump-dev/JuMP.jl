using JuMP
using FactCheck
using CPLEX
import ECOS

facts("SOC dual vector return test") do

    m = Model(solver=ECOS.ECOSSolver())

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

end

facts("LP dual vs SOC dual test / MAX") do

    m1 = Model(solver=CplexSolver())
    @defVar(m1, x1[1:2] >= 0)
    @setObjective(m1, Max, x1[1] + 2x1[2])
    @addConstraint(m1, c11, 3x1[1] + x1[2] <= 4)
    @addConstraint(m1, c12, x1[1] + 2x1[2] >= 1)
    @addConstraint(m1, c13, -x1[1] + x1[2] == 0.5)


    m2 = Model(solver=ECOS.ECOSSolver())
    @defVar(m2, x2[1:2] >= 0)
    @setObjective(m2, Max, x2[1] + 2x2[2])
    @addConstraint(m2, c21, 3x2[1] + x2[2] <= 4)
    @addConstraint(m2, c22, x2[1] + 2x2[2] >= 1)
    @addConstraint(m2, c23, -x2[1] + x2[2] == 0.5)
    @addConstraint(m2, norm(x2[1]) <= x2[2])

    solve(m1)
    solve(m2)
    #@show m2.conicconstrDuals

    @fact getDual(c11) --> roughly(getDual(c21),1e-6)
    @fact getDual(c12) --> roughly(getDual(c22),1e-6)
    @fact getDual(c13) --> roughly(getDual(c23),1e-6)

end

facts("LP vs SOC dual test / MIN") do

    m1 = Model(solver=CplexSolver())
    @defVar(m1, x1[1:2] >= 0)
    @setObjective(m1, Min, -x1[1] - 2x1[2])
    @addConstraint(m1, c11, 3x1[1] + x1[2] <= 4)
    @addConstraint(m1, c12, x1[1] + 2x1[2] >= 1)
    @addConstraint(m1, c13, -x1[1] + x1[2] == 0.5)


    m2 = Model(solver=ECOS.ECOSSolver())
    @defVar(m2, x2[1:2] >= 0)
    @setObjective(m2, Min, -x2[1] - 2x2[2])
    @addConstraint(m2, c21, 3x2[1] + x2[2] <= 4)
    @addConstraint(m2, c22, x2[1] + 2x2[2] >= 1)
    @addConstraint(m2, c23, -x2[1] + x2[2] == 0.5)
    @addConstraint(m2, norm(x2[1]) <= x2[2])

    solve(m1)
    solve(m2)
    #@show m2.conicconstrDuals

    @fact getDual(c11) --> roughly(getDual(c21),1e-6)
    @fact getDual(c12) --> roughly(getDual(c22),1e-6)
    @fact getDual(c13) --> roughly(getDual(c23),1e-6)

end

facts("LP vs SOC reduced costs test") do

    m1 = Model(solver=CplexSolver())
    @defVar(m1, x11)
    @defVar(m1, 1 <= x12 <= 3)
    @defVar(m1, 0 >= x13 >= -5)
    @setObjective(m1, Min, -x11 - 2x12 + x13)
    @addConstraint(m1, c11, 3x11 + x12 <= 4)
    @addConstraint(m1, c12, x11 + 2x12 >= 1)
    @addConstraint(m1, c13, -x11 + x12 == 0.5)
    @addConstraint(m1, c14, x11 + x12 >= x13)


    m2 = Model(solver=ECOS.ECOSSolver())
    @defVar(m2, x21)
    @defVar(m2, 1 <= x22 <= 3)
    @defVar(m2, 0 >= x23 >= -5)
    @setObjective(m2, Min, -x21 - 2x22 + x23)
    @addConstraint(m2, c21, 3x21 + x22 <= 4)
    @addConstraint(m2, c22, x21 + 2x22 >= 1)
    @addConstraint(m2, c23, -x21 + x22 == 0.5)
    @addConstraint(m2, c24, x21 + x22 >= x23)
    @addConstraint(m2, norm(x21) <= x22)

    solve(m1)
    solve(m2)
    @show getDual(x11)
    @show getDual(x12)
    @show getDual(x13)
    @show getValue(x23)

    @fact getDual(x11) --> roughly(getDual(x21),1e-6)
    @fact getDual(x12) --> roughly(getDual(x22),1e-6)
    @fact getDual(x13) --> roughly(getDual(x23),1e-6)
end

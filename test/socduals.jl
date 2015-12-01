using JuMP
using FactCheck
using CPLEX
import ECOS

TOL = 1e-6

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

    @fact getDual(c11) --> roughly(-getDual(c21),TOL)
    @fact getDual(c12) --> roughly(-getDual(c22),TOL)
    @fact getDual(c13) --> roughly(-getDual(c23),TOL)

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

    @fact getDual(c11) --> roughly(-getDual(c21),TOL)
    @fact getDual(c12) --> roughly(-getDual(c22),TOL)
    @fact getDual(c13) --> roughly(-getDual(c23),TOL)

end

facts("LP vs SOC reduced costs test") do


    m1 = Model(solver=CplexSolver())
    m2 = Model(solver=ECOS.ECOSSolver())
    
    @defVar(m1, x1 >= 0)
    @defVar(m1, y1 >= 0)
    @addConstraint(m1, x1 + y1 == 1)
    @setObjective(m1, Max, y1)

    @defVar(m2, x2 >= 0)
    @defVar(m2, y2 >= 0)
    @addConstraint(m2, x2 + y2 == 1)
    @setObjective(m2, Max, y2)
    @addConstraint(m2, norm(x2) <= y2)
    
    solve(m1)
    solve(m2)

    @fact dot([getDual(y1), getDual(x1)],[getValue(y1); getValue(x1)]) --> greater_than_or_equal(-TOL)
    @fact dot([getDual(y2), getDual(x2)],[getValue(y2); getValue(x2)]) --> greater_than_or_equal(-TOL)

    @fact getValue(x1) --> roughly(getValue(x2),TOL)
    @fact getValue(y1) --> roughly(getValue(y2),TOL)

    #@fact getDual(x1) --> roughly(-getDual(x2),TOL)
    #@fact getDual(y1) --> roughly(-getDual(y2),TOL)
end

facts("LP vs SOC reduced costs test") do


    m1 = Model(solver=CplexSolver())
    m2 = Model(solver=ECOS.ECOSSolver())
    
    @defVar(m1, x1 >= 0)
    @defVar(m1, y1 <= 5)
    @defVar(m1, 3 >= z1 >= 2)
    @addConstraint(m1, x1 + y1 == 1)
    @addConstraint(m1, y1 + z1 >= 3)
    @setObjective(m1, Max, y1)

    @defVar(m2, x2 >= 0)
    @defVar(m2, y2 <= 5)
    @defVar(m2, 3 >= z2 >= 2)
    @addConstraint(m2, x2 + y2 == 1)
    @addConstraint(m2, y2 + z2 >= 3)
    @setObjective(m2, Max, y2)
    @addConstraint(m2, norm(x2) <= y2)
    @addConstraint(m2, norm(y2) <= z2)
    
    solve(m1)
    solve(m2)

    @fact dot([getDual(y1), getDual(x1)],[getValue(y1); getValue(x1)]) --> greater_than_or_equal(-TOL)
    @fact dot([getDual(y2), getDual(x2)],[getValue(y2); getValue(x2)]) --> greater_than_or_equal(-TOL)


    @fact dot([getDual(y1), getDual(z1)],[getValue(y1); getValue(z1)]) --> greater_than_or_equal(-TOL)
    @fact dot([getDual(y2), getDual(z2)],[getValue(y2); getValue(z2)]) --> greater_than_or_equal(-TOL)

    @fact getValue(x1) --> roughly(getValue(x2),TOL)
    @fact getValue(y1) --> roughly(getValue(y2),TOL)

    #@fact getDual(x1) --> roughly(-getDual(x2),TOL)
    #@fact getDual(y1) --> roughly(-getDual(y2),TOL)
end


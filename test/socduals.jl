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

facts("LP dual vs SOC dual test") do

    m1 = Model(solver=CplexSolver())
    @defVar(m1, x[1:2] >= 0)
    @setObjective(m1, Max, x[1] + 2x[2])
    @addConstraint(m1, c1, 3x[1] + x[2] <= 4)


    m2 = Model(solver=ECOS.ECOSSolver())
    @defVar(m2, x[1:2] >= 0)
    @setObjective(m2, Max, x[1] + 2x[2])
    @addConstraint(m2, c2, 3x[1] + x[2] <= 4)
    @addConstraint(m2, norm(x[1]) <= x[2])

    solve(m1)
    solve(m2)
    #@show m2.conicconstrDuals

    @fact getDual(c1) --> roughly(getDual(c2))

end

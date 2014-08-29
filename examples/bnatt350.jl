#############################################################################
# JuMP
# An algebraic modelling langauge for Julia
# See http://github.com/JuliaOpt/JuMP.jl
#############################################################################
# bnatt350.jl
#
# Solve a nontrivial MIP problem by injecting the optimal solution. The LP
# relaxation has the optimal objective value, so optimality is proven
# immediately. The JuMP model is populated from an .MPS file by using
# MathProgBase to construct another model and then copying the data over to
# a JuMP model. The problem comes from the MIPLIB 2010 benchmark collection.
#############################################################################

using JuMP
using MathProgBase

## Uncomment the following three lines to solve using Gurobi
# using Gurobi
# mod = Model(solver=GurobiSolver(PreCrush=1, Cuts=0, Presolve=0, DisplayInterval=1))
# m_internal = MathProgBase.model(GurobiSolver())

using CPLEX
mod = Model(solver=CplexSolver(CPX_PARAM_PRELINEAR=0, CPX_PARAM_PREIND=0, CPX_PARAM_ADVIND=0, CPX_PARAM_MIPSEARCH=1,CPX_PARAM_MIPCBREDLP=0))
m_internal = MathProgBase.model(CplexSolver())

## Uncomment the following three lines to solve using GLPK
# using GLPKMathProgInterface
# mod = Model(solver=GLPKSolverMIP())
# m_internal = MathProgBase.model(GLPKSolverMIP())

# Load the model from .MPS file
MathProgBase.loadproblem!(m_internal, "data/bnatt350.mps")

# grab MathProgBase data
c = MathProgBase.getobj(m_internal)
A = MathProgBase.getconstrmatrix(m_internal)
m, n = size(A)
xlb = MathProgBase.getvarLB(m_internal)
xub = MathProgBase.getvarUB(m_internal)
l = MathProgBase.getconstrLB(m_internal)
u = MathProgBase.getconstrUB(m_internal)
vtypes = MathProgBase.getvartype(m_internal)

# populate JuMP model with data from internal model
@defVar(mod, x[1:n])
for i in 1:n
    setLower(x[i], xlb[i])
    setUpper(x[i], xub[i])
    (vtypes[i] == 'I' || vtypes[i] == 'B') ? mod.colCat[x[i].col] = :Integer : nothing # change vartype to integer when appropriate
end
At = A' # transpose to get useful row-wise sparse representation
for i in 1:At.n
    @addConstraint( mod, l[i] <= sum{ At.nzval[idx]*x[At.rowval[idx]], idx = At.colptr[i]:(At.colptr[i+1]-1) } <= u[i] )
end
@setObjective(mod, Min, sum{ c[i]*x[i], i=1:n })

function myheuristic(cb)
    fp = open("data/bnatt350.sol", "r")
    readline(fp)
    for i in 1:n
        line = chomp(readline(fp))
        spl = filter(x->!isempty(x), split(line, " "))
        setSolutionValue!(cb, x[i], round(float(spl[2]),4))
    end
    addSolution(cb)
end  # End of callback function

setHeuristicCallback(mod, myheuristic)
stat = solve(mod)
println("Solve status: ", stat)
println("Objective value: ", getObjectiveValue(mod))


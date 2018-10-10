#  Copyright 2017, Iain Dunning, Joey Huchette, Miles Lubin, and contributors
#  This Source Code Form is subject to the terms of the Mozilla Public
#  License, v. 2.0. If a copy of the MPL was not distributed with this
#  file, You can obtain one at http://mozilla.org/MPL/2.0/.
#############################################################################
# JuMP
# An algebraic modeling langauge for Julia
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
# m_internal = MathProgBase.LinearQuadraticModel(GurobiSolver())

using CPLEX
mod = Model(solver=CplexSolver(CPX_PARAM_PRELINEAR=0, CPX_PARAM_PREIND=0, CPX_PARAM_ADVIND=0, CPX_PARAM_MIPSEARCH=1,CPX_PARAM_MIPCBREDLP=0))
m_internal = MathProgBase.LinearQuadraticModel(CplexSolver())

## Uncomment the following three lines to solve using GLPK
# using GLPKMathProgInterface
# mod = Model(solver=GLPKSolverMIP())
# m_internal = MathProgBase.LinearQuadraticModel(GLPKSolverMIP())

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
@variable(mod, x[1:n])
for i in 1:n
    setlowerbound(x[i], xlb[i])
    setupperbound(x[i], xub[i])
    mod.colCat[x[i].col] = vtypes[i]
end
@constraint( mod, l .<= A*x .<= u)
@objective(mod, Min, dot(c,x))

function myheuristic(cb)
    fp = open("data/bnatt350.sol", "r")
    readline(fp)
    for i in 1:n
        line = chomp(readline(fp))
        spl = filter(x->!isempty(x), split(line, " "))
        setsolutionvalue(cb, x[i], round(float(spl[2]),4))
    end
    addsolution(cb)
end  # End of callback function

addheuristiccallback(mod, myheuristic)
stat = solve(mod)
println("Solve status: ", stat)
println("Objective value: ", getobjectivevalue(mod))

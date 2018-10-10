#  Copyright 2017, Iain Dunning, Joey Huchette, Miles Lubin, and contributors
#  This Source Code Form is subject to the terms of the Mozilla Public
#  License, v. 2.0. If a copy of the MPL was not distributed with this
#  file, You can obtain one at http://mozilla.org/MPL/2.0/.
#############################################################################
# JuMP
# An algebraic modeling langauge for Julia
# See http://github.com/JuliaOpt/JuMP.jl
#############################################################################
# noswot.jl
#
# Solve a nontrivial MIP problem with user cuts. The jump model is populated
# from an .MPS file by using MathProgBase to construct another model and
# copying the data over to the JuMP model. The problem comes from the
# MIPLIB 3.0 benchmark collection.
#############################################################################

using JuMP
using MathProgBase

using Gurobi
mod = Model(solver=GurobiSolver(PreCrush=1, Cuts=0, Presolve=0, Heuristics=0.0, DisplayInterval=1))
m_internal = MathProgBase.LinearQuadraticModel(GurobiSolver())

## Uncomment the following three lines to solve using CPLEX
# using CPLEX
# mod = Model(solver=CplexSolver(CPX_PARAM_PRELINEAR=0, CPX_PARAM_PREIND=0, CPX_PARAM_ADVIND=0, CPX_PARAM_MIPSEARCH=1,CPX_PARAM_MIPCBREDLP=0))
# m_internal = MathProgBase.LinearQuadraticModel(CplexSolver())

# Load the model from .MPS file
MathProgBase.loadproblem!(m_internal, "data/noswot.mps")

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
    if vtypes[i] == 'I'
        setcategory(x[i], :Int)
    end
end
@constraint(mod, l .<= A*x .<= u)
@objective(mod, Min, dot(c,x))

function mycutgenerator(cb) # valid cuts
    x_val = getvalue(x)
    println("in callback")
    @usercut(cb, x[62]-x[63] <= 0)
    @usercut(cb, x[63]-x[64] <= 0)
    @usercut(cb, x[64]-x[65] <= 0)
    @usercut(cb, 2.08x[52] + 2.98x[62] +   3.47x[72] + 2.24x[82] +  2.08x[92] + 0.25x[51] + 0.25x[61] + 0.25x[71] + 0.25x[81] + 0.25x[91] <= 20.25)
    @usercut(cb, 2.08x[54] + 2.98x[64] +   3.47x[74] + 2.24x[84] +  2.08x[94] + 0.25x[53] + 0.25x[63] + 0.25x[73] + 0.25x[83] + 0.25x[93] <= 20.25)
    @usercut(cb, 2.08x[56] + 2.98x[66] + 3.4722x[76] + 2.24x[86] +  2.08x[96] + 0.25x[55] + 0.25x[65] + 0.25x[75] + 0.25x[85] + 0.25x[95] <= 20.25)
    @usercut(cb, 2.08x[58] + 2.98x[68] +   3.47x[78] + 2.24x[88] +  2.08x[98] + 0.25x[57] + 0.25x[67] + 0.25x[77] + 0.25x[87] + 0.25x[97] <= 20.25)
    @usercut(cb, 2.08x[60] + 2.98x[70] +   3.47x[80] + 2.24x[90] + 2.08x[100] + 0.25x[59] + 0.25x[69] + 0.25x[79] + 0.25x[89] + 0.25x[99] <= 16.25)
end  # End of callback function

# # Tell JuMP/CPLEX to use our callback function
addcutcallback(mod, mycutgenerator)

stat = solve(mod)
println("Solve status: ", stat)
println("Objective value: ", getobjectivevalue(mod))

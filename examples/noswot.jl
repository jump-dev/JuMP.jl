#############################################################################
# JuMP
# An algebraic modelling langauge for Julia
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
m_internal = MathProgBase.model(GurobiSolver())

## Uncomment the following three lines to solve using CPLEX
# using CPLEX
# mod = Model(solver=CplexSolver(CPX_PARAM_PRELINEAR=0, CPX_PARAM_PREIND=0, CPX_PARAM_ADVIND=0, CPX_PARAM_MIPSEARCH=1,CPX_PARAM_MIPCBREDLP=0))
# m_internal = MathProgBase.model(CplexSolver())

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
@defVar(mod, x[1:n])
for i in 1:n
    setLower(x[i], xlb[i])
    setUpper(x[i], xub[i])
    vtypes[i] == 'I' ? mod.colCat[x[i].col] = :Int : nothing # change vartype to integer when appropriate
end
At = A' # transpose to get useful row-wise sparse representation
for i in 1:At.n
    @addConstraint( mod, l[i] <= sum{ At.nzval[idx]*x[At.rowval[idx]], idx = At.colptr[i]:(At.colptr[i+1]-1) } <= u[i] )
end
@setObjective(mod, Min, sum{ c[i]*x[i], i=1:n })

function mycutgenerator(cb) # valid cuts
    x_val = getValue(x)
    println("in callback")
    @addUserCut(cb, x[62]-x[63] <= 0)
    @addUserCut(cb, x[63]-x[64] <= 0)
    @addUserCut(cb, x[64]-x[65] <= 0)
    @addUserCut(cb, 2.08x[52] + 2.98x[62] +   3.47x[72] + 2.24x[82] +  2.08x[92] + 0.25x[51] + 0.25x[61] + 0.25x[71] + 0.25x[81] + 0.25x[91] <= 20.25)
    @addUserCut(cb, 2.08x[54] + 2.98x[64] +   3.47x[74] + 2.24x[84] +  2.08x[94] + 0.25x[53] + 0.25x[63] + 0.25x[73] + 0.25x[83] + 0.25x[93] <= 20.25)
    @addUserCut(cb, 2.08x[56] + 2.98x[66] + 3.4722x[76] + 2.24x[86] +  2.08x[96] + 0.25x[55] + 0.25x[65] + 0.25x[75] + 0.25x[85] + 0.25x[95] <= 20.25)
    @addUserCut(cb, 2.08x[58] + 2.98x[68] +   3.47x[78] + 2.24x[88] +  2.08x[98] + 0.25x[57] + 0.25x[67] + 0.25x[77] + 0.25x[87] + 0.25x[97] <= 20.25)
    @addUserCut(cb, 2.08x[60] + 2.98x[70] +   3.47x[80] + 2.24x[90] + 2.08x[100] + 0.25x[59] + 0.25x[69] + 0.25x[79] + 0.25x[89] + 0.25x[99] <= 16.25)
end  # End of callback function

# # Tell JuMP/CPLEX to use our callback function
addCutCallback(mod, mycutgenerator)

stat = solve(mod)
println("Solve status: ", stat)
println("Objective value: ", getObjectiveValue(mod))

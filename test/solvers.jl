#  Copyright 2015, Iain Dunning, Joey Huchette, Miles Lubin, and contributors
#  This Source Code Form is subject to the terms of the Mozilla Public
#  License, v. 2.0. If a copy of the MPL was not distributed with this
#  file, You can obtain one at http://mozilla.org/MPL/2.0/.
#############################################################################
# JuMP
# An algebraic modelling langauge for Julia
# See http://github.com/JuliaOpt/JuMP.jl
#############################################################################
# test/solvers.jl
# Detect and load solvers
# Should be run as part of runtests.jl
#############################################################################

function try_import(name::Symbol)
    try
        @eval import $name
        return true
    catch e
        return false
    end
end

# Load available solvers
grb = try_import(:Gurobi)
cpx = try_import(:CPLEX)
mos = try_import(:Mosek)
cbc = try_import(:Cbc)
if cbc; import Clp; end
glp = try_import(:GLPKMathProgInterface)
ipt = try_import(:Ipopt)
nlo = try_import(:NLopt)
kni = try_import(:KNITRO)
eco = try_import(:ECOS)
osl = try_import(:CoinOptServices)
scs = try_import(:SCS)
nlw = try_import(:AmplNLWriter)

# Create solver lists
# LP solvers
lp_solvers = Any[]
grb && push!(lp_solvers, Gurobi.GurobiSolver(OutputFlag=0))
cpx && push!(lp_solvers, CPLEX.CplexSolver(CPX_PARAM_SCRIND=0))
mos && push!(lp_solvers, Mosek.MosekSolver(LOG=0))
cbc && push!(lp_solvers, Clp.ClpSolver())
glp && push!(lp_solvers, GLPKMathProgInterface.GLPKSolverLP())
ipt && push!(lp_solvers, Ipopt.IpoptSolver(print_level=0))
eco && push!(lp_solvers, ECOS.ECOSSolver(verbose=false))
osl && push!(lp_solvers, CoinOptServices.OsilSolver())
scs && push!(lp_solvers, SCS.SCSSolver(eps=1e-6,verbose=0))
# MILP solvers
ip_solvers = Any[]
grb && push!(ip_solvers, Gurobi.GurobiSolver(OutputFlag=0))
cpx && push!(ip_solvers, CPLEX.CplexSolver(CPX_PARAM_SCRIND=0))
mos && push!(ip_solvers, Mosek.MosekSolver(LOG=0))
cbc && push!(ip_solvers, Cbc.CbcSolver(logLevel=0))
glp && push!(ip_solvers, GLPKMathProgInterface.GLPKSolverMIP())
osl && push!(ip_solvers, CoinOptServices.OsilSolver())
# Support semi-continuous, semi-integer solvers
semi_solvers = Any[]
grb && push!(semi_solvers, Gurobi.GurobiSolver(OutputFlag=0))
cpx && push!(semi_solvers, CPLEX.CplexSolver(CPX_PARAM_SCRIND=0))
# SOS solvers
sos_solvers = Any[]
grb && push!(sos_solvers, Gurobi.GurobiSolver(OutputFlag=0))
cpx && push!(sos_solvers, CPLEX.CplexSolver(CPX_PARAM_SCRIND=0))
cbc && push!(sos_solvers, Cbc.CbcSolver())
# Conic solvers with duals
conic_solvers_with_duals = Any[]
eco && push!(conic_solvers_with_duals, ECOS.ECOSSolver(verbose=false))
scs && push!(conic_solvers_with_duals, SCS.SCSSolver(eps=1e-6,verbose=0))
#mos && push!(conic_solvers_with_duals, Mosek.MosekSolver(LOG=0))
# Callback solvers
lazy_solvers, cut_solvers, heur_solvers, info_solvers = Any[], Any[], Any[], Any[]
if grb
    push!(lazy_solvers, Gurobi.GurobiSolver(OutputFlag=0, Presolve=0))
    push!( cut_solvers, Gurobi.GurobiSolver(PreCrush=1, Cuts=0, Presolve=0, Heuristics=0.0, OutputFlag=0))
    push!(heur_solvers, Gurobi.GurobiSolver(Cuts=0, Presolve=0, Heuristics=0.0, OutputFlag=0))
    # push!(heur_solvers, Gurobi.GurobiSolver(Heuristics=0.0, OutputFlag=0))
    push!(info_solvers, Gurobi.GurobiSolver(OutputFlag=0))
end
if cpx
    push!(lazy_solvers, CPLEX.CplexSolver(CPX_PARAM_PRELINEAR=0, CPX_PARAM_PREIND=0, CPX_PARAM_ADVIND=0, CPX_PARAM_MIPSEARCH=1,CPX_PARAM_MIPCBREDLP=0,CPX_PARAM_SCRIND=0))
    push!( cut_solvers, CPLEX.CplexSolver(CPX_PARAM_PRELINEAR=0, CPX_PARAM_PREIND=0, CPX_PARAM_ADVIND=0, CPX_PARAM_MIPSEARCH=1,CPX_PARAM_MIPCBREDLP=0,CPX_PARAM_SCRIND=0))
    push!(heur_solvers, CPLEX.CplexSolver(CPX_PARAM_PRELINEAR=0, CPX_PARAM_PREIND=0, CPX_PARAM_ADVIND=0, CPX_PARAM_MIPSEARCH=1,CPX_PARAM_MIPCBREDLP=0,CPX_PARAM_SCRIND=0))
    push!(info_solvers, CPLEX.CplexSolver(CPX_PARAM_PRELINEAR=0, CPX_PARAM_PREIND=0, CPX_PARAM_ADVIND=0, CPX_PARAM_MIPSEARCH=1,CPX_PARAM_MIPCBREDLP=0,CPX_PARAM_SCRIND=0))
end
if glp
    push!(lazy_solvers, GLPKMathProgInterface.GLPKSolverMIP())
    push!( cut_solvers, GLPKMathProgInterface.GLPKSolverMIP())
    push!(heur_solvers, GLPKMathProgInterface.GLPKSolverMIP())
end
# Quadratic support
quad_solvers = Any[]
grb && push!(quad_solvers, Gurobi.GurobiSolver(QCPDual=1,OutputFlag=0))
cpx && push!(quad_solvers, CPLEX.CplexSolver(CPX_PARAM_SCRIND=0))
mos && push!(quad_solvers, Mosek.MosekSolver(LOG=0))
quad_mip_solvers = copy(quad_solvers)
osl && push!(quad_solvers, CoinOptServices.OsilSolver(CoinOptServices.OSOption("sb","yes",solver="ipopt")))
soc_solvers = copy(quad_solvers)
ipt && push!(quad_solvers, Ipopt.IpoptSolver(print_level=0))
eco && push!(soc_solvers, ECOS.ECOSSolver(verbose=false))
scs && push!(soc_solvers, SCS.SCSSolver(eps=1e-6,verbose=0))
osl && push!(quad_mip_solvers, CoinOptServices.OsilBonminSolver(CoinOptServices.OSOption("sb","yes",category="ipopt")))
osl && push!(quad_mip_solvers, CoinOptServices.OsilCouenneSolver())
rsoc_solvers = Any[]
mos && push!(rsoc_solvers, Mosek.MosekSolver(LOG=0))
# Nonlinear solvers
nlp_solvers = Any[]
ipt && push!(nlp_solvers, Ipopt.IpoptSolver(print_level=0))
nlo && push!(nlp_solvers, NLopt.NLoptSolver(algorithm=:LD_SLSQP))
kni && push!(nlp_solvers, KNITRO.KnitroSolver(objrange=1e16,outlev=0))
osl && push!(nlp_solvers, CoinOptServices.OsilSolver(CoinOptServices.OSOption("sb","yes",solver="ipopt")))
nlw && osl && push!(nlp_solvers, AmplNLWriter.BonminNLSolver(["bonmin.nlp_log_level=0"; "bonmin.bb_log_level=0"]))
convex_nlp_solvers = copy(nlp_solvers)
mos && push!(convex_nlp_solvers, Mosek.MosekSolver(LOG=0))
# Mixed-Integer Nonlinear solvers
minlp_solvers = Any[]
kni && push!(minlp_solvers, KNITRO.KnitroSolver(outlev=0))
osl && push!(minlp_solvers, CoinOptServices.OsilBonminSolver(CoinOptServices.OSOption("sb","yes",category="ipopt")))
osl && push!(minlp_solvers, CoinOptServices.OsilCouenneSolver())
nlw && osl && push!(minlp_solvers, AmplNLWriter.BonminNLSolver(["bonmin.nlp_log_level=0"; "bonmin.bb_log_level=0"]))
nlw && osl && push!(minlp_solvers, AmplNLWriter.CouenneNLSolver())
# Semidefinite solvers
sdp_solvers = Any[]
mos && push!(sdp_solvers, Mosek.MosekSolver(LOG=0))
scs && push!(sdp_solvers, SCS.SCSSolver(eps=1e-6,verbose=0))

const error_map = Dict()
grb && (error_map[Gurobi.GurobiSolver] = Gurobi.GurobiError)
cpx && (error_map[CPLEX.CplexSolver] = CPLEX.CplexError)

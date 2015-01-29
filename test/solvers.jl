#############################################################################
# JuMP
# An algebraic modelling langauge for Julia
# See http://github.com/JuliaOpt/JuMP.jl
#############################################################################
# test/solvers.jl
# Detect and load solvers
# Should be run as part of runtests.jl
#############################################################################

# Detect solvers
grb = isdir(Pkg.dir("Gurobi"))
cpx = isdir(Pkg.dir("CPLEX"))
mos = isdir(Pkg.dir("Mosek"))
cbc = isdir(Pkg.dir("Cbc"))
glp = isdir(Pkg.dir("GLPKMathProgInterface"))
ipt = isdir(Pkg.dir("Ipopt"))
nlo = isdir(Pkg.dir("NLopt"))
kni = isdir(Pkg.dir("KNITRO"))
eco = isdir(Pkg.dir("ECOS"))
# Load them
if grb; import Gurobi; end
if cpx; import CPLEX; end
if mos; import Mosek; end
if cbc; import Cbc; import Clp; end
if glp; import GLPKMathProgInterface; end
if ipt; import Ipopt; end
if nlo; import NLopt; end
if kni; import KNITRO; end
if eco; import ECOS; end
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
# MILP solvers
ip_solvers = Any[]
grb && push!(ip_solvers, Gurobi.GurobiSolver(OutputFlag=0))
cpx && push!(ip_solvers, CPLEX.CplexSolver(CPX_PARAM_SCRIND=0))
mos && push!(ip_solvers, Mosek.MosekSolver(LOG=0))
cbc && push!(ip_solvers, Cbc.CbcSolver(logLevel=0))
glp && push!(ip_solvers, GLPKMathProgInterface.GLPKSolverMIP())
# Support semi-continuous, semi-integer solvers
semi_solvers = Any[]
grb && push!(semi_solvers, Gurobi.GurobiSolver(OutputFlag=0))
cpx && push!(semi_solvers, CPLEX.CplexSolver(CPX_PARAM_SCRIND=0))
# SOS solvers
sos_solvers = Any[]
grb && push!(sos_solvers, Gurobi.GurobiSolver(OutputFlag=0))
cpx && push!(sos_solvers, CPLEX.CplexSolver(CPX_PARAM_SCRIND=0))
cbc && push!(sos_solvers, Cbc.CbcSolver())
# Callback solvers
lazy_solvers, cut_solvers, heur_solvers = Any[], Any[], Any[]
if grb
    push!(lazy_solvers, Gurobi.GurobiSolver(OutputFlag=0))
    push!( cut_solvers, Gurobi.GurobiSolver(PreCrush=1, Cuts=0, Presolve=0, Heuristics=0.0, OutputFlag=0))
    push!(heur_solvers, Gurobi.GurobiSolver(Cuts=0, Presolve=0, Heuristics=0.0, OutputFlag=0))
end
if cpx
    push!(lazy_solvers, CPLEX.CplexSolver(CPX_PARAM_PRELINEAR=0, CPX_PARAM_PREIND=0, CPX_PARAM_ADVIND=0, CPX_PARAM_MIPSEARCH=1,CPX_PARAM_MIPCBREDLP=0,CPX_PARAM_SCRIND=0))
    push!( cut_solvers, CPLEX.CplexSolver(CPX_PARAM_PRELINEAR=0, CPX_PARAM_PREIND=0, CPX_PARAM_ADVIND=0, CPX_PARAM_MIPSEARCH=1,CPX_PARAM_MIPCBREDLP=0,CPX_PARAM_SCRIND=0))
    push!(heur_solvers, CPLEX.CplexSolver(CPX_PARAM_PRELINEAR=0, CPX_PARAM_PREIND=0, CPX_PARAM_ADVIND=0, CPX_PARAM_MIPSEARCH=1,CPX_PARAM_MIPCBREDLP=0,CPX_PARAM_SCRIND=0))
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
soc_solvers = copy(quad_solvers)
ipt && push!(quad_solvers, Ipopt.IpoptSolver(print_level=0))
eco && push!(soc_solvers, ECOS.ECOSSolver(verbose=false))
# Nonlinear solvers
nlp_solvers = Any[]
ipt && push!(nlp_solvers, Ipopt.IpoptSolver(print_level=0))
nlo && push!(nlp_solvers, NLopt.NLoptSolver(algorithm=:LD_SLSQP))
kni && push!(nlp_solvers, KNITRO.KnitroSolver(objrange=1e16,outlev=0))
convex_nlp_solvers = copy(nlp_solvers)
mos && push!(convex_nlp_solvers, Mosek.MosekSolver())
# Mixed-Integer Nonlinear solvers
minlp_solvers = Any[]
kni && push!(minlp_solvers, KNITRO.KnitroSolver(outlev=0))

#############################################################################
# JuMP
# An algebraic modelling langauge for Julia
# See http://github.com/JuliaOpt/JuMP.jl
#############################################################################
# test/nonlinear.jl
# To run standlone:
#   julia -L test/nonlinear.jl -e "run_nl_tests(load_nl_solvers())"
#############################################################################

using JuMP
using Base.Test

function test_hs071(nl_solver)
    # hs071
    # Polynomial objective and constraints
    # min x1 * x4 * (x1 + x2 + x3) + x3
    # st  x1 * x2 * x3 * x4 >= 25
    #     x1^2 + x2^2 + x3^2 + x4^2 = 40
    #     1 <= x1, x2, x3, x4 <= 5
    # Start at (1,5,5,1)
    # End at (1.000..., 4.743..., 3.821..., 1.379...)
    m = Model(solver=nl_solver)

    @defVar(m, 1 <= x[1:4] <= 5)

    @setNLObjective(m, Min, x[1]*x[4]*(x[1]+x[2]+x[3]) + x[3])

    @addNLConstraint(m, x[1]*x[2]*x[3]*x[4] >= 25)
    @addNLConstraint(m, sum{x[i]^2,i=1:4} == 40)

    setValue(x[1],1.0)
    setValue(x[2],5.0)
    setValue(x[3],5.0)
    setValue(x[4],1.0)

    status = solve(m)

    @test status == :Optimal
    @test_approx_eq_eps getValue(x[1]) 1.00000000 1e-5
    @test_approx_eq_eps getValue(x[2]) 4.74299963 1e-5
    @test_approx_eq_eps getValue(x[3]) 3.82114998 1e-5
    @test_approx_eq_eps getValue(x[4]) 1.37940829 1e-5
end

function test_hs071_linobj(nl_solver)
    # hs071, with epigraph formulation
    # Polynomial objective and constraints
    # min t
    # st  t >= x1 * x4 * (x1 + x2 + x3) + x3
    #     x1 * x2 * x3 * x4 >= 25
    #     x1^2 + x2^2 + x3^2 + x4^2 = 40
    #     1 <= x1, x2, x3, x4 <= 5
    # Start at (1,5,5,1)
    # End at (1.000..., 4.743..., 3.821..., 1.379...)
    m = Model(solver=nl_solver)
    @defVar(m, 1 <= x[1:4] <= 5)
    @defVar(m, t)
    @setObjective(m, Min, t)
    @addNLConstraint(m, t >= x[1]*x[4]*(x[1]+x[2]+x[3]) + x[3])
    @addNLConstraint(m, x[1]*x[2]*x[3]*x[4] >= 25)
    @addNLConstraint(m, sum{x[i]^2,i=1:4} == 40)

    setValue(x[1],1.0)
    setValue(x[2],5.0)
    setValue(x[3],5.0)
    setValue(x[4],1.0)
    setValue(t, 100)

    status = solve(m)

    @test status == :Optimal
    @test_approx_eq_eps getValue(x[1]) 1.00000000 1e-5
    @test_approx_eq_eps getValue(x[2]) 4.74299963 1e-5
    @test_approx_eq_eps getValue(x[3]) 3.82114998 1e-5
    @test_approx_eq_eps getValue(x[4]) 1.37940829 1e-5
end

function test_nl_quadobj(nl_solver)
    # Solve a problem with quadratic objective with linear
    # constraints, but force it to use the nonlinear code.
    m = Model(solver=nl_solver)
    @defVar(m, 0.5 <= x <=  2)
    @defVar(m, 0.0 <= y <= 30)
    @setObjective(m, Min, (x+y)^2)
    @addNLConstraint(m, x + y >= 1)
    status = solve(m)
    
    @test status == :Optimal
    @test_approx_eq_eps m.objVal 1.0 1e-6
    @test_approx_eq_eps (getValue(x)+getValue(y)) 1.0 1e-6
end

function test_nl_quadcon(nl_solver)
    # Solve a problem with linear objective with quadratic
    # constraints, but force it to use the nonlinear code.
    m = Model(solver=nl_solver)
    @defVar(m, -2 <= x <= 2)
    @defVar(m, -2 <= y <= 2)
    @setNLObjective(m, Min, x - y)
    @addConstraint(m, x + x^2 + x*y + y^2 <= 1)
    status = solve(m)

    @test status == :Optimal
    @test_approx_eq_eps getObjectiveValue(m) -1-4/sqrt(3) 1e-6
    @test_approx_eq_eps (getValue(x) + getValue(y)) -1/3 1e-3
end

function test_nl_maxobj(nl_solver)
    # Solve a simple problem with a maximization objective
    m = Model(solver=nl_solver)
    @defVar(m, -2 <= x <= 2); setValue(x, -1.8)
    @defVar(m, -2 <= y <= 2); setValue(y,  1.5)
    @setNLObjective(m, Max, y - x)
    @addConstraint(m, x + x^2 + x*y + y^2 <= 1)
    status = solve(m)

    @test status == :Optimal
    @test_approx_eq_eps getObjectiveValue(m) 1+4/sqrt(3) 1e-6
    @test_approx_eq_eps (getValue(x) + getValue(y)) -1/3 1e-3
end

function test_nl_infeas(nl_solver)
    # (Attempt to) solve an infeasible problem
    m = Model(solver=nl_solver)
    n = 10
    @defVar(m, 0 <= x[i=1:n] <= 1)
    @setNLObjective(m, Max, x[n])
    for i in 1:n-1
        @addNLConstraint(m, x[i+1]-x[i] == 0.15)
    end
    status = solve(m)
    @test status == :Infeasible
end

function test_nl_unbnd(nl_solver)
    # (Attempt to) solve an unbounded problem
    m = Model(solver=nl_solver)
    @defVar(m, x >= 0)
    @setNLObjective(m, Max, x)
    @addNLConstraint(m, x >= 5)
    status = solve(m)
    @show status
    @test status == :Unbounded
end

#############################################################################
# Test that output is produced in correct MPB form
type DummyNLPSolver <: MathProgBase.AbstractMathProgSolver
end
type DummyNLPModel <: MathProgBase.AbstractMathProgModel
end
MathProgBase.model(s::DummyNLPSolver) = DummyNLPModel()
function MathProgBase.loadnonlinearproblem!(m::DummyNLPModel, numVar, numConstr, x_l, x_u, g_lb, g_ub, sense, d::MathProgBase.AbstractNLPEvaluator)
    MathProgBase.initialize(d, [:ExprGraph])
    objexpr = MathProgBase.obj_expr(d)
    @test (objexpr == :(x[1]^x[2])) || (objexpr == :(-1.0*x[1]+1.0*x[2]))
    @assert MathProgBase.isconstrlinear(d,1)
    @test MathProgBase.constr_expr(d,1) == :(2.0*x[1] + 1.0*x[2] <= 1.0)
    @test MathProgBase.constr_expr(d,2) == :(2.0*x[1]*x[1] + 1.0*x[2] + -2.0 >= 0)
    @test MathProgBase.constr_expr(d,3) == :(sin(x[1]) * cos(x[2]) - 5 == 0.0)
    @test MathProgBase.constr_expr(d,4) == :(1.0*x[1]*x[1] - 1.0 == 0.0)
    @test MathProgBase.constr_expr(d,5) == :(2.0*x[1]*x[1] - 2.0 == 0.0)
end
MathProgBase.setwarmstart!(m::DummyNLPModel,x) = nothing
MathProgBase.optimize!(m::DummyNLPModel) = nothing
MathProgBase.status(m::DummyNLPModel) = :Optimal
MathProgBase.getobjval(m::DummyNLPModel) = NaN
MathProgBase.getsolution(m::DummyNLPModel) = [1.0,1.0]
function test_nl_mpb()
    m = Model(solver=DummyNLPSolver())
    @defVar(m, x)
    @defVar(m, y)
    @setObjective(m, Min, -x+y)
    @addConstraint(m, 2x+y <= 1)
    @addConstraint(m, 2x^2+y >= 2)
    @addNLConstraint(m, sin(x)*cos(y) == 5)
    @addNLConstraint(m, nlconstr[i=1:2], i*x^2 == i)
    solve(m)

    @setNLObjective(m, Min, x^y)
    solve(m)
end

#############################################################################
function run_nl_tests(nl_solvers)
    if length(nl_solvers) == 0
        warn("  No nonlinear solvers available, skipping most tests!")
    end
    for (solver_name, nl_solver) in nl_solvers
        println("  Running with $solver_name")
        test_hs071(nl_solver)
        test_hs071_linobj(nl_solver)
        test_nl_quadobj(nl_solver)
        test_nl_quadcon(nl_solver)
        test_nl_maxobj(nl_solver)
        test_nl_infeas(nl_solver)
        test_nl_unbnd(nl_solver)
    end
    println("  Testing MPB interface")
    test_nl_mpb()
end

function load_nl_solvers()
    nl_solvers = Any[]
    if Pkg.installed("Ipopt") != nothing
        eval(Expr(:import,:Ipopt))
        push!(nl_solvers, ("Ipopt",Ipopt.IpoptSolver(print_level=0)))
    end
    if Pkg.installed("NLopt") != nothing
        eval(Expr(:import,:NLopt))
        push!(nl_solvers, ("NLopt",NLopt.NLoptSolver(algorithm=:LD_SLSQP)))
    end
    if Pkg.installed("KNITRO") != nothing
        eval(Expr(:import,:KNITRO))
        push!(nl_solvers, ("KNITRO",KNITRO.KnitroSolver()))
    end
    return nl_solvers
end

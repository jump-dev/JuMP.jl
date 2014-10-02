using JuMP
using Base.Test

# hs071
let
    # min x1 * x4 * (x1 + x2 + x3) + x3
    # st  x1 * x2 * x3 * x4 >= 25
    #     x1^2 + x2^2 + x3^2 + x4^2 = 40
    #     1 <= x1, x2, x3, x4 <= 5
    # Start at (1,5,5,1)
    # End at (1.000..., 4.743..., 3.821..., 1.379...)

    m = Model()

    @defVar(m, 1 <= x[1:4] <= 5)

    @setNLObjective(m, Min, x[1]*x[4]*(x[1]+x[2]+x[3]) + x[3])

    @addNLConstraint(m, x[1]*x[2]*x[3]*x[4] >= 25)
    @addNLConstraint(m, sum{x[i]^2,i=1:4} == 40)

    setValue(x[1],1.0)
    setValue(x[2],5.0)
    setValue(x[3],5.0)
    setValue(x[4],1.0)

    status = solve(m)

    @assert status == :Optimal
    @test_approx_eq_eps getValue(x[1]) 1.00000000 1e-5
    @test_approx_eq_eps getValue(x[2]) 4.74299963 1e-5
    @test_approx_eq_eps getValue(x[3]) 3.82114998 1e-5
    @test_approx_eq_eps getValue(x[4]) 1.37940829 1e-5
end

# hs071, linear objective
let
    # min x1 * x4 * (x1 + x2 + x3) + x3
    # st  x1 * x2 * x3 * x4 >= 25
    #     x1^2 + x2^2 + x3^2 + x4^2 = 40
    #     1 <= x1, x2, x3, x4 <= 5
    # Start at (1,5,5,1)
    # End at (1.000..., 4.743..., 3.821..., 1.379...)

    m = Model()

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

    @assert status == :Optimal
    @test_approx_eq_eps getValue(x[1]) 1.00000000 1e-5
    @test_approx_eq_eps getValue(x[2]) 4.74299963 1e-5
    @test_approx_eq_eps getValue(x[3]) 3.82114998 1e-5
    @test_approx_eq_eps getValue(x[4]) 1.37940829 1e-5
end

# quadratic objective
let
    m = Model()
    @defVar(m, 0.5 <= x <= 2)
    @defVar(m, 0 <= y <= 30)

    setObjective(m, :Min, (x+y)*(x+y))

    @addNLConstraint(m, x + y >= 1)

    status = solve(m)
    
    @assert status == :Optimal
    @test_approx_eq_eps m.objVal 1.0 1e-6
    @test_approx_eq_eps (getValue(x)+getValue(y)) 1.0 1e-6
end

# quadratic constraints
let
    m = Model()
    @defVar(m, -2 <= x <= 2)
    @defVar(m, -2 <= y <= 2)

    @setNLObjective(m, Min, x - y)
    addConstraint(m, x + x*x + x*y + y*y <= 1)

    status = solve(m)

    @assert status == :Optimal
    @test_approx_eq_eps m.objVal -1-4/sqrt(3) 1e-6
    @test_approx_eq_eps (getValue(x) + getValue(y)) -1/3 1e-3
end

# maximization sense
let
    m = Model()
    @defVar(m, -2 <= x <= 2)
    @defVar(m, -2 <= y <= 2)

    @setNLObjective(m, Max, y - x)
    addConstraint(m, x + x*x + x*y + y*y <= 1)

    status = solve(m)

    @assert status == :Optimal
    @test_approx_eq_eps m.objVal 1+4/sqrt(3) 1e-6
    @test_approx_eq_eps (getValue(x) + getValue(y)) -1/3 1e-3
end

# accessing expressions
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

end

MathProgBase.setwarmstart!(m::DummyNLPModel,x) = nothing
MathProgBase.optimize!(m::DummyNLPModel) = nothing
MathProgBase.status(m::DummyNLPModel) = :Optimal
MathProgBase.getobjval(m::DummyNLPModel) = NaN
MathProgBase.getsolution(m::DummyNLPModel) = [1.0,1.0]

let
    m = Model(solver=DummyNLPSolver())
    @defVar(m, x)
    @defVar(m, y)
    @setObjective(m, Min, -x+y)
    @addConstraint(m, 2x+y <= 1)
    @addConstraint(m, 2x^2+y >= 2)
    @addNLConstraint(m, sin(x)*cos(y) == 5)
    solve(m)

    @setNLObjective(m, Min, x^y)
    solve(m)
end

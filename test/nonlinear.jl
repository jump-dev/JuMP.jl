#  Copyright 2015, Iain Dunning, Joey Huchette, Miles Lubin, and contributors
#  This Source Code Form is subject to the terms of the Mozilla Public
#  License, v. 2.0. If a copy of the MPL was not distributed with this
#  file, You can obtain one at http://mozilla.org/MPL/2.0/.
#############################################################################
# JuMP
# An algebraic modelling langauge for Julia
# See http://github.com/JuliaOpt/JuMP.jl
#############################################################################
# test/nonlinear.jl
# Test general nonlinear
# Must be run as part of runtests.jl, as it needs a list of solvers.
#############################################################################
using JuMP, FactCheck

facts("[nonlinear] Test HS071 solves correctly") do
for nlp_solver in nlp_solvers
context("With solver $(typeof(nlp_solver))") do
    # hs071
    # Polynomial objective and constraints
    # min x1 * x4 * (x1 + x2 + x3) + x3
    # st  x1 * x2 * x3 * x4 >= 25
    #     x1^2 + x2^2 + x3^2 + x4^2 = 40
    #     1 <= x1, x2, x3, x4 <= 5
    # Start at (1,5,5,1)
    # End at (1.000..., 4.743..., 3.821..., 1.379...)
    m = Model(solver=nlp_solver)
    initval = [1,5,5,1]
    @defVar(m, 1 <= x[i=1:4] <= 5, start=initval[i])
    @setNLObjective(m, Min, x[1]*x[4]*(x[1]+x[2]+x[3]) + x[3])
    @addNLConstraint(m, x[1]*x[2]*x[3]*x[4] >= 25)
    @addNLConstraint(m, sum{x[i]^2,i=1:4} == 40)
    @fact MathProgBase.numconstr(m) --> 2
    status = solve(m)

    @fact status --> :Optimal
    @fact getValue(x)[:] --> roughly(
        [1.000000, 4.742999, 3.821150, 1.379408], 1e-5)
end; end; end


facts("[nonlinear] Test HS071 solves correctly, epigraph") do
for nlp_solver in nlp_solvers
context("With solver $(typeof(nlp_solver))") do
        # hs071, with epigraph formulation
        # Linear objective, nonlinear constraints
        # min t
        # st  t >= x1 * x4 * (x1 + x2 + x3) + x3
        #     ...
        m = Model(solver=nlp_solver)
        start = [1.0, 5.0, 5.0, 1.0]
        @defVar(m, 1 <= x[i=1:4] <= 5, start = start[i])
        @defVar(m, t, start = 100)
        @setObjective(m, Min, t)
        @addNLConstraint(m, t >= x[1]*x[4]*(x[1]+x[2]+x[3]) + x[3])
        @addNLConstraint(m, x[1]*x[2]*x[3]*x[4] >= 25)
        @addNLConstraint(m, sum{x[i]^2,i=1:4} == 40)
        status = solve(m)

        @fact status --> :Optimal
        @fact getValue(x)[:] --> roughly(
            [1.000000, 4.742999, 3.821150, 1.379408], 1e-5)
end; end; end

facts("[nonlinear] Test ifelse") do
for nlp_solver in nlp_solvers
(contains("$(typeof(nlp_solver))", "OsilSolver") || contains("$(typeof(nlp_solver))", "NLoptSolver")) && continue
context("With solver $(typeof(nlp_solver))") do
        m = Model(solver=nlp_solver)
        @defVar(m, x, start = 2)
        # minimizer at smooth point, solvers should be okay
        @setNLObjective(m, Min, ifelse( x <= 1, x^2, x) )
        status = solve(m)

        @fact status --> :Optimal
        @fact getValue(x) --> roughly(0.0, 1e-5)
end; end; end

facts("[nonlinear] Accepting fixed variables") do
for nlp_solver in convex_nlp_solvers
context("With solver $(typeof(nlp_solver))") do
    m = Model(solver=nlp_solver)
    @defVar(m, x == 0)
    @defVar(m, y ≥ 0)
    @setObjective(m, Min, y)
    @addNLConstraint(m, y ≥ x^2)
    for α in 1:4
        setValue(x, α)
        solve(m)
        @fact getValue(y) --> roughly(α^2, 1e-6)
    end
end; end; end

facts("[nonlinear] Test QP solve through NL pathway") do
for nlp_solver in convex_nlp_solvers
context("With solver $(typeof(nlp_solver))") do
    # Solve a problem with quadratic objective with linear
    # constraints, but force it to use the nonlinear code.
    m = Model(solver=nlp_solver)
    @defVar(m, 0.5 <= x <=  2)
    @defVar(m, 0.0 <= y <= 30)
    @setObjective(m, Min, (x+y)^2)
    param = [1.0]
    @addNLConstraint(m, x + y >= param[1])
    status = solve(m)

    @fact status --> :Optimal
    @fact m.objVal --> roughly(1.0, 1e-6)
    @fact getValue(x)+getValue(y) --> roughly(1.0, 1e-6)

    # sneaky problem modification
    param[1] = 10
    @fact m.internalModelLoaded --> true
    status = solve(m)
    @fact m.objVal --> roughly(10.0^2, 1e-6)
    @fact getValue(x)+getValue(y) --> roughly(10.0, 1e-6)

end; end; end


facts("[nonlinear] Test quad con solve through NL pathway") do
for nlp_solver in convex_nlp_solvers
context("With solver $(typeof(nlp_solver))") do
    # Solve a problem with linear objective with quadratic
    # constraints, but force it to use the nonlinear code.
    m = Model(solver=nlp_solver)
    @defVar(m, -2 <= x <= 2)
    @defVar(m, -2 <= y <= 2)
    @setNLObjective(m, Min, x - y)
    @addConstraint(m, x + x^2 + x*y + y^2 <= 1)
    status = solve(m)

    @fact status --> :Optimal
    @fact getObjectiveValue(m) --> roughly(-1-4/sqrt(3), 1e-6)
    @fact getValue(x) + getValue(y) --> roughly(-1/3, 1e-3)
end; end; end

facts("[nonlinear] Test two-sided nonlinear constraints") do
for nlp_solver in convex_nlp_solvers
context("With solver $(typeof(nlp_solver))") do
    m = Model(solver=nlp_solver)
    @defVar(m, x)
    @setNLObjective(m, Max, x)
    l = -1
    u = 1
    @addNLConstraint(m, l <= x <= u)
    status = solve(m)

    @fact status --> :Optimal
    @fact getObjectiveValue(m) --> roughly(u, 1e-6)

    @setNLObjective(m, Min, x)
    status = solve(m)

    @fact status --> :Optimal
    @fact getObjectiveValue(m) --> roughly(l, 1e-6)
end; end; end

facts("[nonlinear] Test mixed integer nonlinear problems") do
for minlp_solver in minlp_solvers
context("With solver $(typeof(minlp_solver))") do
    ## Solve test problem 1 (Synthesis of processing system) in
     # M. Duran & I.E. Grossmann, "An outer approximation algorithm for
     # a class of mixed integer nonlinear programs", Mathematical
     # Programming 36, pp. 307-339, 1986.  The problem also appears as
     # problem synthes1 in the MacMINLP test set.
    m = Model(solver=minlp_solver)
    x_U = [2,2,1]
    @defVar(m, x_U[i] >= x[i=1:3] >= 0)
    @defVar(m, y[4:6], Bin)
    @setNLObjective(m, Min, 10 + 10*x[1] - 7*x[3] + 5*y[4] + 6*y[5] + 8*y[6] - 18*log(x[2]+1) - 19.2*log(x[1]-x[2]+1))
    @addNLConstraints(m, begin
        0.8*log(x[2] + 1) + 0.96*log(x[1] - x[2] + 1) - 0.8*x[3] >= 0
        log(x[2] + 1) + 1.2*log(x[1] - x[2] + 1) - x[3] - 2*y[6] >= -2
        x[2] - x[1] <= 0
        x[2] - 2*y[4] <= 0
        x[1] - x[2] - 2*y[5] <= 0
        y[4] + y[5] <= 1
    end)
    status = solve(m)

    @fact status --> :Optimal
    @fact getObjectiveValue(m) --> roughly(6.00976, 1e-5)
    @fact getValue(x)[:] --> roughly([1.30098, 0.0, 1.0], 1e-5)
    @fact getValue(y)[:] --> roughly([0.0, 1.0, 0.0], 1e-5)
end; end; end

facts("[nonlinear] Test continuous relaxation of minlp test problem") do
for nlp_solver in nlp_solvers
context("With solver $(typeof(nlp_solver))") do
    ## Solve continuous relaxation of test problem 1 (Synthesis of processing system) in
     # M. Duran & I.E. Grossmann, "An outer approximation algorithm for
     # a class of mixed integer nonlinear programs", Mathematical
     # Programming 36, pp. 307-339, 1986.  The problem also appears as
     # problem synthes1 in the MacMINLP test set.
     # Introduce auxiliary nonnegative variable for the x[1]-x[2]+1 term
    m = Model(solver=nlp_solver)
    x_U = [2,2,1]
    @defVar(m, x_U[i] >= x[i=1:3] >= 0)
    @defVar(m, 1 >= y[4:6] >= 0)
    @defVar(m, z >= 0, start=1)
    @setNLObjective(m, Min, 10 + 10*x[1] - 7*x[3] + 5*y[4] + 6*y[5] + 8*y[6] - 18*log(x[2]+1) - 19.2*log(z))
    @addNLConstraints(m, begin
        0.8*log(x[2] + 1) + 0.96*log(z) - 0.8*x[3] >= 0
        log(x[2] + 1) + 1.2*log(z) - x[3] - 2*y[6] >= -2
        x[2] - x[1] <= 0
        x[2] - 2*y[4] <= 0
        x[1] - x[2] - 2*y[5] <= 0
        y[4] + y[5] <= 1
        x[1] - x[2] + 1 == z
    end)
    status = solve(m)

    @fact status --> :Optimal
    @fact getObjectiveValue(m) --> roughly(0.7593, 5e-5)
    @fact getValue(x)[:] --> roughly([1.1465, 0.54645, 1.0], 2e-4)
    @fact getValue(y)[:] --> roughly([0.2732, 0.3, 0.0], 2e-4)
    @fact getValue(z) --> roughly(1.6, 2e-4)
end; end; end

facts("[nonlinear] Test maximization objective") do
for nlp_solver in convex_nlp_solvers
context("With solver $(typeof(nlp_solver))") do
    # Solve a simple problem with a maximization objective
    m = Model(solver=nlp_solver)
    @defVar(m, -2 <= x <= 2); setValue(x, -1.8)
    @defVar(m, -2 <= y <= 2); setValue(y,  1.5)
    @setNLObjective(m, Max, y - x)
    @addConstraint(m, x + x^2 + x*y + y^2 <= 1)

    @fact solve(m) --> :Optimal
    @fact getObjectiveValue(m) --> roughly(1+4/sqrt(3), 1e-6)
    @fact getValue(x) + getValue(y) --> roughly(-1/3, 1e-3)
end; end; end

facts("[nonlinear] Test maximization objective (embedded expressions)") do
for nlp_solver in convex_nlp_solvers
context("With solver $(typeof(nlp_solver))") do
    m = Model(solver=nlp_solver)
    @defVar(m, -2 <= x <= 2); setValue(x, -1.8)
    @defVar(m, -2 <= y <= 2); setValue(y,  1.5)
    @setNLObjective(m, Max, y - x)
    @defNLExpr(quadexpr, x + x^2 + x*y + y^2)
    @addNLConstraint(m, quadexpr <= 1)

    @fact solve(m) --> :Optimal
    @fact getObjectiveValue(m) --> roughly(1+4/sqrt(3), 1e-6)
    @fact getValue(x) + getValue(y) --> roughly(-1/3, 1e-3)
    @fact getValue(quadexpr) --> roughly(1, 1e-5)
end; end; end


facts("[nonlinear] Test infeasibility detection") do
for nlp_solver in convex_nlp_solvers
contains(string(typeof(nlp_solver)),"NLoptSolver") && continue
context("With solver $(typeof(nlp_solver))") do
    # (Attempt to) solve an infeasible problem
    m = Model(solver=nlp_solver)
    n = 10
    @defVar(m, 0 <= x[i=1:n] <= 1)
    @setNLObjective(m, Max, x[n])
    for i in 1:n-1
        @addNLConstraint(m, x[i+1]-x[i] == 0.15)
    end
    @fact solve(m, suppress_warnings=true) --> :Infeasible
end; end; end


facts("[nonlinear] Test unboundedness detection") do
for nlp_solver in convex_nlp_solvers
contains(string(typeof(nlp_solver)),"NLoptSolver") && continue
context("With solver $(typeof(nlp_solver))") do
    # (Attempt to) solve an unbounded problem
    m = Model(solver=nlp_solver)
    @defVar(m, x >= 0)
    @setNLObjective(m, Max, x)
    @addNLConstraint(m, x >= 5)
    @fact solve(m, suppress_warnings=true) --> :Unbounded
end; end; end

facts("[nonlinear] Test entropy maximization") do
for nlp_solver in convex_nlp_solvers
context("With solver $(typeof(nlp_solver))") do
    m = Model(solver=nlp_solver)
    N = 3
    @defVar(m, x[1:N] >= 0, start = 1)
    @defNLExpr(entropy[i=1:N], -x[i]*log(x[i]))
    @setNLObjective(m, Max, sum{entropy[i], i = 1:N})
    @addConstraint(m, sum(x) == 1)

    @fact solve(m) --> :Optimal
    @fact norm(getValue(x)[:] - [1/3,1/3,1/3]) --> roughly(0.0, 1e-4)
end; end; end

facts("[nonlinear] Test entropy maximization (reformulation)") do
for nlp_solver in convex_nlp_solvers
context("With solver $(typeof(nlp_solver))") do
    m = Model(solver=nlp_solver)
    idx = [1,2,3,4]
    @defVar(m, x[idx] >= 0, start = 1)
    @defVar(m, z[1:4], start = 0)
    @defNLExpr(entropy[i=idx], -x[i]*log(x[i]))
    @setNLObjective(m, Max, sum{z[i], i = 1:2} + sum{z[i]/2, i=3:4})
    @addNLConstraint(m, z_constr1[i=1], z[i] <= entropy[i])
    @addNLConstraint(m, z_constr1[i=2], z[i] <= entropy[i]) # duplicate expressions
    @addNLConstraint(m, z_constr2[i=3:4], z[i] <= 2*entropy[i])
    @addConstraint(m, sum(x) == 1)

    @fact solve(m) --> :Optimal
    @fact norm([getValue(x[i]) for i in idx] - [1/4,1/4,1/4,1/4]) --> roughly(0.0, 1e-4)
    @fact getValue(entropy[1]) --> roughly(-(1/4)*log(1/4), 1e-4)
    @defNLExpr(zexpr[i=1:4], z[i])
    @fact getValue(zexpr[1]) --> roughly(-(1/4)*log(1/4), 1e-4)
end; end; end

facts("[nonlinear] Test nonlinear duals") do
for nlp_solver in nlp_solvers
applicable(MathProgBase.getconstrduals, MathProgBase.NonlinearModel(nlp_solver)) || continue
context("With solver $(typeof(nlp_solver))") do
    modA = Model(solver=nlp_solver)
    @defVar(modA, x >= 0)
    @defVar(modA, y <= 5)
    @defVar(modA, 2 <= z <= 4)
    @defVar(modA, 0 <= r[i=3:6] <= i)
    @setNLObjective(modA, Min, -((x + y)/2.0 + 3.0)/3.0 - z - r[3])
    @addConstraint(modA, cons1, x+y >= 2)
    @addConstraint(modA, cons2, sum{r[i],i=3:5} <= (2 - x)/2.0)
    @addNLConstraint(modA, cons3, 7.0*y <= z + r[6]/1.9)

    # Solution
    @fact solve(modA) --> :Optimal
    @fact getObjectiveValue(modA) --> roughly(-5.8446115, 1e-6)
    @fact getValue(x)       --> roughly(0.9774436, 1e-6)
    @fact getValue(y)       --> roughly(1.0225563, 1e-6)
    @fact getValue(z)       --> roughly(4.0, 1e-6)
    @fact getValue(r)[3]    --> roughly(0.5112781, 1e-6)
    @fact getValue(r)[4]    --> roughly(0.0, 1e-6)
    @fact getValue(r)[5]    --> roughly(0.0, 1e-6)
    @fact getValue(r)[6]    --> roughly(6.0, 1e-6)

    # Reduced costs
    @fact getDual(x)    --> roughly( 0.0, 1e-6)
    @fact getDual(y)    --> roughly( 0.0, 1e-6)
    @fact getDual(z)    --> roughly(-1.0714286, 1e-6)
    @fact getDual(r)[3] --> roughly( 0.0, 1e-6)
    @fact getDual(r)[4] --> roughly(1.0, 1e-6)
    @fact getDual(r)[5] --> roughly(1.0, 1e-6)
    @fact getDual(r)[6] --> roughly(-0.03759398, 1e-6)

    # Row duals
    @fact getDual(cons1) --> roughly( 0.333333, 1e-6)
    @fact getDual(cons2) --> roughly(-1.0, 1e-6)
    @fact getDual(cons3) --> roughly(-0.0714286, 1e-6)
end; end; end

facts("[nonlinear] Test nonlinear duals (Max)") do
for nlp_solver in nlp_solvers
applicable(MathProgBase.getconstrduals, MathProgBase.NonlinearModel(nlp_solver)) || continue
context("With solver $(typeof(nlp_solver))") do
    modA = Model(solver=nlp_solver)
    @defVar(modA, x >= 0)
    @defVar(modA, y <= 5)
    @defVar(modA, 2 <= z <= 4)
    @defVar(modA, 0 <= r[i=3:6] <= i)
    @setNLObjective(modA, Max, ((x + y)/2.0 + 3.0)/3.0 + z + r[3])
    @addConstraint(modA, cons1, x+y >= 2)
    @addConstraint(modA, cons2, sum{r[i],i=3:5} <= (2 - x)/2.0)
    @addNLConstraint(modA, cons3, 7.0*y <= z + r[6]/1.9)

    # Solution
    @fact solve(modA) --> :Optimal
    @fact getObjectiveValue(modA) --> roughly(5.8446115, 1e-6)
    @fact getValue(x)       --> roughly(0.9774436, 1e-6)
    @fact getValue(y)       --> roughly(1.0225563, 1e-6)
    @fact getValue(z)       --> roughly(4.0, 1e-6)
    @fact getValue(r)[3]    --> roughly(0.5112781, 1e-6)
    @fact getValue(r)[4]    --> roughly(0.0, 1e-6)
    @fact getValue(r)[5]    --> roughly(0.0, 1e-6)
    @fact getValue(r)[6]    --> roughly(6.0, 1e-6)

    # Reduced costs
    @fact getDual(x)    --> roughly( 0.0, 1e-6)
    @fact getDual(y)    --> roughly( 0.0, 1e-6)
    @fact getDual(z)    --> roughly(1.0714286, 1e-6)
    @fact getDual(r)[3] --> roughly( 0.0, 1e-6)
    @fact getDual(r)[4] --> roughly(-1.0, 1e-6)
    @fact getDual(r)[5] --> roughly(-1.0, 1e-6)
    @fact getDual(r)[6] --> roughly(0.03759398, 1e-6)

    # Row duals
    @fact getDual(cons1) --> roughly(-0.333333, 1e-6)
    @fact getDual(cons2) --> roughly(1.0, 1e-6)
    @fact getDual(cons3) --> roughly(0.0714286, 1e-6)
end; end; end


#############################################################################
# Test that output is produced in correct MPB form
type DummyNLPSolver <: MathProgBase.AbstractMathProgSolver
end
type DummyNLPModel <: MathProgBase.AbstractNonlinearModel
end
MathProgBase.NonlinearModel(s::DummyNLPSolver) = DummyNLPModel()
function MathProgBase.loadproblem!(m::DummyNLPModel, numVar, numConstr, x_l, x_u, g_lb, g_ub, sense, d::MathProgBase.AbstractNLPEvaluator)
    MathProgBase.initialize(d, [:ExprGraph])
    objexpr = MathProgBase.obj_expr(d)
    facts("[nonlinear] Test NL MPB interface ($objexpr)") do
        @fact objexpr --> anyof(:(x[1]^x[2]), :(-1.0*x[1]+1.0*x[2]))
        @fact MathProgBase.isconstrlinear(d,1) --> true
        @fact MathProgBase.isconstrlinear(d,3) --> true
        @fact MathProgBase.constr_expr(d,1) --> :(2.0*x[1] + 1.0*x[2] <= 1.0)
        @fact MathProgBase.constr_expr(d,2) --> :(2.0*x[1] + 1.0*x[2] <= 0.0)
        @fact MathProgBase.constr_expr(d,3) --> :(-5.0 <= 2.0*x[1] + 1.0*x[2] <= 5.0)
        if numConstr > 3
            @fact MathProgBase.constr_expr(d,4) --> :(2.0*x[1]*x[1] + 1.0*x[2] + -2.0 >= 0)
            @fact MathProgBase.constr_expr(d,5) --> :(sin(x[1]) * cos(x[2]) - 5 == 0.0)
            @fact MathProgBase.constr_expr(d,6) --> :(1.0*x[1]^2 - 1.0 == 0.0)
            @fact MathProgBase.constr_expr(d,7) --> :(2.0*x[1]^2 - 2.0 == 0.0)
            @fact MathProgBase.constr_expr(d,8) --> :(-0.5 <= sin(x[1]) <= 0.5)
        end
    end
end
MathProgBase.setwarmstart!(m::DummyNLPModel,x) = nothing
MathProgBase.optimize!(m::DummyNLPModel) = nothing
MathProgBase.status(m::DummyNLPModel) = :Optimal
MathProgBase.getobjval(m::DummyNLPModel) = NaN
MathProgBase.getsolution(m::DummyNLPModel) = [1.0,1.0]
MathProgBase.setvartype!(m::DummyNLPModel,vartype) = @fact any(vartype .== :Fixed) --> false
function test_nl_mpb()
    m = Model(solver=DummyNLPSolver())
    @defVar(m, x == 1)
    @defVar(m, y, Bin)
    @setObjective(m, Min, -x+y)
    @addConstraint(m, 2x+y <= 1)
    @addConstraint(m, 2x+y <= 0)
    @addConstraint(m, -5 <= 2x+y <= 5)
    #solve(m) # FIXME maybe?
    lb,ub = getConstraintBounds(m)
    @fact lb --> [-Inf,-Inf,-5.0]
    @fact ub --> [1.0,-0.0,5.0]

    @addConstraint(m, 2x^2+y >= 2)
    @addNLConstraint(m, sin(x)*cos(y) == 5)
    @addNLConstraint(m, nlconstr[i=1:2], i*x^2 == i)
    @addNLConstraint(m, -0.5 <= sin(x) <= 0.5)
    solve(m)

    @setNLObjective(m, Min, x^y)
    solve(m)

    lb,ub = getConstraintBounds(m)
    @fact lb --> [-Inf,-Inf,-5.0,0.0,0.0,0.0,0.0,-0.5]
    @fact ub --> [1.0,-0.0,5.0,Inf,0.0,0.0,0.0,0.5]
end
test_nl_mpb()

facts("[nonlinear] Expression graph for linear problem") do
    m = Model()
    @defVar(m, x)
    @addConstraint(m, 0 <= x <= 1)
    @setObjective(m, Max, x)
    d = JuMPNLPEvaluator(m)
    MathProgBase.initialize(d, [:ExprGraph])
    @fact MathProgBase.obj_expr(d) --> :(+(1.0 * x[1]))
end

facts("[nonlinear] Hessians through MPB") do
    # Issue 435
    m = Model()
    @defVar(m, a, start = 1)
    @defVar(m, b, start = 2)
    @defVar(m, c, start = 3)

    @defNLExpr(foo, a * b + c^2)

    @setNLObjective(m, Min, foo)
    d = JuMPNLPEvaluator(m)
    MathProgBase.initialize(d, [:Hess])
    I,J = MathProgBase.hesslag_structure(d)
    V = zeros(length(I))
    MathProgBase.eval_hesslag(d, V, m.colVal, 1.0, Float64[])
    hess_raw = sparse(I,J,V)
    # Convert from lower triangular
    hess_sparse = hess_raw + hess_raw' - sparse(diagm(diag(hess_raw)))
    @fact hess_sparse --> roughly([0.0 1.0 0.0; 1.0 0.0 0.0; 0.0 0.0 2.0])
end

facts("[nonlinear] Hess-vec through MPB") do
    m = Model()
    @defVar(m, a, start = 1)
    @defVar(m, b, start = 2)
    @defVar(m, c, start = 3)

    @setNLObjective(m, Min, a*b + c^2)
    @addConstraint(m, c*b <= 1)
    @addNLConstraint(m, a^2/2 <= 1)
    d = JuMPNLPEvaluator(m)
    MathProgBase.initialize(d, [:HessVec])
    h = zeros(3)
    v = [2.4,3.5,1.2]
    MathProgBase.eval_hesslag_prod(d, h, m.colVal, v, 1.0, [2.0,3.0])
    correct = [3.0 1.0 0.0; 1.0 0.0 2.0; 0.0 2.0 2.0]*v
    @fact h --> roughly(correct)
end

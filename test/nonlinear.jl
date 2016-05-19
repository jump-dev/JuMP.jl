#  Copyright 2016, Iain Dunning, Joey Huchette, Miles Lubin, and contributors
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
#############################################################################
using JuMP, FactCheck

# If solvers not loaded, load them (i.e running just these tests)
!isdefined(:nlp_solvers) && include("solvers.jl")

facts("[nonlinear] Test getvalue on arrays") do
    m = Model()
    @variable(m, x, start = π/2)
    @NLexpression(m, f1, sin(x))
    @NLexpression(m, f2, sin(2x))
    @NLexpression(m, f3[i=1:2], sin(i*x))

    @fact getvalue(f1) --> roughly(1, 1e-5)
    @fact getvalue(f2) --> roughly(0, 1e-5)

    @fact getvalue([f1, f2]) --> getvalue(f3)

    v = [1.0, 2.0]
    @NLparameter(m, vparam[i=1:2] == v[i])
    @fact getvalue(vparam) --> v
    v[1] = 3.0
    setvalue(vparam, v)
    @fact getvalue(vparam) --> v
end

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
    @variable(m, 1 <= x[i=1:4] <= 5, start=initval[i])
    @NLobjective(m, Min, x[1]*x[4]*(x[1]+x[2]+x[3]) + x[3])
    @NLconstraint(m, x[1]*x[2]*x[3]*x[4] >= 25)
    @NLconstraint(m, sum{x[i]^2,i=1:4} == 40)
    @fact MathProgBase.numconstr(m) --> 2
    status = solve(m)

    @fact status --> :Optimal
    @fact getvalue(x)[:] --> roughly(
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
        @variable(m, 1 <= x[i=1:4] <= 5, start = start[i])
        @variable(m, t, start = 100)
        @objective(m, Min, t)
        @NLconstraint(m, t >= x[1]*x[4]*(x[1]+x[2]+x[3]) + x[3])
        @NLconstraint(m, x[1]*x[2]*x[3]*x[4] >= 25)
        @NLconstraint(m, sum{x[i]^2,i=1:4} == 40)
        status = solve(m)

        @fact status --> :Optimal
        @fact getvalue(x)[:] --> roughly(
            [1.000000, 4.742999, 3.821150, 1.379408], 1e-5)
end; end; end

facts("[nonlinear] Test ifelse") do
for nlp_solver in nlp_solvers
(contains("$(typeof(nlp_solver))", "OsilSolver") || contains("$(typeof(nlp_solver))", "NLoptSolver")) && continue
context("With solver $(typeof(nlp_solver))") do
        m = Model(solver=nlp_solver)
        @variable(m, x, start = 2)
        # minimizer at smooth point, solvers should be okay
        @NLobjective(m, Min, ifelse( x <= 1, x^2, x) )
        status = solve(m)

        @fact status --> :Optimal
        @fact getvalue(x) --> roughly(0.0, 1e-5)
end; end; end

facts("[nonlinear] Accepting fixed variables") do
for nlp_solver in convex_nlp_solvers
for simplify in [true, false]
context("With solver $(typeof(nlp_solver)), simplify = $simplify") do
    m = Model(solver=nlp_solver, simplify_nonlinear_expressions=simplify)
    @variable(m, x == 0)
    @variable(m, y ≥ 0)
    @objective(m, Min, y)
    @NLconstraint(m, y ≥ x^2)
    EnableNLPResolve()
    for α in 1:4
        setvalue(x, α)
        solve(m)
        @fact getvalue(y) --> roughly(α^2, 1e-6)
    end
    DisableNLPResolve()
end; end; end; end

facts("[nonlinear] Test QP solve through NL pathway") do
for nlp_solver in convex_nlp_solvers
context("With solver $(typeof(nlp_solver))") do
    # Solve a problem with quadratic objective with linear
    # constraints, but force it to use the nonlinear code.
    m = Model(solver=nlp_solver)
    @variable(m, 0.5 <= x <=  2)
    @variable(m, 0.0 <= y <= 30)
    @NLparameter(m, param == 1.0)
    @objective(m, Min, (x+y)^2)
    @NLconstraint(m, x + y >= param)
    status = solve(m)

    @fact status --> :Optimal
    @fact m.objVal --> roughly(1.0, 1e-6)
    @NLexpression(m, lhs, x+y)
    @fact getvalue(x)+getvalue(y) --> roughly(1.0, 1e-6)
    @fact getvalue(lhs) --> roughly(1.0, 1e-6)

    setvalue(param,10)
    @fact m.internalModelLoaded --> true
    status = solve(m)
    @fact m.objVal --> roughly(10.0^2, 1e-6)
    @fact getvalue(x)+getvalue(y) --> roughly(10.0, 1e-6)
    @fact getvalue(lhs) --> roughly(10.0, 1e-6)

end; end; end


facts("[nonlinear] Test quad con solve through NL pathway") do
for nlp_solver in convex_nlp_solvers
context("With solver $(typeof(nlp_solver))") do
    # Solve a problem with linear objective with quadratic
    # constraints, but force it to use the nonlinear code.
    m = Model(solver=nlp_solver)
    @variable(m, -2 <= x <= 2)
    @variable(m, -2 <= y <= 2)
    @NLobjective(m, Min, x - y)
    @constraint(m, x + x^2 + x*y + y^2 <= 1)
    status = solve(m)

    @fact status --> :Optimal
    @fact getobjectivevalue(m) --> roughly(-1-4/sqrt(3), 1e-6)
    @fact getvalue(x) + getvalue(y) --> roughly(-1/3, 1e-3)
end; end; end

facts("[nonlinear] Test resolve with parameter") do
for nlp_solver in convex_nlp_solvers
for simplify in [true,false]
context("With solver $(typeof(nlp_solver)), simplify = $simplify") do
    m = Model(solver=nlp_solver, simplify_nonlinear_expressions=simplify)
    @variable(m, z)
    @NLparameter(m, x == 1.0)
    @NLobjective(m, Min, (z-x)^2)
    status = solve(m)
    @fact status --> :Optimal
    @fact getvalue(z) --> roughly(1.0, 1e-3)

    setvalue(x, 5.0)
    status = solve(m)
    @fact status --> :Optimal
    @fact getvalue(z) --> roughly(5.0, 1e-3)
end; end; end; end

facts("[nonlinear] Test two-sided nonlinear constraints") do
for nlp_solver in convex_nlp_solvers
context("With solver $(typeof(nlp_solver))") do
    m = Model(solver=nlp_solver)
    @variable(m, x)
    @NLobjective(m, Max, x)
    l = -1
    u = 1
    @NLconstraint(m, l <= x <= u)
    status = solve(m)

    @fact status --> :Optimal
    @fact getobjectivevalue(m) --> roughly(u, 1e-6)

    @NLobjective(m, Min, x)
    status = solve(m)

    @fact status --> :Optimal
    @fact getobjectivevalue(m) --> roughly(l, 1e-6)
end; end; end

facts("[nonlinear] Quadratic equality constraints") do
for nlp_solver in nlp_solvers
context("With solver $(typeof(nlp_solver))") do
    m = Model(solver=nlp_solver)
    @variable(m, 0 <= x[1:2] <= 1)
    @constraint(m, x[1]^2 + x[2]^2 == 1/2)
    @NLobjective(m, Max, x[1] - x[2])
    status = solve(m)

    @fact status --> :Optimal
    @fact getvalue(x) --> roughly([sqrt(1/2), 0], 1e-6)
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
    @variable(m, x_U[i] >= x[i=1:3] >= 0)
    @variable(m, y[4:6], Bin)
    @NLobjective(m, Min, 10 + 10*x[1] - 7*x[3] + 5*y[4] + 6*y[5] + 8*y[6] - 18*log(x[2]+1) - 19.2*log(x[1]-x[2]+1))
    @NLconstraints(m, begin
        0.8*log(x[2] + 1) + 0.96*log(x[1] - x[2] + 1) - 0.8*x[3] >= 0
        log(x[2] + 1) + 1.2*log(x[1] - x[2] + 1) - x[3] - 2*y[6] >= -2
        x[2] - x[1] <= 0
        x[2] - 2*y[4] <= 0
        x[1] - x[2] - 2*y[5] <= 0
        y[4] + y[5] <= 1
    end)
    status = solve(m)

    @fact status --> :Optimal
    @fact getobjectivevalue(m) --> roughly(6.00976, 1e-5)
    @fact getvalue(x)[:] --> roughly([1.30098, 0.0, 1.0], 1e-5)
    @fact getvalue(y)[:] --> roughly([0.0, 1.0, 0.0], 1e-5)
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
    @variable(m, x_U[i] >= x[i=1:3] >= 0)
    @variable(m, 1 >= y[4:6] >= 0)
    @variable(m, z >= 0, start=1)
    @NLobjective(m, Min, 10 + 10*x[1] - 7*x[3] + 5*y[4] + 6*y[5] + 8*y[6] - 18*log(x[2]+1) - 19.2*log(z))
    @NLconstraints(m, begin
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
    @fact getobjectivevalue(m) --> roughly(0.7593, 5e-5)
    @fact getvalue(x)[:] --> roughly([1.1465, 0.54645, 1.0], 2e-4)
    @fact getvalue(y)[:] --> roughly([0.2732, 0.3, 0.0], 2e-4)
    @fact getvalue(z) --> roughly(1.6, 2e-4)
end; end; end

facts("[nonlinear] Test maximization objective") do
for nlp_solver in convex_nlp_solvers
context("With solver $(typeof(nlp_solver))") do
    # Solve a simple problem with a maximization objective
    m = Model(solver=nlp_solver)
    @variable(m, -2 <= x <= 2); setvalue(x, -1.8)
    @variable(m, -2 <= y <= 2); setvalue(y,  1.5)
    @NLobjective(m, Max, y - x)
    @constraint(m, x + x^2 + x*y + y^2 <= 1)

    @fact solve(m) --> :Optimal
    @fact getobjectivevalue(m) --> roughly(1+4/sqrt(3), 1e-6)
    @fact getvalue(x) + getvalue(y) --> roughly(-1/3, 1e-3)
end; end; end

facts("[nonlinear] Test maximization objective (embedded expressions)") do
for nlp_solver in convex_nlp_solvers
for simplify in [true,false]
context("With solver $(typeof(nlp_solver)), simplify = $simplify") do
    m = Model(solver=nlp_solver, simplify_nonlinear_expressions=simplify)
    @variable(m, -2 <= x <= 2); setvalue(x, -1.8)
    @variable(m, -2 <= y <= 2); setvalue(y,  1.5)
    @NLobjective(m, Max, y - x)
    @NLexpression(m, quadexpr, x + x^2 + x*y + y^2)
    @NLconstraint(m, quadexpr <= 1)

    @fact solve(m) --> :Optimal
    @fact getobjectivevalue(m) --> roughly(1+4/sqrt(3), 1e-6)
    @fact getvalue(x) + getvalue(y) --> roughly(-1/3, 1e-3)
    @fact getvalue(quadexpr) --> roughly(1, 1e-5)
    @NLexpression(m, quadexpr2, x + x^2 + x*y + y^2)
    @fact getvalue(quadexpr2) --> roughly(1, 1e-5)
end; end; end; end


facts("[nonlinear] Test infeasibility detection") do
for nlp_solver in convex_nlp_solvers
contains(string(typeof(nlp_solver)),"NLoptSolver") && continue
context("With solver $(typeof(nlp_solver))") do
    # (Attempt to) solve an infeasible problem
    m = Model(solver=nlp_solver)
    n = 10
    @variable(m, 0 <= x[i=1:n] <= 1)
    @NLobjective(m, Max, x[n])
    for i in 1:n-1
        @NLconstraint(m, x[i+1]-x[i] == 0.15)
    end
    @fact solve(m, suppress_warnings=true) --> :Infeasible
end; end; end


facts("[nonlinear] Test unboundedness detection") do
for nlp_solver in convex_nlp_solvers
contains(string(typeof(nlp_solver)),"NLoptSolver") && continue
context("With solver $(typeof(nlp_solver))") do
    # (Attempt to) solve an unbounded problem
    m = Model(solver=nlp_solver)
    @variable(m, x >= 0)
    @NLobjective(m, Max, x)
    @NLconstraint(m, x >= 5)
    @fact solve(m, suppress_warnings=true) --> :Unbounded
end; end; end

facts("[nonlinear] Test entropy maximization") do
for nlp_solver in convex_nlp_solvers
context("With solver $(typeof(nlp_solver))") do
    m = Model(solver=nlp_solver)
    N = 3
    @variable(m, x[1:N] >= 0, start = 1)
    @NLexpression(m, entropy[i=1:N], -x[i]*log(x[i]))
    @NLobjective(m, Max, sum{entropy[i], i = 1:N})
    @constraint(m, sum(x) == 1)

    @fact solve(m) --> :Optimal
    @fact norm(getvalue(x)[:] - [1/3,1/3,1/3]) --> roughly(0.0, 1e-4)
end; end; end

facts("[nonlinear] Test entropy maximization (reformulation)") do
for nlp_solver in convex_nlp_solvers
context("With solver $(typeof(nlp_solver))") do
    m = Model(solver=nlp_solver)
    idx = [1,2,3,4]
    @variable(m, x[idx] >= 0, start = 1)
    @variable(m, z[1:4], start = 0)
    @NLexpression(m, entropy[i=idx], -x[i]*log(x[i]))
    @NLobjective(m, Max, sum{z[i], i = 1:2} + sum{z[i]/2, i=3:4})
    @NLconstraint(m, z_constr1[i=1], z[i] <= entropy[i])
    @NLconstraint(m, z_constr1[i=2], z[i] <= entropy[i]) # duplicate expressions
    @NLconstraint(m, z_constr2[i=3:4], z[i] <= 2*entropy[i])
    @constraint(m, sum(x) == 1)

    @fact solve(m) --> :Optimal
    @fact norm([getvalue(x[i]) for i in idx] - [1/4,1/4,1/4,1/4]) --> roughly(0.0, 1e-4)
    @fact getvalue(entropy[1]) --> roughly(-(1/4)*log(1/4), 1e-4)
    @NLexpression(m, zexpr[i=1:4], z[i])
    @fact getvalue(zexpr[1]) --> roughly(-(1/4)*log(1/4), 1e-4)
end; end; end

facts("[nonlinear] Test derivatives of x^4, x < 0") do
for nlp_solver in convex_nlp_solvers
context("With solver $(typeof(nlp_solver))") do
    m = Model(solver=nlp_solver)
    @variable(m, x >= -1, start = -0.5)
    @NLobjective(m, Min, x^4)
    status = solve(m)

    @fact status --> :Optimal
    @fact getvalue(x) --> roughly(0.0, 1e-2)
end; end; end

facts("[nonlinear] Test nonlinear duals") do
for nlp_solver in nlp_solvers
for simplify in [true,false]
applicable(MathProgBase.getconstrduals, MathProgBase.NonlinearModel(nlp_solver)) || continue
context("With solver $(typeof(nlp_solver)), simplify = $simplify") do
    modA = Model(solver=nlp_solver, simplify_nonlinear_expressions=simplify)
    @variable(modA, x >= 0)
    @variable(modA, y <= 5)
    @variable(modA, 2 <= z <= 4)
    @variable(modA, 0 <= r[i=3:6] <= i)
    @NLobjective(modA, Min, -((x + y)/2.0 + 3.0)/3.0 - z - r[3])
    @constraint(modA, cons1, x+y >= 2)
    @constraint(modA, cons2, sum{r[i],i=3:5} <= (2 - x)/2.0)
    @NLconstraint(modA, cons3, 7.0*y <= z + r[6]/1.9)

    # Solution
    @fact solve(modA) --> :Optimal
    @fact getobjectivevalue(modA) --> roughly(-5.8446115, 1e-6)
    @fact getvalue(x)       --> roughly(0.9774436, 1e-6)
    @fact getvalue(y)       --> roughly(1.0225563, 1e-6)
    @fact getvalue(z)       --> roughly(4.0, 1e-6)
    @fact getvalue(r)[3]    --> roughly(0.5112781, 1e-6)
    @fact getvalue(r)[4]    --> roughly(0.0, 1e-6)
    @fact getvalue(r)[5]    --> roughly(0.0, 1e-6)
    @fact getvalue(r)[6]    --> roughly(6.0, 1e-6)

    # Reduced costs
    @fact getdual(x)    --> roughly( 0.0, 1e-6)
    @fact getdual(y)    --> roughly( 0.0, 1e-6)
    @fact getdual(z)    --> roughly(-1.0714286, 1e-6)
    @fact getdual(r)[3] --> roughly( 0.0, 1e-6)
    @fact getdual(r)[4] --> roughly(1.0, 1e-6)
    @fact getdual(r)[5] --> roughly(1.0, 1e-6)
    @fact getdual(r)[6] --> roughly(-0.03759398, 1e-6)

    # Row duals
    @fact getdual(cons1) --> roughly( 0.333333, 1e-6)
    @fact getdual(cons2) --> roughly(-1.0, 1e-6)
    @fact getdual(cons3) --> roughly(-0.0714286, 1e-6)
end; end; end; end

facts("[nonlinear] Test nonlinear duals (Max)") do
for nlp_solver in nlp_solvers
applicable(MathProgBase.getconstrduals, MathProgBase.NonlinearModel(nlp_solver)) || continue
context("With solver $(typeof(nlp_solver))") do
    modA = Model(solver=nlp_solver)
    @variable(modA, x >= 0)
    @variable(modA, y <= 5)
    @variable(modA, 2 <= z <= 4)
    @variable(modA, 0 <= r[i=3:6] <= i)
    @NLobjective(modA, Max, ((x + y)/2.0 + 3.0)/3.0 + z + r[3])
    @constraint(modA, cons1, x+y >= 2)
    @constraint(modA, cons2, sum{r[i],i=3:5} <= (2 - x)/2.0)
    cons3 = @NLconstraint(modA, 7.0*y <= z + r[6]/1.9)

    # Solution
    @fact solve(modA) --> :Optimal
    @fact getobjectivevalue(modA) --> roughly(5.8446115, 1e-6)
    @fact getvalue(x)       --> roughly(0.9774436, 1e-6)
    @fact getvalue(y)       --> roughly(1.0225563, 1e-6)
    @fact getvalue(z)       --> roughly(4.0, 1e-6)
    @fact getvalue(r)[3]    --> roughly(0.5112781, 1e-6)
    @fact getvalue(r)[4]    --> roughly(0.0, 1e-6)
    @fact getvalue(r)[5]    --> roughly(0.0, 1e-6)
    @fact getvalue(r)[6]    --> roughly(6.0, 1e-6)

    # Reduced costs
    @fact getdual(x)    --> roughly( 0.0, 1e-6)
    @fact getdual(y)    --> roughly( 0.0, 1e-6)
    @fact getdual(z)    --> roughly(1.0714286, 1e-6)
    @fact getdual(r)[3] --> roughly( 0.0, 1e-6)
    @fact getdual(r)[4] --> roughly(-1.0, 1e-6)
    @fact getdual(r)[5] --> roughly(-1.0, 1e-6)
    @fact getdual(r)[6] --> roughly(0.03759398, 1e-6)

    # Row duals
    @fact getdual(cons1) --> roughly(-0.333333, 1e-6)
    @fact getdual(cons2) --> roughly(1.0, 1e-6)
    @fact getdual(cons3) --> roughly(0.0714286, 1e-6)
end; end; end

facts("[nonlinear] Test changing objectives") do
for nlp_solver in nlp_solvers
context("With solver $(typeof(nlp_solver))") do
    m = Model(solver=nlp_solver)
    @variable(m, x >= 0)
    @variable(m, y >= 0)
    @objective(m, Max, x+y)
    @NLconstraint(m, x+2y <= 1)
    @fact solve(m) --> :Optimal
    @fact getvalue(x) --> roughly(1.0,1e-4)
    @fact getvalue(y) --> roughly(0.0,1e-4)
    @fact getobjectivevalue(m) --> roughly(1.0,1e-4)

    @objective(m, Max, 2x+y)
    @fact solve(m) --> :Optimal
    @fact getvalue(x) --> roughly(1.0,1e-4)
    @fact getvalue(y) --> roughly(0.0,1e-4)
    @fact getobjectivevalue(m) --> roughly(2.0,1e-4)
end; end; end

facts("[nonlinear] Test Hessian chunking code") do
for nlp_solver in nlp_solvers
context("With solver $(typeof(nlp_solver))") do
    m = Model(solver=nlp_solver)
    @variable(m, x[1:18] >= 1, start = 1.2)
    @NLobjective(m, Min, prod{x[i],i=1:18})
    @fact solve(m) --> :Optimal
    @fact getvalue(x) --> roughly(ones(18),1e-4)
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
        if VERSION >= v"0.5.0-dev+3231"
            @fact MathProgBase.constr_expr(d,2) --> :(2.0*x[1] + 1.0*x[2] <= -0.0)
        else
            @fact MathProgBase.constr_expr(d,2) --> :(2.0*x[1] + 1.0*x[2] <= 0.0)
        end
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
    @variable(m, x == 1)
    @variable(m, y, Bin)
    @objective(m, Min, -x+y)
    @constraint(m, 2x+y <= 1)
    @constraint(m, 2x+y <= 0)
    @constraint(m, -5 <= 2x+y <= 5)
    #solve(m) # FIXME maybe?
    lb,ub = JuMP.constraintbounds(m)
    @fact lb --> [-Inf,-Inf,-5.0]
    @fact ub --> [1.0,-0.0,5.0]

    @constraint(m, 2x^2+y >= 2)
    @NLconstraint(m, sin(x)*cos(y) == 5)
    @NLconstraint(m, nlconstr[i=1:2], i*x^2 == i)
    @NLconstraint(m, -0.5 <= sin(x) <= 0.5)
    solve(m)

    @NLobjective(m, Min, x^y)
    solve(m)

    lb,ub = JuMP.constraintbounds(m)
    @fact lb --> [-Inf,-Inf,-5.0,0.0,0.0,0.0,0.0,-0.5]
    @fact ub --> [1.0,-0.0,5.0,Inf,0.0,0.0,0.0,0.5]
end
test_nl_mpb()

facts("[nonlinear] Expression graph for linear problem") do
    m = Model()
    @variable(m, x)
    @constraint(m, 0 <= x <= 1)
    @objective(m, Max, x)
    d = JuMP.NLPEvaluator(m)
    MathProgBase.initialize(d, [:ExprGraph])
    @fact MathProgBase.obj_expr(d) --> :(+(1.0 * x[1]))
end

facts("[nonlinear] Expression graph for ifelse") do
    m = Model()
    @variable(m, x, start = 2)
    @NLobjective(m, Min, ifelse( x <= 1, x^2, x) )
    d = JuMP.NLPEvaluator(m)
    MathProgBase.initialize(d, [:ExprGraph])
    @fact MathProgBase.obj_expr(d) --> :(ifelse( x[1] <= 1, x[1]^2, x[1]))
end

facts("[nonlinear] Hessians through MPB") do
    # Issue 435
    m = Model()
    @variable(m, a, start = 1)
    @variable(m, b, start = 2)
    @variable(m, c, start = 3)

    @NLexpression(m, foo, a * b + c^2)

    @NLobjective(m, Min, foo)
    d = JuMP.NLPEvaluator(m)
    MathProgBase.initialize(d, [:Hess])
    I,J = MathProgBase.hesslag_structure(d)
    V = zeros(length(I))
    MathProgBase.eval_hesslag(d, V, m.colVal, 1.0, Float64[])
    hess_raw = sparse(I,J,V)
    # Convert from lower triangular
    hess_sparse = hess_raw + hess_raw' - sparse(diagm(diag(hess_raw)))
    @fact hess_sparse --> roughly([0.0 1.0 0.0; 1.0 0.0 0.0; 0.0 0.0 2.0])

    # make sure we don't get NaNs in this case
    @NLobjective(m, Min, a * b + 3*c^2)
    d = JuMP.NLPEvaluator(m)
    MathProgBase.initialize(d, [:Hess])
    setvalue(c, -1.0)
    V = zeros(length(I))
    MathProgBase.eval_hesslag(d, V, m.colVal, 1.0, Float64[])
    hess_raw = sparse(I,J,V)
    hess_sparse = hess_raw + hess_raw' - sparse(diagm(diag(hess_raw)))
    @fact hess_sparse --> roughly([0.0 1.0 0.0; 1.0 0.0 0.0; 0.0 0.0 6.0])

    # Initialize again
    MathProgBase.initialize(d, [:Hess])
    V = zeros(length(I))
    MathProgBase.eval_hesslag(d, V, m.colVal, 1.0, Float64[])
    hess_raw = sparse(I,J,V)
    hess_sparse = hess_raw + hess_raw' - sparse(diagm(diag(hess_raw)))
    @fact hess_sparse --> roughly([0.0 1.0 0.0; 1.0 0.0 0.0; 0.0 0.0 6.0])
end

facts("[nonlinear] Hess-vec through MPB") do
    m = Model()
    @variable(m, a, start = 1)
    @variable(m, b, start = 2)
    @variable(m, c, start = 3)

    @NLobjective(m, Min, a*b + c^2)
    @constraint(m, c*b <= 1)
    @NLconstraint(m, a^2/2 <= 1)
    d = JuMP.NLPEvaluator(m)
    MathProgBase.initialize(d, [:HessVec])
    h = ones(3) # test that input values are overwritten
    v = [2.4,3.5,1.2]
    MathProgBase.eval_hesslag_prod(d, h, m.colVal, v, 1.0, [2.0,3.0])
    correct = [3.0 1.0 0.0; 1.0 0.0 2.0; 0.0 2.0 2.0]*v
    @fact h --> roughly(correct)
end

facts("[nonlinear] NaN corner case (#695)") do

    m = Model()
    x0 = 0.0
    y0 = 0.0

    @variable(m, x >= -1, start = 1.0)
    @variable(m, y, start = 2.0)

    @NLobjective(m, Min, (x - x0) /(sqrt(y0) + sqrt(y)))

    d = JuMP.NLPEvaluator(m)
    MathProgBase.initialize(d, [:HessVec])
    h = ones(2)
    v = [2.4,3.5]
    MathProgBase.eval_hesslag_prod(d, h, m.colVal, v, 1.0, Float64[])
    correct = [0.0 -1/(2*2^(3/2)); -1/(2*2^(3/2)) 3/(4*2^(5/2))]*v
    @fact h --> roughly(correct)
end

if length(convex_nlp_solvers) > 0
    facts("[nonlinear] Error on NLP resolve") do
        m = Model(solver=convex_nlp_solvers[1])
        @variable(m, x, start = 1)
        @NLobjective(m, Min, x^2)
        solve(m)
        setvalue(x, 2)
        @fact_throws ErrorException solve(m)
        EnableNLPResolve()
        status = solve(m)
        @fact status --> :Optimal
    end
end

mysquare(x) = x^2
function myf(x,y)
    return (x-1)^2+(y-2)^2
end

if length(convex_nlp_solvers) > 0
    facts("[nonlinear] User-defined functions") do
        JuMP.register(:myf, 2, myf, autodiff=true)
        JuMP.register(:myf_2, 2, myf, (g,x,y) -> (g[1] = 2(x-1); g[2] = 2(y-2)))
        JuMP.register(:mysquare, 1, mysquare, autodiff=true)
        JuMP.register(:mysquare_2, 1, mysquare, x-> 2x, autodiff=true)
        JuMP.register(:mysquare_3, 1, mysquare, x-> 2x, x -> 2.0)

        m = Model(solver=convex_nlp_solvers[1])

        @variable(m, x[1:2] >= 0.5)
        @NLobjective(m, Min, myf(x[1],mysquare(x[2])))

        d = JuMP.NLPEvaluator(m)
        MathProgBase.initialize(d, [:Grad])
        gradout = zeros(2)
        xval = [1,sqrt(2.0)]
        @fact MathProgBase.eval_f(d, xval) --> roughly(0.0,1e-10)
        MathProgBase.eval_grad_f(d, gradout, xval)
        @fact gradout --> roughly([0.0,0.0],1e-10)

        @fact solve(m) --> :Optimal

        @fact getvalue(x) --> roughly(xval)

        @NLobjective(m, Min, myf_2(x[1],mysquare_2(x[2])))

        d = JuMP.NLPEvaluator(m)
        MathProgBase.initialize(d, [:Grad])
        gradout = zeros(2)
        xval = [1,sqrt(2.0)]
        @fact MathProgBase.eval_f(d, xval) --> roughly(0.0,1e-10)
        MathProgBase.eval_grad_f(d, gradout, xval)
        @fact gradout --> roughly([0.0,0.0],1e-10)
        setvalue(x[1],0.5)
        setvalue(x[2],0.5)
        @fact solve(m) --> :Optimal

        @fact getvalue(x) --> roughly(xval)

        # Test just univariate functions because hessians are disabled
        # if any multivariate functions are present.
        @NLobjective(m, Min, mysquare(x[1]-1) + mysquare_2(x[2]-2) + mysquare_3(x[1]))
        @fact solve(m) --> :Optimal
        @fact getvalue(x) --> roughly([0.5,2.0],1e-4)

    end
end

#  Copyright 2017, Iain Dunning, Joey Huchette, Miles Lubin, and contributors
#  This Source Code Form is subject to the terms of the Mozilla Public
#  License, v. 2.0. If a copy of the MPL was not distributed with this
#  file, You can obtain one at http://mozilla.org/MPL/2.0/.
#############################################################################
# JuMP
# An algebraic modeling langauge for Julia
# See http://github.com/JuliaOpt/JuMP.jl
#############################################################################
# test/nonlinear.jl
# Test general nonlinear
#############################################################################
using Compat, Compat.LinearAlgebra, Compat.SparseArrays, Compat.Test
using MathProgBase, JuMP
# If solvers not loaded, load them (i.e running just these tests)
!isdefined(@__MODULE__, :nlp_solvers) && include("solvers.jl")

mutable struct DummyNLPSolver <: MathProgBase.AbstractMathProgSolver
end
mutable struct DummyNLPModel <: MathProgBase.AbstractNonlinearModel
end


@testset "Nonlinear" begin

    @testset "getvalue on arrays" begin
        m = Model()
        @variable(m, x, start = π/2)
        @NLexpression(m, f1, sin(x))
        @NLexpression(m, f2, sin(2x))
        @NLexpression(m, f3[i=1:2], sin(i*x))

        @test isapprox(getvalue(f1), 1, atol=1e-5)
        @test isapprox(getvalue(f2), 0, atol=1e-5)

        @test getvalue([f1, f2]) == getvalue(f3)

        v = [1.0, 2.0]
        @NLparameter(m, vparam[i=1:2] == v[i])
        @test getvalue(vparam) == v
        v[1] = 3.0
        setvalue(vparam, v)
        @test getvalue(vparam) == v
    end

    @testset "HS071 solves correctly with $nlp_solver" for nlp_solver in nlp_solvers
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
        @NLconstraint(m, sum(x[i]^2 for i=1:4) == 40)
        @test JuMP.numnlconstr(m) == 2
        @test MathProgBase.numconstr(m) == 2
        status = solve(m)

        @test status == :Optimal
        @test isapprox(getvalue(x),
            [1.000000, 4.742999, 3.821150, 1.379408], atol=1e-3)
    end

    @testset "HS071 solves correctly (no macros) with $nlp_solver" for nlp_solver in nlp_solvers
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
        JuMP.setNLobjective(m, :Min, :($(x[1])*$(x[4])*($(x[1])+$(x[2])+$(x[3])) + $(x[3])))
        JuMP.addNLconstraint(m, :($(x[1])*$(x[2])*$(x[3])*$(x[4]) >= 25))
        JuMP.addNLconstraint(m, :($(x[1])^2+$(x[2])^2+$(x[3])^2+$(x[4])^2 == 40))
        @test MathProgBase.numconstr(m) == 2
        @test_throws ErrorException JuMP.addNLconstraint(m, :(x[1]^2+x[2]^2+x[3]^2+x[4]^2 == 40))
        status = solve(m)

        @test status == :Optimal
        @test isapprox(getvalue(x),
            [1.000000, 4.742999, 3.821150, 1.379408], atol=1e-3)
    end

    @testset "HS071 solves correctly (epigraph) with $nlp_solver" for nlp_solver in nlp_solvers
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

        @test status == :Optimal
        @test isapprox(getvalue(x),
            [1.000000, 4.742999, 3.821150, 1.379408], atol=1e-3)
    end

    @testset "ifelse with $nlp_solver" for nlp_solver in nlp_solvers
        (occursin("OsilSolver", "$(typeof(nlp_solver))") || occursin("NLoptSolver", "$(typeof(nlp_solver))") || occursin("BaronSolver", "$(typeof(nlp_solver))")) && continue
        m = Model(solver=nlp_solver)
        @variable(m, x, start = 2)
        # minimizer at smooth point, solvers should be okay
        @NLobjective(m, Min, ifelse( x <= 1, x^2, x) )
        status = solve(m)

        @test status == :Optimal
        @test isapprox(getvalue(x), 0.0, atol=1e-5)
    end

    @testset "Accepting fixed variables with $nlp_solver (simplify = $simplify)" for nlp_solver in convex_nlp_solvers, simplify in [true,false]
        m = Model(solver=nlp_solver, simplify_nonlinear_expressions=simplify)
        @variable(m, x == 0)
        @variable(m, y ≥ 0)
        @objective(m, Min, y)
        @NLconstraint(m, y ≥ x^2)
        for α in 1:4
            JuMP.fix(x, α)
            solve(m)
            @test isapprox(getvalue(y), α^2, atol=1e-6)
        end
    end

    @testset "QP solve through NL pathway with $nlp_solver" for nlp_solver in nlp_solvers
        # Solve a problem with quadratic objective with linear
        # constraints, but force it to use the nonlinear code.
        m = Model(solver=nlp_solver)
        @variable(m, 0.5 <= x <=  2)
        @variable(m, 0.0 <= y <= 30)
        @NLparameter(m, param == 1.0)
        @objective(m, Min, (x+y)^2)
        @NLconstraint(m, x + y >= param)
        status = solve(m)

        @test status == :Optimal
        @test isapprox(getobjectivevalue(m), 1.0, atol=1e-6)
        @NLexpression(m, lhs, x+y)
        @test isapprox(getvalue(x)+getvalue(y), 1.0, atol=1e-6)
        @test isapprox(getvalue(lhs), 1.0, atol=1e-6)

        setvalue(param,10)
        @test m.internalModelLoaded == true
        status = solve(m)
        @test isapprox(getobjectivevalue(m), 10.0^2, atol=1e-6)
        @test isapprox(getvalue(x)+getvalue(y), 10.0, atol=1e-6)
        @test isapprox(getvalue(lhs), 10.0, atol=1e-6)
    end


    @testset "Quad con solve through NL pathway" for nlp_solver in convex_nlp_solvers
        # Solve a problem with linear objective with quadratic
        # constraints, but force it to use the nonlinear code.
        m = Model(solver=nlp_solver)
        @variable(m, -2 <= x <= 2)
        @variable(m, -2 <= y <= 2)
        @NLobjective(m, Min, x - y)
        @constraint(m, x + x^2 + x*y + y^2 <= 1)
        status = solve(m)

        @test status == :Optimal
        @test isapprox(getobjectivevalue(m), -1-4/sqrt(3), atol=1e-6)
        @test isapprox(getvalue(x) + getvalue(y), -1/3, atol=1e-3)
    end

    @testset "Resolve with parameter with $nlp_solver (simplify = $simplify)" for nlp_solver in convex_nlp_solvers, simplify in [true,false]
        m = Model(solver=nlp_solver, simplify_nonlinear_expressions=simplify)
        @variable(m, z)
        @NLparameter(m, x == 1.0)
        @NLobjective(m, Min, (z-x)^2)
        status = solve(m)
        @test status == :Optimal
        @test isapprox(getvalue(z), 1.0, atol=1e-3)

        setvalue(x, 5.0)
        status = solve(m)
        @test status == :Optimal
        @test isapprox(getvalue(z), 5.0, atol=1e-3)
    end

    @testset "Two-sided nonlinear constraints with $nlp_solver" for nlp_solver in convex_nlp_solvers
        m = Model(solver=nlp_solver)
        @variable(m, x)
        @NLobjective(m, Max, x)
        l = -1
        u = 1
        @NLconstraint(m, l <= x <= u)
        status = solve(m)

        @test status == :Optimal
        @test isapprox(getobjectivevalue(m), u, atol=1e-6)

        @NLobjective(m, Min, x)
        status = solve(m)

        @test status == :Optimal
        @test isapprox(getobjectivevalue(m), l, atol=1e-6)
    end

    @testset "Two-sided nonlinear constraints (no macros) with $nlp_solver" for nlp_solver in convex_nlp_solvers
        m = Model(solver=nlp_solver)
        @variable(m, x)
        JuMP.setNLobjective(m, :Max, x)
        l = -1
        u = 1
        JuMP.addNLconstraint(m, :($l <= $x <= $u))
        status = solve(m)

        @test status == :Optimal
        @test isapprox(getobjectivevalue(m), u, atol=1e-6)

        JuMP.setNLobjective(m, :Min, x)
        status = solve(m)

        @test status == :Optimal
        @test isapprox(getobjectivevalue(m), l, atol=1e-6)
    end

    @testset "Quadratic equality constraints with $nlp_solver" for nlp_solver in nlp_solvers
        m = Model(solver=nlp_solver)
        @variable(m, 0 <= x[1:2] <= 1)
        @constraint(m, x[1]^2 + x[2]^2 == 1/2)
        @NLobjective(m, Max, x[1] - x[2])
        status = solve(m)

        @test status == :Optimal
        @test isapprox(getvalue(x), [sqrt(1/2), 0], atol=1e-6)
    end

    if ipt
        @testset "Passing starting solutions through QP pathway with Ipopt" begin
            # https://discourse.julialang.org/t/create-quadratic-objective-will-objective-and-nlobjective-lead-to-different-solution-using-ipopt/1666
            m = Model(solver=Ipopt.IpoptSolver(print_level=0))
            @variable(m, 0<= x1 <= 1, start=1)
            @variable(m, 0<= x2 <=1, start=1)
            @variable(m, 0<= x3 <=1)
            @variable(m, 0<= x4 <=1, start=1)
            @variable(m, 0<= x5 <=1)
            @constraint(m, 20*x1 + 12*x2 + 11*x3 + 7*x4 + 4*x5 <= 40)
            @objective(m, Min, 42*x1 - 0.5*(100*x1*x1 + 100*x2*x2 + 100*x3*x3 + 100*x4*x4 + 100*x5*x5) + 44*x2 + 45*x3 + 47*x4 + 47.5*x5)
            status=solve(m)

            @test status == :Optimal
            @test isapprox(getobjectivevalue(m), -17, atol=1e-4)
        end
    end

    @testset "Mixed integer nonlinear problems with $minlp_solver" for minlp_solver in minlp_solvers
        # Solve test problem 1 (Synthesis of processing system) in
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

        @test status == :Optimal
        @test isapprox(getobjectivevalue(m), 6.00976, atol=1e-5)
        @test isapprox(getvalue(x), [1.30098, 0.0, 1.0], atol=1e-5)
        @test isapprox(getvalue(y)[:], [0.0, 1.0, 0.0], atol=1e-5)
    end

    @testset "Continuous relaxation of minlp test problem with $nlp_solver" for nlp_solver in nlp_solvers
        # Solve continuous relaxation of test problem 1 (Synthesis of processing system) in
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

        @test status == :Optimal
        @test isapprox(getobjectivevalue(m), 0.7593, atol=5e-5)
        @test isapprox(getvalue(x), [1.1465, 0.54645, 1.0], atol=2e-4)
        @test isapprox(getvalue(y)[:], [0.2732, 0.3, 0.0], atol=2e-4)
        @test isapprox(getvalue(z), 1.6, atol=2e-4)
    end

    @testset "Maximization objective with $nlp_solver" for nlp_solver in convex_nlp_solvers
        # Solve a simple problem with a maximization objective
        m = Model(solver=nlp_solver)
        @variable(m, -2 <= x <= 2); setvalue(x, -1.8)
        @variable(m, -2 <= y <= 2); setvalue(y,  1.5)
        @NLobjective(m, Max, y - x)
        @constraint(m, x + x^2 + x*y + y^2 <= 1)

        @test solve(m) == :Optimal
        @test isapprox(getobjectivevalue(m), 1+4/sqrt(3), atol=1e-6)
        @test isapprox(getvalue(x) + getvalue(y), -1/3, atol=1e-3)
    end

    @testset "Maximization objective (embedded expressions) with $nlp_solver (simplify = $simplify)" for nlp_solver in convex_nlp_solvers, simplify in [true,false]
        m = Model(solver=nlp_solver, simplify_nonlinear_expressions=simplify)
        @variable(m, -2 <= x <= 2); setvalue(x, -1.8)
        @variable(m, -2 <= y <= 2); setvalue(y,  1.5)
        @NLobjective(m, Max, y - x)
        @NLexpression(m, quadexpr, x + x^2 + x*y + y^2)
        @NLconstraint(m, quadexpr <= 1)

        @test solve(m) == :Optimal
        @test isapprox(getobjectivevalue(m), 1+4/sqrt(3), atol=1e-6)
        @test isapprox(getvalue(x) + getvalue(y), -1/3, atol=1e-3)
        @test isapprox(getvalue(quadexpr), 1, atol=1e-5)
        @NLexpression(m, quadexpr2, x + x^2 + x*y + y^2)
        @test isapprox(getvalue(quadexpr2), 1, atol=1e-5)
    end


    @testset "Infeasibility detection with $nlp_solver" for nlp_solver in convex_nlp_solvers
        occursin("NLoptSolver", string(typeof(nlp_solver))) && continue
        occursin("MosekSolver", string(typeof(nlp_solver))) && continue
        # (Attempt to) solve an infeasible problem
        m = Model(solver=nlp_solver)
        n = 10
        @variable(m, 0 <= x[i=1:n] <= 1)
        @NLobjective(m, Max, x[n])
        for i in 1:n-1
            @NLconstraint(m, x[i+1]-x[i] == 0.15)
        end
        @test solve(m, suppress_warnings=true) == :Infeasible
    end


    @testset "Unboundedness detection with $nlp_solver" for nlp_solver in convex_nlp_solvers
        occursin("NLoptSolver", string(typeof(nlp_solver))) && continue
        # (Attempt to) solve an unbounded problem
        m = Model(solver=nlp_solver)
        @variable(m, x >= 0)
        @NLobjective(m, Max, x)
        @NLconstraint(m, x >= 5)
        @test solve(m, suppress_warnings=true) == :Unbounded
    end

    @testset "Entropy maximization with $nlp_solver" for nlp_solver in convex_nlp_solvers
        m = Model(solver=nlp_solver)
        N = 3
        @variable(m, x[1:N] >= 0, start = 1)
        @NLexpression(m, entropy[i=1:N], -x[i]*log(x[i]))
        @NLobjective(m, Max, sum(entropy[i] for i = 1:N))
        @constraint(m, sum(x) == 1)

        @test solve(m) == :Optimal
        @test isapprox(getvalue(x), [1/3,1/3,1/3], atol=1e-3)
    end

    @testset "Entropy maximization (reformulation) with $nlp_solver" for nlp_solver in convex_nlp_solvers
        m = Model(solver=nlp_solver)
        idx = [1,2,3,4]
        @variable(m, x[idx] >= 0, start = 1)
        @variable(m, z[1:4], start = 0)
        @NLexpression(m, entropy[i=idx], -x[i]*log(x[i]))
        @NLobjective(m, Max, sum(z[i] for i = 1:2) + sum(z[i]/2 for i=3:4))
        @NLconstraint(m, z_constr1[i=1], z[i] <= entropy[i])
        @NLconstraint(m, z_constr1_dup[i=2], z[i] <= entropy[i]) # duplicate expressions
        @NLconstraint(m, z_constr2[i=3:4], z[i] <= 2*entropy[i])
        @constraint(m, sum(x) == 1)

        @test solve(m) == :Optimal
        @test isapprox([getvalue(x[i]) for i in idx], [1/4,1/4,1/4,1/4], atol=1e-3)
        @test isapprox(getvalue(entropy[1]), -(1/4)*log(1/4), atol=1e-4)
        @NLexpression(m, zexpr[i=1:4], z[i])
        @test isapprox(getvalue(zexpr[1]), -(1/4)*log(1/4), atol=1e-4)
    end

    @testset "Derivatives of x^4, x < 0 with $nlp_solver" for nlp_solver in convex_nlp_solvers
        m = Model(solver=nlp_solver)
        @variable(m, x >= -1, start = -0.5)
        @NLobjective(m, Min, x^4)
        status = solve(m)

        @test status == :Optimal
        @test isapprox(getvalue(x), 0.0, atol=1e-2)
    end

    @testset "Nonlinear duals with $nlp_solver (simplify = $simplify)" for nlp_solver in nlp_solvers, simplify in [true,false]
        applicable(MathProgBase.getconstrduals, MathProgBase.NonlinearModel(nlp_solver)) || continue
        modA = Model(solver=nlp_solver, simplify_nonlinear_expressions=simplify)
        @variable(modA, x >= 0)
        @variable(modA, y <= 5)
        @variable(modA, 2 <= z <= 4)
        @variable(modA, 0 <= r[i=3:6] <= i)
        @NLobjective(modA, Min, -((x + y)/2.0 + 3.0)/3.0 - z - r[3])
        @constraint(modA, cons1, x+y >= 2)
        @constraint(modA, cons2, sum(r[i] for i=3:5) <= (2 - x)/2.0)
        @NLconstraint(modA, cons3, 7.0*y <= z + r[6]/1.9)

        # Getter/setters
        @test MathProgBase.numconstr(modA) == 3
        @test JuMP.numnlconstr(modA) == 1

        # Solution
        @test solve(modA) == :Optimal
        @test isapprox(getobjectivevalue(modA), -5.8446115, atol=1e-6)
        @test isapprox(getvalue(x), 0.9774436, atol=1e-6)
        @test isapprox(getvalue(y), 1.0225563, atol=1e-6)
        @test isapprox(getvalue(z), 4.0, atol=1e-6)
        @test isapprox(getvalue(r)[3], 0.5112781, atol=1e-6)
        @test isapprox(getvalue(r)[4], 0.0, atol=1e-6)
        @test isapprox(getvalue(r)[5], 0.0, atol=1e-6)
        @test isapprox(getvalue(r)[6], 6.0, atol=1e-6)

        # Reduced costs
        @test isapprox(getdual(x), 0.0, atol=1e-6)
        @test isapprox(getdual(y), 0.0, atol=1e-6)
        @test isapprox(getdual(z), -1.0714286, atol=1e-6)
        @test isapprox(getdual(r)[3], 0.0, atol=1e-6)
        @test isapprox(getdual(r)[4], 1.0, atol=1e-6)
        @test isapprox(getdual(r)[5], 1.0, atol=1e-6)
        @test isapprox(getdual(r)[6], -0.03759398, atol=1e-6)

        # Row duals
        @test isapprox(getdual(cons1), 0.333333, atol=1e-6)
        @test isapprox(getdual(cons2), -1.0, atol=1e-6)
        @test isapprox(getdual(cons3), -0.0714286, atol=1e-6)
    end

    @testset "Nonlinear duals (Max) with $nlp_solver" for nlp_solver in nlp_solvers
        applicable(MathProgBase.getconstrduals, MathProgBase.NonlinearModel(nlp_solver)) || continue
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
        @test solve(modA) == :Optimal
        @test isapprox(getobjectivevalue(modA), 5.8446115, atol=1e-6)
        @test isapprox(getvalue(x), 0.9774436, atol=1e-6)
        @test isapprox(getvalue(y), 1.0225563, atol=1e-6)
        @test isapprox(getvalue(z), 4.0, atol=1e-6)
        @test isapprox(getvalue(r)[3], 0.5112781, atol=1e-6)
        @test isapprox(getvalue(r)[4], 0.0, atol=1e-6)
        @test isapprox(getvalue(r)[5], 0.0, atol=1e-6)
        @test isapprox(getvalue(r)[6], 6.0, atol=1e-6)

        # Reduced costs
        @test isapprox(getdual(x), 0.0, atol=1e-6)
        @test isapprox(getdual(y), 0.0, atol=1e-6)
        @test isapprox(getdual(z), 1.0714286, atol=1e-6)
        @test isapprox(getdual(r)[3], 0.0, atol=1e-6)
        @test isapprox(getdual(r)[4], -1.0, atol=1e-6)
        @test isapprox(getdual(r)[5], -1.0, atol=1e-6)
        @test isapprox(getdual(r)[6], 0.03759398, atol=1e-6)

        # Row duals
        @test isapprox(getdual(cons1), -0.333333, atol=1e-6)
        @test isapprox(getdual(cons2), 1.0, atol=1e-6)
        @test isapprox(getdual(cons3), 0.0714286, atol=1e-6)
    end

    @testset "Changing objectives with $nlp_solver" for nlp_solver in nlp_solvers
        m = Model(solver=nlp_solver)
        @variable(m, x >= 0)
        @variable(m, y >= 0)
        @objective(m, Max, x+y)
        @NLconstraint(m, x+2y <= 1)
        @test solve(m) == :Optimal
        @test isapprox(getvalue(x), 1.0, atol=1e-4)
        @test isapprox(getvalue(y), 0.0, atol=1e-4)
        @test isapprox(getobjectivevalue(m), 1.0, atol=1e-4)

        @objective(m, Max, 2x+y)
        @test solve(m) == :Optimal
        @test isapprox(getvalue(x), 1.0, atol=1e-4)
        @test isapprox(getvalue(y), 0.0, atol=1e-4)
        @test isapprox(getobjectivevalue(m), 2.0, atol=1e-4)
    end

    @testset "Hessian chunking code with $nlp_solver" for nlp_solver in nlp_solvers
        m = Model(solver=nlp_solver)
        @variable(m, x[1:18] >= 1, start = 1.2)
        @NLobjective(m, Min, prod(x[i] for i=1:18))
        @test solve(m) == :Optimal
        @test isapprox(getvalue(x), ones(18), atol=1e-4)
    end

    #############################################################################
    # Test that output is produced in correct MPB form
    MathProgBase.NonlinearModel(s::DummyNLPSolver) = DummyNLPModel()
    function MathProgBase.loadproblem!(m::DummyNLPModel, numVar, numConstr, x_l, x_u, g_lb, g_ub, sense, d::MathProgBase.AbstractNLPEvaluator)
        MathProgBase.initialize(d, [:ExprGraph])
        objexpr = MathProgBase.obj_expr(d)
        @testset "NL MPB interface ($objexpr)" begin
            @test objexpr == :(x[1]^x[2]) || objexpr == :(-1.0*x[1]+1.0*x[2])
            @test MathProgBase.isconstrlinear(d,1)
            @test MathProgBase.isconstrlinear(d,3)
            @test MathProgBase.constr_expr(d,1) == :(2.0*x[1] + 1.0*x[2] <= 1.0)
            @test MathProgBase.constr_expr(d,2) == :(2.0*x[1] + 1.0*x[2] <= -0.0)
            @test MathProgBase.constr_expr(d,3) == :(-5.0 <= 2.0*x[1] + 1.0*x[2] <= 5.0)
            if numConstr > 3
                @test MathProgBase.constr_expr(d,4) == :(2.0*x[1]*x[1] + 1.0*x[2] + -2.0 >= 0)
                @test MathProgBase.constr_expr(d,5) == :(sin(x[1]) * cos(x[2]) - 5 == 0.0)
                @test MathProgBase.constr_expr(d,6) == :(1.0*x[1]^2 - 1.0 == 0.0)
                @test MathProgBase.constr_expr(d,7) == :(2.0*x[1]^2 - 2.0 == 0.0)
                @test MathProgBase.constr_expr(d,8) == :(-0.5 <= sin(x[1]) <= 0.5)
                @test MathProgBase.constr_expr(d,9) == :(ψ(x[1]) + Ψ(x[1],x[2]) - 3.0 <= 0.0)
            end
        end
    end
    MathProgBase.setwarmstart!(m::DummyNLPModel,x) = nothing
    MathProgBase.optimize!(m::DummyNLPModel) = nothing
    MathProgBase.status(m::DummyNLPModel) = :Optimal
    MathProgBase.getobjval(m::DummyNLPModel) = NaN
    MathProgBase.getsolution(m::DummyNLPModel) = [1.0,1.0]
    MathProgBase.setvartype!(m::DummyNLPModel,vartype) = @test any(vartype .== :Fixed) == false

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
        @test lb == [-Inf,-Inf,-5.0]
        @test ub == [1.0,-0.0,5.0]

        @constraint(m, 2x^2+y >= 2)
        ψ(x) = 1
        Ψ(x,y) = 2
        JuMP.register(m, :ψ, 1, ψ, autodiff=true)
        JuMP.register(m, :Ψ, 2, Ψ, autodiff=true)

        @NLconstraint(m, sin(x)*cos(y) == 5)
        @NLconstraint(m, nlconstr[i=1:2], i*x^2 == i)
        @NLconstraint(m, -0.5 <= sin(x) <= 0.5)
        @NLconstraint(m, ψ(x) + Ψ(x,y) <= 3)
        solve(m)

        @NLobjective(m, Min, x^y)
        solve(m)

        lb,ub = JuMP.constraintbounds(m)
        @test lb == [-Inf,-Inf,-5.0,0.0,0.0,0.0,0.0,-0.5,-Inf]
        @test ub == [1.0,-0.0,5.0,Inf,0.0,0.0,0.0,0.5,0.0]
    end
    test_nl_mpb()

    @testset "Simplified expression graphs" begin
        m = Model(simplify_nonlinear_expressions=true)
        # this is expert behavior, expression simplification is experimental
        @variable(m, x == 2)
        @variable(m, y)
        @NLobjective(m, Min, x^2 + y^2)
        @NLexpression(m, ex, exp(x))
        @NLconstraint(m, ex - y == 0)
        @NLconstraint(m, ex + 1 == 0)

        d = JuMP.NLPEvaluator(m)
        MathProgBase.initialize(d, [:ExprGraph])
        @test MathProgBase.obj_expr(d) == :(4.0 + x[2] ^ 2.0)
        @test MathProgBase.constr_expr(d,1) == :(($(exp(2)) - x[2]) - 0.0 == 0.0)
        @test MathProgBase.constr_expr(d,2) == :($(exp(2) + 1) == 0.0)
    end

    @testset "Expression graph for linear problem" begin
        m = Model()
        @variable(m, x)
        @constraint(m, 0 <= x <= 1)
        @objective(m, Max, x)
        d = JuMP.NLPEvaluator(m)
        MathProgBase.initialize(d, [:ExprGraph])
        @test MathProgBase.obj_expr(d) == :(+(1.0 * x[1]))
    end

    @testset "Expression graph for ifelse" begin
        m = Model()
        @variable(m, x, start = 2)
        @NLobjective(m, Min, ifelse( x <= 1, x^2, x) )
        d = JuMP.NLPEvaluator(m)
        MathProgBase.initialize(d, [:ExprGraph])
        @test MathProgBase.obj_expr(d) == :(ifelse( x[1] <= 1, x[1]^2, x[1]))
    end

    @testset "Expression graphs for corner cases" begin
        m = Model()
        @variable(m, x, start = 2)
        @constraint(m, 0 <= 1)
        @NLconstraint(m, x <= sum(0 for i in []) + prod(1 for i in []))
        d = JuMP.NLPEvaluator(m)
        MathProgBase.initialize(d, [:ExprGraph])
        @test MathProgBase.constr_expr(d,1) == :(0 <= 1.0)
        @test MathProgBase.constr_expr(d,2) == :(x[1] - (0 + 1) <= 0.0)
    end

    @testset "Hessians through MPB" begin
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
        hess_sparse = hess_raw + hess_raw' - sparse(Diagonal(diag(hess_raw)))
        @test isapprox(hess_sparse, [0.0 1.0 0.0; 1.0 0.0 0.0; 0.0 0.0 2.0])

        # test view
        x = zeros(length(m.colVal) + 10)
        x[6:length(m.colVal)+5] = m.colVal
        xv = @view x[6:length(m.colVal)+5]
        Vv = zeros(length(I))
        MathProgBase.eval_hesslag(d, Vv, xv, 1.0, Float64[])
        @test Vv == V

        V2 = zeros(length(I) + 10)
        Vv = @view V2[6:length(I)+5]
        MathProgBase.eval_hesslag(d, Vv, m.colVal, 1.0, Float64[])
        @test Vv == V

        # make sure we don't get NaNs in this case
        @NLobjective(m, Min, a * b + 3*c^2)
        d = JuMP.NLPEvaluator(m)
        MathProgBase.initialize(d, [:Hess])
        setvalue(c, -1.0)
        V = zeros(length(I))
        MathProgBase.eval_hesslag(d, V, m.colVal, 1.0, Float64[])
        hess_raw = sparse(I,J,V)
        hess_sparse = hess_raw + hess_raw' - sparse(Diagonal(diag(hess_raw)))
        @test isapprox(hess_sparse, [0.0 1.0 0.0; 1.0 0.0 0.0; 0.0 0.0 6.0])

        # Initialize again
        MathProgBase.initialize(d, [:Hess])
        V = zeros(length(I))
        MathProgBase.eval_hesslag(d, V, m.colVal, 1.0, Float64[])
        hess_raw = sparse(I,J,V)
        hess_sparse = hess_raw + hess_raw' - sparse(Diagonal(diag(hess_raw)))
        @test isapprox(hess_sparse, [0.0 1.0 0.0; 1.0 0.0 0.0; 0.0 0.0 6.0])
    end

    @testset "Hess-vec through MPB" begin
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
        @test isapprox(h, correct)

        # test view
        h2 = zeros(10)
        hv = @view h2[3:5]
        MathProgBase.eval_hesslag_prod(d, hv, m.colVal, v, 1.0, [2.0,3.0])
        @test hv == h

        x = zeros(length(m.colVal) + 10)
        x[6:length(m.colVal)+5] = m.colVal
        xv = @view x[6:length(m.colVal)+5]
        vv = zeros(13)
        vv[6:8] = v
        vw = @view vv[6:8]
        h2 = zeros(3)
        MathProgBase.eval_hesslag_prod(d, h2, xv, vw, 1.0, [2.0,3.0])
        @test h2 == h
    end

    @testset "Hess-vec through MPB with subexpressions" begin
        m = Model()
        @variable(m, a, start = 1)
        @variable(m, b, start = 2)
        @variable(m, c, start = 3)

        @NLexpression(m, ab, a*b)
        @NLobjective(m, Min, ab + c^2)
        @constraint(m, c*b <= 1)
        @NLconstraint(m, a^2/2 <= 1)
        d = JuMP.NLPEvaluator(m)
        MathProgBase.initialize(d, [:HessVec])
        h = ones(3) # test that input values are overwritten
        v = [2.4,3.5,1.2]
        MathProgBase.eval_hesslag_prod(d, h, m.colVal, v, 1.0, [2.0,3.0])
        correct = [3.0 1.0 0.0; 1.0 0.0 2.0; 0.0 2.0 2.0]*v
        @test isapprox(h, correct)
    end

    @testset "NaN corner case (#695)" begin

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
        @test isapprox(h, correct)
    end

    mysquare(x) = x^2
    function myf(x,y)
        return (x-1)^2+(y-2)^2
    end

    if length(convex_nlp_solvers) > 0
        @testset "User-defined functions" begin
            m = Model(solver=convex_nlp_solvers[1])

            JuMP.register(m, :myf, 2, myf, autodiff=true)
            JuMP.register(m, :myf_2, 2, myf, (g,x,y) -> (g[1] = 2(x-1); g[2] = 2(y-2)))
            JuMP.register(m, :mysquare, 1, mysquare, autodiff=true)
            JuMP.register(m, :mysquare_2, 1, mysquare, x-> 2x, autodiff=true)
            JuMP.register(m, :mysquare_3, 1, mysquare, x-> 2x, x -> 2.0)

            @variable(m, x[1:2] >= 0.5)
            @NLobjective(m, Min, myf(x[1],mysquare(x[2])))

            d = JuMP.NLPEvaluator(m)
            MathProgBase.initialize(d, [:Grad])
            gradout = zeros(2)
            xval = [1,sqrt(2.0)]
            @test isapprox(MathProgBase.eval_f(d, xval), 0.0, atol=1e-10)
            MathProgBase.eval_grad_f(d, gradout, xval)
            @test isapprox(gradout, [0.0,0.0], atol=1e-10)

            @test solve(m) == :Optimal

            @test isapprox(getvalue(x), xval)

            @NLobjective(m, Min, myf_2(x[1],mysquare_2(x[2])))

            d = JuMP.NLPEvaluator(m)
            MathProgBase.initialize(d, [:Grad])
            gradout = zeros(2)
            xval = [1,sqrt(2.0)]
            @test isapprox(MathProgBase.eval_f(d, xval), 0.0, atol=1e-10)
            MathProgBase.eval_grad_f(d, gradout, xval)
            @test isapprox(gradout, [0.0,0.0], atol=1e-10)
            setvalue(x[1],0.5)
            setvalue(x[2],0.5)
            @test solve(m) == :Optimal

            @test isapprox(getvalue(x), xval)

            # Test just univariate functions because hessians are disabled
            # if any multivariate functions are present.
            @NLobjective(m, Min, mysquare(x[1]-1) + mysquare_2(x[2]-2) + mysquare_3(x[1]))
            @test solve(m) == :Optimal
            @test isapprox(getvalue(x), [0.5,2.0], atol=1e-4)

            # Test #927
            m = Model(solver=convex_nlp_solvers[1])
            JuMP.register(m, :myf, 2, myf, autodiff=true)
            @variable(m, x)
            @NLobjective(m, Min, myf(x,x))
            @test solve(m) == :Optimal
            @test isapprox(getvalue(x), 1.5, atol=1e-4)

        end

        @testset "Anonymous nonlinear expression" begin
            m = Model(solver=convex_nlp_solvers[1])
            @variable(m, -1 <= x <= 1)
            obj = @NLexpression(m, x^4 + 1)
            @NLobjective(m, Min, obj)
            @test solve(m) == :Optimal
            @test isapprox(getvalue(x), 0, atol=1e-4)
        end
    end
end

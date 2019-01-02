#  Copyright 2017, Iain Dunning, Joey Huchette, Miles Lubin, and contributors
#  This Source Code Form is subject to the terms of the Mozilla Public
#  License, v. 2.0. If a copy of the MPL was not distributed with this
#  file, You can obtain one at http://mozilla.org/MPL/2.0/.
#############################################################################
# JuMP
# An algebraic modeling language for Julia
# See http://github.com/JuliaOpt/JuMP.jl
#############################################################################
# test/nonlinear.jl
# Test general nonlinear
#############################################################################
using JuMP
using Base.Test

# If solvers not loaded, load them (i.e running just these tests)
!isdefined(:nlp_solvers) && include("solvers.jl")

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
        occursin("NLoptSolver",string(typeof(nlp_solver))) && continue
        occursin("MosekSolver",string(typeof(nlp_solver))) && continue
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
        occursin("NLoptSolver",string(typeof(nlp_solver))) && continue
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

    @testset "Expression graph for linear problem" begin
        m = Model()
        @variable(m, x)
        @constraint(m, 0 <= x <= 1)
        @objective(m, Max, x)
        d = JuMP.NLPEvaluator(m)
        MathProgBase.initialize(d, [:ExprGraph])
        @test MathProgBase.obj_expr(d) == :(+(1.0 * x[1]))
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

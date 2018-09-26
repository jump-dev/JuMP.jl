#  Copyright 2017, Iain Dunning, Joey Huchette, Miles Lubin, and contributors
#  This Source Code Form is subject to the terms of the Mozilla Public
#  License, v. 2.0. If a copy of the MPL was not distributed with this
#  file, You can obtain one at http://mozilla.org/MPL/2.0/.
#############################################################################
# JuMP
# An algebraic modeling langauge for Julia
# See http://github.com/JuliaOpt/JuMP.jl
#############################################################################
# test/probmod.jl
# Testing problem modification
# Must be run as part of runtests.jl, as it needs a list of solvers.
#############################################################################
using JuMP, Compat.Test, Compat
# If solvers not loaded, load them (i.e running just these tests)
!isdefined(@__MODULE__, :lp_solvers) && include("solvers.jl")

const TOL = 1e-4

@testset "Problem modification" begin
    @testset "Problem modification basics with $solver" for solver in lp_solvers
        # max 1.1x + 1.0y
        # st     x +    y <= 3
        #     0 <= x <= 3
        #     1 <= y <= 3
        # x* = 2, y* = 1
        m = Model(solver=solver)
        @variable(m, 0 <= x <= 3)
        @variable(m, 1 <= y <= 3)
        @objective(m, :Max, 1.1x + 1.0y)
        maincon = @constraint(m, x + y <= 3)
        @test solve(m) == :Optimal
        @test m.internalModelLoaded == true
        @test isapprox(getvalue(x), 2.0, atol=TOL)
        @test isapprox(getvalue(y), 1.0, atol=TOL)

        # Test adding a variable
        # max 1.1x + 1.0y + 100.0z
        # st     x +    y +      z <= 3
        #     0 <= x <= 3
        #     1 <= y <= 3
        #     0 <= z <= 5
        # x* = 0, y* = 1, z* = 2
        @variable(m, 0 <= z <= 5, objective=100.0, inconstraints=[maincon], coefficients=[1.0])
        @test solve(m) == :Optimal
        @test isapprox(getvalue(x), 0.0, atol=TOL)
        @test isapprox(getvalue(y), 1.0, atol=TOL)
        @test isapprox(getvalue(z), 2.0, atol=TOL)
        @test isapprox(getobjectivevalue(m), 201.0, atol=TOL)


        # Test changing bounds
        # max 1.1x + 1.0y + 100.0z
        # st     x +    y +      z <= 3
        #     0 <= x <= 3
        #     0 <= y <= 3
        #     0 <= z <= 2
        # x* = 1, y* = 0, z* = 2
        setlowerbound(y, 0.0)
        setupperbound(z, 2.0)
        @test solve(m) == :Optimal
        @test isapprox(getvalue(x), 1.0, atol=TOL)
        @test isapprox(getvalue(y), 0.0, atol=TOL)
        @test isapprox(getvalue(z), 2.0, atol=TOL)
        m.internalModelLoaded = false

        # Test fixing variable
        # max 1.1x + 1.0y + 100.0z
        # st     x +    y +      z <= 3
        #     0.5 <= x <= 0.5
        #     0 <= y <= 3
        #     0 <= z <= 2
        # x* = 0.5, y* = 0.5, z* = 2
        JuMP.fix(x, 0.5)
        @test solve(m) == :Optimal
        @test isapprox(getvalue(x), 0.5, atol=TOL)
        @test isapprox(getvalue(y), 0.5, atol=TOL)
        @test isapprox(getvalue(z), 2.0, atol=TOL)
        m.internalModelLoaded = false
    end


    @testset "Problem modification part two with $solver" for solver in ip_solvers
        # max 1.1x + 1.0y + 100.0z
        # st     x +    y +      z <= 3
        #     0 <= x <= 3
        #     1 <= y <= 3
        #     0 <= z <= 0
        m = Model(solver=solver)
        @variable(m, 0 <= x <= 3)
        @variable(m, 0 <= y <= 3)
        @variable(m, 0 <= z <= 0)
        @objective(m, :Max, 1.1x + 1.0y + 100.0z)
        maincon = @constraint(m, x + y + z <= 3)
        @test solve(m) == :Optimal

        # Test changing problem type
        # max 1.1x + 1.0y + 100.0z
        # st     x +    y +      z <= 3
        #     0 <= x <= 3
        #     0 <= y <= 3
        #     0 <= z <= 1.5, Integer
        # x* = 2, y* = 0, z* = 1
        setupperbound(z, 1.5)
        m.colCat[3] = :Int
        @test solve(m) == :Optimal
        @test isapprox(getvalue(x), 2.0, atol=TOL)
        @test isapprox(getvalue(y), 0.0, atol=TOL)
        @test isapprox(getvalue(z), 1.0, atol=TOL)

        # Test changing constraint bound (<= constraint)
        # max 1.1x + 1.0y + 100.0z
        # st     x +    y +      z <= 2
        #     0 <= x <= 3
        #     0 <= y <= 3
        #     0 <= z <= 1.5, Integer
        # x* = 1, y* = 0, z* = 1
        JuMP.setRHS(maincon, 2.0)
        @test solve(m) == :Optimal
        @test isapprox(getvalue(x), 1.0, atol=TOL)
        @test isapprox(getvalue(y), 0.0, atol=TOL)
        @test isapprox(getvalue(z), 1.0, atol=TOL)

        # Test adding a constraint
        # max 1.1x + 1.0y + 100.0z
        # st     x +    y +      z <= 2
        #        x        +      z <= 0
        #     0 <= x <= 3
        #     0 <= y <= 3
        #     0 <= z <= 1.5, Integer
        # x* = 0, y* = 2, z* = 0
        xz0ref = @constraint(m, x + z <=0)
        @test solve(m) == :Optimal
        @test isapprox(getvalue(x), 0.0, atol=TOL)
        @test isapprox(getvalue(y), 2.0, atol=TOL)
        @test isapprox(getvalue(z), 0.0, atol=TOL)

        # Test changing constraint bound (>= constraint)
        # max 1.1x + 1.0y + 100.0z
        # st     x +    y +      z <= 2
        #        x        +      z <= 1
        #        x +    y          >= 1
        #     0 <= x <= 3
        #     0 <= y <= 3
        #     0 <= z <= 1.5, Integer
        # x* = 1, y* = 0, z* = 1
        JuMP.setRHS(xz0ref, 2.0)
        xyg0ref = @constraint(m, x + y >= 0)
        @test solve(m) == :Optimal
        JuMP.setRHS(xyg0ref, 1)
        @test solve(m) == :Optimal
        @test isapprox(getvalue(x), 1.0, atol=TOL)
        @test isapprox(getvalue(y), 0.0, atol=TOL)
        @test isapprox(getvalue(z), 1.0, atol=TOL)
    end


    @testset "Adding a range constraint and modifying it" begin
        m = Model()
        @variable(m, x)
        rangeref = @constraint(m, -10 <= x <= 10)
        @test_throws ErrorException JuMP.setRHS(rangeref, 11)
    end


    @testset "Adding a 'decoupled' variable (#205) with $solver" for solver in lp_solvers
        m = Model(solver=solver)
        @variable(m, x >= 0)
        @objective(m, Min, x)
        solve(m)
        @variable(m, y >= 0)
        @constraint(m, x + y == 1)
        @objective(m, Min, 2x+y)
        solve(m)
        @test isapprox(getvalue(x), 0.0, atol=TOL)
        @test isapprox(getvalue(y), 1.0, atol=TOL)
    end


    @testset "JuMP.build with $solver" for solver in lp_solvers
        m = Model(solver=solver)
        @variable(m, x >= 0)
        @variable(m, y >= 0)
        @constraint(m, x + y == 1)
        @objective(m, Max, y)
        JuMP.build(m)
        @test internalmodel(m) !== nothing
        @test m.internalModelLoaded == true
        stat = solve(m)
        @test stat == :Optimal
        @test isapprox(getvalue(x), 0.0, atol=TOL)
        @test isapprox(getvalue(y), 1.0, atol=TOL)
        @test isapprox(getobjectivevalue(m), 1.0, atol=TOL)
        @test isapprox(getdual(x), -1.0, atol=TOL)
        @test isapprox(getdual(y),  0.0, atol=TOL)
    end


    @testset "JuMP.build (MIP) with $solver" for solver in ip_solvers
        m = Model(solver=solver)
        @variable(m, x >= 0, Int)
        @variable(m, y, Bin)
        @constraint(m, x + y == 1)
        @objective(m, Max, y)
        JuMP.build(m)
        @test internalmodel(m) !== nothing
        @test m.internalModelLoaded == true
        stat = solve(m)
        @test stat == :Optimal
        @test isapprox(getvalue(x), 0.0, atol=TOL)
        @test isapprox(getvalue(y), 1.0, atol=TOL)
        @test isapprox(getobjectivevalue(m), 1.0, atol=TOL)
    end


    @testset "Adding empty constraints with $solver" for solver in lp_solvers
        m = Model(solver=solver)
        @variable(m, 0 <= x <= 9)
        @objective(m, Max, x)
        @constraint(m, x <= 5)
        @constraint(m, 0 <= 1)
        solve(m)
        @constraint(m, 0 <= 1)
        @test solve(m) == :Optimal
        @test isapprox(getvalue(x), 5.0, atol=TOL)
    end


    @testset "Bound modification on binaries with $solver" for solver in ip_solvers
        # Test semantics for modifying bounds on binary variables:
        # Variables should be restricted to the intersection of
        # {0,1} and their bounds.
        mod = Model(solver=solver)
        @variable(mod, x, Bin)
        @objective(mod, Max, x)
        solve(mod)
        @test isapprox(getvalue(x), 1.0)
        setupperbound(x, 2.0)
        solve(mod)
        @test isapprox(getvalue(x), 1.0)
        setupperbound(x, 0.0)
        solve(mod)
        @test isapprox(getvalue(x), 0.0)
        # same thing, other direction
        mod = Model(solver=solver)
        @variable(mod, x, Bin)
        @objective(mod, Min, x)
        solve(mod)
        @test isapprox(getvalue(x), 0.0)
        setlowerbound(x, -1.0)
        solve(mod)
        @test isapprox(getvalue(x), 0.0)
        setlowerbound(x, 1.0)
        solve(mod)
        @test isapprox(getvalue(x), 1.0)
    end

    @testset "Switching from quadratic to linear objective with $solver" for solver in quad_solvers
        m = Model(solver=solver)
        @variable(m, x >=0)
        @variable(m, y >=0)

        @constraint(m, x + 2*y <= 4)
        @constraint(m, 2*x + y <= 4)

        @objective(m, Min, x^2) # Quadratic objective function
        @test solve(m) == :Optimal
        @test isapprox(getobjectivevalue(m), 0.0, atol=TOL)

        @objective(m, Max, 3*x + 4*y)  # Change to linear objective function
        @test solve(m) == :Optimal
        @test isapprox(getobjectivevalue(m), 9+1/3, atol=TOL)
    end

    @testset "Applicable regressions" begin

    function methods_test(solverobj, supp)
        mod = Model(solver=solverobj)
        @variable(mod, x >= 0)
        @constraint(mod, 2x == 2)
        solve(mod, suppress_warnings=true)
        internal_mod = internalmodel(mod)
        for (it,(meth, args)) in enumerate(mpb_methods)
            if supp[it]
                @test applicable(meth, internal_mod, args...)
                @test hasmethod(meth, map(typeof, tuple(internal_mod, args...)))
            end
        end
    end

    mpb_methods =
        [(MathProgBase.addquadconstr!, (Cint[1],Float64[1.0],Cint[1],Cint[1],Float64[1],'>',1.0)),
         (MathProgBase.setquadobjterms!, (Cint[1], Cint[1], Float64[1.0])),
         (MathProgBase.addconstr!,   ([1],[1.0],1.0,1.0)),
         (MathProgBase.addsos1!,     ([1],[1.0])),
         (MathProgBase.addsos2!,     ([1],[1.0])),
         (MathProgBase.addvar!,      ([1],[1.0],1.0,1.0,1.0)),
         (MathProgBase.setvarLB!,    ([1.0],)),
         (MathProgBase.setvarUB!,    ([1.0],)),
         (MathProgBase.setconstrLB!, ([1.0],)),
         (MathProgBase.setconstrUB!, ([1.0],)),
         (MathProgBase.setobj!,      ([1.0],)),
         (MathProgBase.setsense!,    (:Min,)),
         (MathProgBase.setvartype!,  ([:Cont],)),
         (MathProgBase.getinfeasibilityray, ()),
         (MathProgBase.getunboundedray, ()),
         (MathProgBase.getreducedcosts, ()),
         (MathProgBase.getconstrduals, ()),
         (MathProgBase.setwarmstart!, ([1.0]))]

    if grb
        supp = (true,true,true,true,true,true,true,true,true,true,true,true,true,true,true,true,true,true)
        @testset "with Gurobi" begin
            methods_test(Gurobi.GurobiSolver(OutputFlag=0), supp)
        end
    end
    if cpx
        supp = (true,true,true,true,true,true,true,true,true,true,true,true,true,true,true,true,true,true)
        @testset "with CPLEX" begin
            methods_test(CPLEX.CplexSolver(), supp)
        end
    end
    if xpr
        supp = (true,true,true,true,true,true,true,true,true,true,true,true,true,true,true,true,true,true)
        @testset "with Xpress" begin
            methods_test(Xpress.XpressSolver(), supp)
        end
    end
    if cbc
        supp = (false,false,true,false,false,true,true,true,true,true,true,true,true,true,true,true,true,false)
        @testset "with Clp" begin
            methods_test(Clp.ClpSolver(), supp)
        end
    end
    if cbc
        supp = (false,false,false,false,false,false,false,false,false,false,false,true,true,false,false,false,false,false)
        @testset "with Cbc" begin
            methods_test(Cbc.CbcSolver(), supp)
        end
    end
    if glp
        supp = (false,false,true,false,false,true,true,true,true,true,true,true,false,true,true,true,true,false)
        @testset "with GLPK" begin
            methods_test(GLPKMathProgInterface.GLPKSolverLP(), supp)
        supp = (false,false,true,false,false, true,true,true,true,true,true,true,true,false,false,false,false,false)
            methods_test(GLPKMathProgInterface.GLPKSolverMIP(), supp)
        end
    end
    if mos
        supp = (true,true,true,false,false,true,true,true,true,true,true,true,true,true,true,true,true,false)
        @testset "with Mosek" begin
            methods_test(Mosek.MosekSolver(), supp)
        end
    end

    end
end

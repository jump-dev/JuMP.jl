# Tests models with quadratic constraints
function qcqp_test(solvername, solverobj)
    println(string("  Running ", solvername))
    let
        println("    Test 1")
        modQ = Model(solver=solverobj)
        @defVar(modQ, -2 <= x <= 2 )
        @defVar(modQ, -2 <= y <= 2 )

        @setObjective(modQ, Min, x - y )
        addConstraint(modQ, x + x*x + x*y + y*y <= 1 )

        status = solve(modQ)
        @test status == :Optimal
        @test_approx_eq_eps modQ.objVal -1-4/sqrt(3) 1e-6
        @test_approx_eq_eps (getValue(x) + getValue(y)) -1/3 1e-3
    end

    let
        println("    Test 2")
        modQ = Model(solver=solverobj)

        @defVar(modQ, -2 <= x <= 2, Int )
        @defVar(modQ, -2 <= y <= 2, Int )

        @setObjective(modQ, Min, x - y )
        addConstraint(modQ, x + x*x + x*y + y*y <= 1 )

        status = solve(modQ)
        @test status == :Optimal
        @test_approx_eq_eps modQ.objVal -3 1e-6
        @test_approx_eq_eps (getValue(x) + getValue(y)) -1 1e-6
    end
end



if Pkg.installed("Gurobi") != nothing  
    using Gurobi
    qcqp_test("Gurobi", GurobiSolver(OutputFlag=0))
end
if Pkg.installed("CPLEXLink") != nothing
    using CPLEXLink
    qcqp_test("CPLEXLink", CplexSolver())
end
if Pkg.installed("Mosek") != nothing
    using Mosek
    qcqp_test("Mosek", MosekSolver())
end

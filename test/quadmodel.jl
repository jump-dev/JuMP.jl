# Test models with quadratic objectives
function qp_test(solvername, solverobj; prob_mod=true)
    println(string("  Running ", solvername))
    let
        println("    Test 1")
        modQ = Model(solver=solverobj)
        @defVar(modQ, 1.1*i <= x[i=1:3] <= 2.5*i, Int)

        setObjective(modQ, :Min, 10*x[1]*x[1] + 3*x[1]*x[2] + 5*x[2]*x[2] + 9*x[3]*x[3])

        @addConstraint(modQ, x[2] <= 1.7*x[3])
        @addConstraint(modQ, x[2] >= 0.5*x[1])

        status = solve(modQ)
        @test status == :Optimal
        @test_approx_eq_eps modQ.objVal 247.0 1e-5
        @test_approx_eq_eps getValue(x) [2.0, 3.0, 4.0] 1e-6
    end

    # test Maximization sense
    let
        println("    Test 2")
        modQ = Model(solver=solverobj)

        @defVar(modQ, 1.1*i <= x[i=1:3] <= 2.5*i, Int)

        setObjective(modQ, :Max, -10*x[1]*x[1] - 3*x[1]*x[2] - 5*x[2]*x[2] - 9*x[3]*x[3])

        @addConstraint(modQ, x[2] <= 1.7*x[3])
        @addConstraint(modQ, x[2] >= 0.5*x[1])

        status = solve(modQ)
        @test status == :Optimal
        @test_approx_eq_eps modQ.objVal -247.0 1e-5
        @test_approx_eq_eps getValue(x) [2.0, 3.0, 4.0] 1e-6
    end

    let
        println("    Test 3")
        modQ = Model(solver=solverobj)

        @defVar(modQ, 0.5 <= x <= 2 )
        @defVar(modQ, 0 <= y <= 30 )

        setObjective(modQ, :Min, (x+y)*(x+y) )
        @addConstraint(modQ, x + y >= 1 )

        status = solve(modQ)
        @test status == :Optimal
        @test_approx_eq_eps modQ.objVal 1.0 1e-6
        @test_approx_eq_eps (getValue(x) + getValue(y)) 1.0 1e-6
    end

    if prob_mod
        # test problem modification (quad constraint)
        let
            println("    Test 4")
            modQ = Model(solver=solverobj)
            @defVar(modQ, x >= 0)
            addConstraint(modQ, x*x <= 1)
            @setObjective(modQ, Max, x)
            stat = solve(modQ)
            @test_approx_eq_eps getObjectiveValue(modQ) 1.0 1e-5

            addConstraint(modQ, 2x*x <= 1)
            @test modQ.internalModelLoaded
            stat = solve(modQ)
            @test_approx_eq getObjectiveValue(modQ) sqrt(0.5)
        end

        # test problem modification (quad objective)
        let
            println("    Test 5")
            modQ = Model(solver=solverobj)
            @defVar(modQ,   0 <= x <= 1)
            @defVar(modQ, 1/2 <= y <= 1)
            setObjective(modQ, :Min, x*x - y)
            stat = solve(modQ)
            @test_approx_eq_eps getObjectiveValue(modQ) -1.0 1e-8

            setObjective(modQ, :Min, y*y - x)
            @test modQ.internalModelLoaded == true
            stat = solve(modQ)
            @test_approx_eq_eps getObjectiveValue(modQ) -0.75 1e-8
        end
    end
end



if Pkg.installed("Gurobi") != nothing  
    using Gurobi
    qp_test("Gurobi", GurobiSolver(OutputFlag=0))
end
if Pkg.installed("CPLEX") != nothing
    using CPLEX
    qp_test("CPLEX", CplexSolver())
end
if Pkg.installed("Mosek") != nothing
    using Mosek
    qp_test("Mosek", MosekSolver(); prob_mod=false) # weird issues with probmod to revisit
end

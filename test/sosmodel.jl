# Test models with SOS constraints
function sos_test(solvername, solverobj)
    println(string("  Running ", solvername))
    let
        println("    Test 1")
        modS = Model(solver=solverobj)

        @defVar(modS, x[1:3], Bin)
        @defVar(modS, y[1:5], Bin)
        @defVar(modS, z)
        @defVar(modS, w)

        setObjective(modS, :Max, z+w)

        a = [1,2,3]
        b = [5,4,7,2,1]

        @addConstraint(modS, z == sum{a[i]*x[i], i=1:3})
        @addConstraint(modS, w == sum{b[i]*y[i], i=1:5})

        addSOS1(modS, [a[i]x[i] for i in 1:3])
        addSOS2(modS, [b[i]y[i] for i in 1:5])

        status = solve(modS)
        @test status == :Optimal
        @test modS.objVal == 15.
        @test getValue(z) == 3.
        @test getValue(w) == 12. 
    end
end



if Pkg.installed("Gurobi") != nothing  
    using Gurobi
    sos_test("Gurobi", GurobiSolver(OutputFlag=0))
end
if Pkg.installed("CPLEX") != nothing
    using CPLEX
    sos_test("CPLEX", CplexSolver())
end
# if Pkg.installed("Cbc") != nothing
#     using Mosek
#     sos_test("Cbc", CbcSolver())
# end

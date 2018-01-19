@testset "Linear Programming" begin
    @testset "LP1" begin
        # simple 2 variable, 1 constraint problem
        # min -x
        # st   x + y <= 1   (x + y - 1 ∈ Nonpositives)
        #       x, y >= 0   (x, y ∈ Nonnegatives)

        m = Model()
        @variable(m, x >= 0.0)
        @variable(m, y >= 0.0)
        @objective(m, Min, -x)

        c = @constraint(m, x + y <= 1)

        JuMP.setlowerboundname(x, "xlb")
        JuMP.setlowerboundname(y, "ylb")
        JuMP.setname(c, "c")

        modelstring = """
        variables: x, y
        minobjective: -1.0*x
        xlb: x >= 0.0
        ylb: y >= 0.0
        c: x + y <= 1.0
        """

        instance = JuMP.JuMPInstance{Float64}()
        MOIU.loadfromstring!(instance, modelstring)
        MOIU.test_instances_equal(m.instance, instance, ["x","y"], ["c", "xlb", "ylb"])

        mocksolver = MOIU.MockSolverInstance(JuMP.JuMPInstance{Float64}())
        JuMP.attach(m, mocksolver)

        MOI.set!(mocksolver, MOI.TerminationStatus(), MOI.Success)
        MOI.set!(mocksolver, MOI.ObjectiveValue(), -1.0)
        MOI.set!(mocksolver, MOI.ResultCount(), 1)
        MOI.set!(mocksolver, MOI.PrimalStatus(), MOI.FeasiblePoint)
        MOI.set!(mocksolver, MOI.DualStatus(), MOI.FeasiblePoint)
        MOI.set!(mocksolver, MOI.VariablePrimal(), JuMP.solverinstanceindex(x), 1.0)
        MOI.set!(mocksolver, MOI.VariablePrimal(), JuMP.solverinstanceindex(y), 0.0)
        MOI.set!(mocksolver, MOI.ConstraintDual(), JuMP.solverinstanceindex(c), -1.0)

        JuMP.solve(m)

        @test JuMP.isattached(m)
        @test JuMP.hasvariableresult(m)

        @test JuMP.terminationstatus(m) == MOI.Success
        @test JuMP.primalstatus(m) == MOI.FeasiblePoint

        @test JuMP.resultvalue(x) == 1.0
        @test JuMP.resultvalue(y) == 0.0
        @test JuMP.resultvalue(x + y) == 1.0
        @test JuMP.objectivevalue(m) == -1.0

        @test JuMP.dualstatus(m) == MOI.FeasiblePoint
        @test JuMP.resultdual(c) ≈ -1 atol=1e-6
    end
end

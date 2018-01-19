@testset "Generation and solve with fake solver" begin
    @testset "LP" begin
        m = Model()
        @variable(m, x <= 2.0)
        @variable(m, y >= 0.0)
        @objective(m, Min, -x)

        c = @constraint(m, x + y <= 1)

        JuMP.setname(JuMP.UpperBoundRef(x), "xub")
        JuMP.setname(JuMP.LowerBoundRef(y), "ylb")
        JuMP.setname(c, "c")

        modelstring = """
        variables: x, y
        minobjective: -1.0*x
        xub: x <= 2.0
        ylb: y >= 0.0
        c: x + y <= 1.0
        """

        instance = JuMP.JuMPInstance{Float64}()
        MOIU.loadfromstring!(instance, modelstring)
        MOIU.test_instances_equal(m.instance, instance, ["x","y"], ["c", "xub", "ylb"])

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
        MOI.set!(mocksolver, MOI.ConstraintDual(), JuMP.solverinstanceindex(JuMP.UpperBoundRef(x)), 0.0)
        MOI.set!(mocksolver, MOI.ConstraintDual(), JuMP.solverinstanceindex(JuMP.LowerBoundRef(y)), 1.0)

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
        @test JuMP.resultdual(c) == -1
        @test JuMP.resultdual(JuMP.UpperBoundRef(x)) == 0.0
        @test JuMP.resultdual(JuMP.LowerBoundRef(y)) == 1.0
    end

    @testset "IP" begin
        m = Model()
        @variable(m, x == 1.0, Int)
        @variable(m, y, Bin)
        @objective(m, Max, x)

        JuMP.setname(JuMP.FixRef(x), "xfix")
        JuMP.setname(JuMP.IntegerRef(x), "xint")
        JuMP.setname(JuMP.BinaryRef(y), "ybin")

        modelstring = """
        variables: x, y
        maxobjective: 1.0*x
        xfix: x == 1.0
        xint: x in Integer()
        ybin: y in ZeroOne()
        """

        instance = JuMP.JuMPInstance{Float64}()
        MOIU.loadfromstring!(instance, modelstring)
        MOIU.test_instances_equal(m.instance, instance, ["x","y"], ["xfix", "xint", "ybin"])

        mocksolver = MOIU.MockSolverInstance(JuMP.JuMPInstance{Float64}())
        JuMP.attach(m, mocksolver)

        MOI.set!(mocksolver, MOI.TerminationStatus(), MOI.Success)
        MOI.set!(mocksolver, MOI.ObjectiveValue(), 1.0)
        MOI.set!(mocksolver, MOI.ResultCount(), 1)
        MOI.set!(mocksolver, MOI.PrimalStatus(), MOI.FeasiblePoint)
        MOI.set!(mocksolver, MOI.VariablePrimal(), JuMP.solverinstanceindex(x), 1.0)
        MOI.set!(mocksolver, MOI.VariablePrimal(), JuMP.solverinstanceindex(y), 0.0)

        JuMP.solve(m)

        @test JuMP.isattached(m)
        @test JuMP.hasvariableresult(m)

        @test JuMP.terminationstatus(m) == MOI.Success
        @test JuMP.primalstatus(m) == MOI.FeasiblePoint

        @test JuMP.resultvalue(x) == 1.0
        @test JuMP.resultvalue(y) == 0.0
        @test JuMP.objectivevalue(m) == 1.0

        @test !JuMP.hasresultdual(JuMP.FixRef(x))
        @test !JuMP.hasresultdual(JuMP.IntegerRef(x))
        @test !JuMP.hasresultdual(JuMP.BinaryRef(y))
    end

    @testset "SOC" begin
        m = Model()
        @variables m begin
            x
            y
            z
        end
        @objective(m, Max, x)
        # TODO: JuMP should generate a VectorOfVariables constraint here
        #varsoc = @constraint(m, [x,y,z] in MOI.SecondOrderCone(3))
        #JuMP.setname(varsoc, "varsoc")
        affsoc = @constraint(m, [x+y,z,1.0] in MOI.SecondOrderCone(3))
        JuMP.setname(affsoc, "affsoc")
        rotsoc = @constraint(m, [x+1,y,z] in MOI.RotatedSecondOrderCone(3))
        JuMP.setname(rotsoc, "rotsoc")

        modelstring = """
        variables: x, y, z
        maxobjective: 1.0*x
        #varsoc: [x,y,z] in SecondOrderCone(3)
        affsoc: [x+y,z,1.0] in SecondOrderCone(3)
        rotsoc: [x+1,y,z] in RotatedSecondOrderCone(3)
        """

        instance = JuMP.JuMPInstance{Float64}()
        MOIU.loadfromstring!(instance, modelstring)
        MOIU.test_instances_equal(m.instance, instance, ["x","y","z"], ["affsoc", "rotsoc"])

        mocksolver = MOIU.MockSolverInstance(JuMP.JuMPInstance{Float64}())
        JuMP.attach(m, mocksolver)

        MOI.set!(mocksolver, MOI.TerminationStatus(), MOI.Success)
        MOI.set!(mocksolver, MOI.ResultCount(), 1)
        MOI.set!(mocksolver, MOI.PrimalStatus(), MOI.FeasiblePoint)
        MOI.set!(mocksolver, MOI.DualStatus(), MOI.FeasiblePoint)
        MOI.set!(mocksolver, MOI.VariablePrimal(), JuMP.solverinstanceindex(x), 1.0)
        MOI.set!(mocksolver, MOI.VariablePrimal(), JuMP.solverinstanceindex(y), 0.0)
        MOI.set!(mocksolver, MOI.VariablePrimal(), JuMP.solverinstanceindex(z), 0.0)
        MOI.set!(mocksolver, MOI.ConstraintDual(), JuMP.solverinstanceindex(affsoc), [1.0,2.0,3.0])

        JuMP.solve(m)

        @test JuMP.isattached(m)
        @test JuMP.hasvariableresult(m)

        @test JuMP.terminationstatus(m) == MOI.Success
        @test JuMP.primalstatus(m) == MOI.FeasiblePoint

        @test JuMP.resultvalue(x) == 1.0
        @test JuMP.resultvalue(y) == 0.0
        @test JuMP.resultvalue(z) == 0.0

        @test JuMP.hasresultdual(affsoc)
        @test JuMP.resultdual(affsoc) == [1.0, 2.0, 3.0]
    end
end

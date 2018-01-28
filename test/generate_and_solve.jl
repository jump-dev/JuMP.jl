#  Copyright 2017, Iain Dunning, Joey Huchette, Miles Lubin, and contributors
#  This Source Code Form is subject to the terms of the Mozilla Public
#  License, v. 2.0. If a copy of the MPL was not distributed with this
#  file, You can obtain one at http://mozilla.org/MPL/2.0/.

# The tests here check JuMP's model generation and communication with solvers.
# Model generation is checked by comparing the internal instance with a serialized
# test instance (in MOIU's lightweight text format).
# Communication with solvers is tested by using a mock solver with solution data
# that we feed to it. Prior to using this testing approach, we would test JuMP
# by calling real solvers, which was flakey and slow.

# Note: No attempt is made to use correct solution data. We're only testing
# that the plumbing works. This could change if JuMP gains the ability to verify
# feasibility independently of a solver.

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
        MOIU.test_instances_equal(m.moibackend.instance, instance, ["x","y"], ["c", "xub", "ylb"])

        mocksolver = MOIU.MockSolverInstance(JuMP.JuMPInstance{Float64}())
        MOIU.resetsolver!(m, mocksolver)
        MOIU.attachsolver!(m)

        MOI.set!(mocksolver, MOI.TerminationStatus(), MOI.Success)
        MOI.set!(mocksolver, MOI.ObjectiveValue(), -1.0)
        MOI.set!(mocksolver, MOI.ResultCount(), 1)
        MOI.set!(mocksolver, MOI.PrimalStatus(), MOI.FeasiblePoint)
        MOI.set!(mocksolver, MOI.DualStatus(), MOI.FeasiblePoint)
        MOI.set!(mocksolver, MOI.VariablePrimal(), JuMP.solverindex(x), 1.0)
        MOI.set!(mocksolver, MOI.VariablePrimal(), JuMP.solverindex(y), 0.0)
        MOI.set!(mocksolver, MOI.ConstraintDual(), JuMP.solverindex(c), -1.0)
        MOI.set!(mocksolver, MOI.ConstraintDual(), JuMP.solverindex(JuMP.UpperBoundRef(x)), 0.0)
        MOI.set!(mocksolver, MOI.ConstraintDual(), JuMP.solverindex(JuMP.LowerBoundRef(y)), 1.0)

        JuMP.solve(m)

        #@test JuMP.isattached(m)
        @test JuMP.hasresultvalues(m)

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

    @testset "LP (Direct mode)" begin
        mocksolver = MOIU.MockSolverInstance(JuMP.JuMPInstance{Float64}())

        m = Model(mode = JuMP.Direct, solver = mocksolver)
        @variable(m, x <= 2.0)
        @variable(m, y >= 0.0)
        @objective(m, Min, -x)

        c = @constraint(m, x + y <= 1)
        MOI.set!(mocksolver, MOI.TerminationStatus(), MOI.Success)
        MOI.set!(mocksolver, MOI.ObjectiveValue(), -1.0)
        MOI.set!(mocksolver, MOI.ResultCount(), 1)
        MOI.set!(mocksolver, MOI.PrimalStatus(), MOI.FeasiblePoint)
        MOI.set!(mocksolver, MOI.DualStatus(), MOI.FeasiblePoint)
        MOI.set!(mocksolver, MOI.VariablePrimal(), JuMP.solverindex(x), 1.0)
        MOI.set!(mocksolver, MOI.VariablePrimal(), JuMP.solverindex(y), 0.0)
        MOI.set!(mocksolver, MOI.ConstraintDual(), JuMP.solverindex(c), -1.0)
        MOI.set!(mocksolver, MOI.ConstraintDual(), JuMP.solverindex(JuMP.UpperBoundRef(x)), 0.0)
        MOI.set!(mocksolver, MOI.ConstraintDual(), JuMP.solverindex(JuMP.LowerBoundRef(y)), 1.0)

        JuMP.solve(m)

        #@test JuMP.isattached(m)
        @test JuMP.hasresultvalues(m)

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

    # TODO: test Manual mode

    @testset "IP" begin
        mocksolver = MOIU.MockSolverInstance(JuMP.JuMPInstance{Float64}())
        # Tests the solver= keyword.
        m = Model(mode = JuMP.Automatic, solver = mocksolver)
        @variable(m, x == 1.0, Int)
        @variable(m, y, Bin)
        @objective(m, Max, x)

        JuMP.setname(JuMP.FixRef(x), "xfix")
        JuMP.setname(JuMP.IntegerRef(x), "xint")
        JuMP.setname(JuMP.BinaryRef(y), "ybin")

        modelstring = """
        variables: x, y
        maxobjective: x
        xfix: x == 1.0
        xint: x in Integer()
        ybin: y in ZeroOne()
        """

        instance = JuMP.JuMPInstance{Float64}()
        MOIU.loadfromstring!(instance, modelstring)
        MOIU.test_instances_equal(m.moibackend.instance, instance, ["x","y"], ["xfix", "xint", "ybin"])

        MOIU.attachsolver!(m)

        MOI.set!(mocksolver, MOI.TerminationStatus(), MOI.Success)
        MOI.set!(mocksolver, MOI.ObjectiveValue(), 1.0)
        MOI.set!(mocksolver, MOI.ResultCount(), 1)
        MOI.set!(mocksolver, MOI.PrimalStatus(), MOI.FeasiblePoint)
        MOI.set!(mocksolver, MOI.VariablePrimal(), JuMP.solverindex(x), 1.0)
        MOI.set!(mocksolver, MOI.VariablePrimal(), JuMP.solverindex(y), 0.0)

        JuMP.solve(m)

        #@test JuMP.isattached(m)
        @test JuMP.hasresultvalues(m)

        @test JuMP.terminationstatus(m) == MOI.Success
        @test JuMP.primalstatus(m) == MOI.FeasiblePoint

        @test JuMP.resultvalue(x) == 1.0
        @test JuMP.resultvalue(y) == 0.0
        @test JuMP.objectivevalue(m) == 1.0

        @test !JuMP.hasresultdual(m, typeof(JuMP.FixRef(x)))
        @test !JuMP.hasresultdual(m, typeof(JuMP.IntegerRef(x)))
        @test !JuMP.hasresultdual(m, typeof(JuMP.BinaryRef(y)))
    end

    @testset "QCQP" begin
        m = Model()
        @variable(m, x)
        @variable(m, y)
        @objective(m, Min, x^2)

        c1 = @constraint(m, 2x*y <= 1)
        c2 = @constraint(m, y^2 == x^2)
        c3 = @constraint(m, 2x + 3y*x >= 2)

        JuMP.setname(c1, "c1")
        JuMP.setname(c2, "c2")
        JuMP.setname(c3, "c3")

        modelstring = """
        variables: x, y
        minobjective: 1*x*x
        c1: 2*x*y <= 1.0
        c2: 1*y*y + -1*x*x == 0.0
        c3: 2x + 3*y*x >= 2.0
        """

        instance = JuMP.JuMPInstance{Float64}()
        MOIU.loadfromstring!(instance, modelstring)
        MOIU.test_instances_equal(m.moibackend.instance, instance, ["x","y"], ["c1", "c2", "c3"])

        mocksolver = MOIU.MockSolverInstance(JuMP.JuMPInstance{Float64}())
        MOIU.resetsolver!(m, mocksolver)
        MOIU.attachsolver!(m)

        MOI.set!(mocksolver, MOI.TerminationStatus(), MOI.Success)
        MOI.set!(mocksolver, MOI.ObjectiveValue(), -1.0)
        MOI.set!(mocksolver, MOI.ResultCount(), 1)
        MOI.set!(mocksolver, MOI.PrimalStatus(), MOI.FeasiblePoint)
        MOI.set!(mocksolver, MOI.DualStatus(), MOI.FeasiblePoint)
        MOI.set!(mocksolver, MOI.VariablePrimal(), JuMP.solverindex(x), 1.0)
        MOI.set!(mocksolver, MOI.VariablePrimal(), JuMP.solverindex(y), 0.0)
        MOI.set!(mocksolver, MOI.ConstraintDual(), JuMP.solverindex(c1), -1.0)
        MOI.set!(mocksolver, MOI.ConstraintDual(), JuMP.solverindex(c2), 2.0)
        MOI.set!(mocksolver, MOI.ConstraintDual(), JuMP.solverindex(c3), 3.0)

        JuMP.solve(m)

        #@test JuMP.isattached(m)
        @test JuMP.hasresultvalues(m)

        @test JuMP.terminationstatus(m) == MOI.Success
        @test JuMP.primalstatus(m) == MOI.FeasiblePoint

        @test JuMP.resultvalue(x) == 1.0
        @test JuMP.resultvalue(y) == 0.0
        @test JuMP.objectivevalue(m) == -1.0

        @test JuMP.dualstatus(m) == MOI.FeasiblePoint
        @test JuMP.resultdual(c1) == -1.0
        @test JuMP.resultdual(c2) == 2.0
        @test JuMP.resultdual(c3) == 3.0
    end

    @testset "SOC" begin
        m = Model()
        @variables m begin
            x
            y
            z
        end
        @objective(m, Max, 1.0*x)
        varsoc = @constraint(m, [x,y,z] in MOI.SecondOrderCone(3))
        JuMP.setname(varsoc, "varsoc")
        affsoc = @constraint(m, [x+y,z,1.0] in MOI.SecondOrderCone(3))
        JuMP.setname(affsoc, "affsoc")
        rotsoc = @constraint(m, [x+1,y,z] in MOI.RotatedSecondOrderCone(3))
        JuMP.setname(rotsoc, "rotsoc")

        modelstring = """
        variables: x, y, z
        maxobjective: 1.0*x
        varsoc: [x,y,z] in SecondOrderCone(3)
        affsoc: [x+y,z,1.0] in SecondOrderCone(3)
        rotsoc: [x+1,y,z] in RotatedSecondOrderCone(3)
        """

        instance = JuMP.JuMPInstance{Float64}()
        MOIU.loadfromstring!(instance, modelstring)
        MOIU.test_instances_equal(m.moibackend.instance, instance, ["x","y","z"], ["varsoc", "affsoc", "rotsoc"])

        mocksolver = MOIU.MockSolverInstance(JuMP.JuMPInstance{Float64}())
        MOIU.resetsolver!(m, mocksolver)
        MOIU.attachsolver!(m)

        MOI.set!(mocksolver, MOI.TerminationStatus(), MOI.Success)
        MOI.set!(mocksolver, MOI.ResultCount(), 1)
        MOI.set!(mocksolver, MOI.PrimalStatus(), MOI.FeasiblePoint)
        MOI.set!(mocksolver, MOI.DualStatus(), MOI.FeasiblePoint)
        MOI.set!(mocksolver, MOI.VariablePrimal(), JuMP.solverindex(x), 1.0)
        MOI.set!(mocksolver, MOI.VariablePrimal(), JuMP.solverindex(y), 0.0)
        MOI.set!(mocksolver, MOI.VariablePrimal(), JuMP.solverindex(z), 0.0)
        MOI.set!(mocksolver, MOI.ConstraintDual(), JuMP.solverindex(varsoc), [-1.0,-2.0,-3.0])
        MOI.set!(mocksolver, MOI.ConstraintDual(), JuMP.solverindex(affsoc), [1.0,2.0,3.0])

        JuMP.solve(m)

        #@test JuMP.isattached(m)
        @test JuMP.hasresultvalues(m)

        @test JuMP.terminationstatus(m) == MOI.Success
        @test JuMP.primalstatus(m) == MOI.FeasiblePoint

        @test JuMP.resultvalue(x) == 1.0
        @test JuMP.resultvalue(y) == 0.0
        @test JuMP.resultvalue(z) == 0.0

        @test JuMP.hasresultdual(m, typeof(varsoc))
        @test JuMP.resultdual(varsoc) == [-1.0, -2.0, -3.0]

        @test JuMP.hasresultdual(m, typeof(affsoc))
        @test JuMP.resultdual(affsoc) == [1.0, 2.0, 3.0]
    end

    @testset "SDP" begin
        m = Model()
        # TODO: PSD variable construction needs to be redone to make it possible
        # to access the corresponding MOI constraint.
        @variable(m, x[1:2,1:2], Symmetric) #, PSD)
        setname(x[1,1], "x11")
        setname(x[1,2], "x12")
        setname(x[2,2], "x22")
        @objective(m, Max, trace(x))
        conpsd = @SDconstraint(m, x ⪰ [1.0 0.0; 0.0 1.0])
        setname(conpsd, "conpsd")

        modelstring = """
        variables: x11, x12, x22
        maxobjective: 1.0*x11 + 1.0*x22
        conpsd: [x11 + -1.0,x12,x22 + -1.0] in PositiveSemidefiniteConeTriangle(2)
        """

        instance = JuMP.JuMPInstance{Float64}()
        MOIU.loadfromstring!(instance, modelstring)
        MOIU.test_instances_equal(m.moibackend.instance, instance, ["x11","x12","x22"], ["conpsd"])

        mocksolver = MOIU.MockSolverInstance(JuMP.JuMPInstance{Float64}())
        MOIU.resetsolver!(m, mocksolver)
        MOIU.attachsolver!(m)

        MOI.set!(mocksolver, MOI.TerminationStatus(), MOI.Success)
        MOI.set!(mocksolver, MOI.ResultCount(), 1)
        MOI.set!(mocksolver, MOI.PrimalStatus(), MOI.FeasiblePoint)
        MOI.set!(mocksolver, MOI.DualStatus(), MOI.FeasiblePoint)
        MOI.set!(mocksolver, MOI.VariablePrimal(), JuMP.solverindex(x[1,1]), 1.0)
        MOI.set!(mocksolver, MOI.VariablePrimal(), JuMP.solverindex(x[1,2]), 2.0)
        MOI.set!(mocksolver, MOI.VariablePrimal(), JuMP.solverindex(x[2,2]), 4.0)
        MOI.set!(mocksolver, MOI.ConstraintDual(), JuMP.solverindex(conpsd), [1.0,2.0,3.0])

        JuMP.solve(m)

        #@test JuMP.isattached(m)
        @test JuMP.hasresultvalues(m)

        @test JuMP.terminationstatus(m) == MOI.Success
        @test JuMP.primalstatus(m) == MOI.FeasiblePoint

        @test JuMP.resultvalue.(x) == [1.0 2.0; 2.0 4.0]
        @test JuMP.hasresultdual(m, typeof(conpsd))
        @test JuMP.resultdual(conpsd) == [1.0,2.0,3.0]

    end
end

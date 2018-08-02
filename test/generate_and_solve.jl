#  Copyright 2017, Iain Dunning, Joey Huchette, Miles Lubin, and contributors
#  This Source Code Form is subject to the terms of the Mozilla Public
#  License, v. 2.0. If a copy of the MPL was not distributed with this
#  file, You can obtain one at http://mozilla.org/MPL/2.0/.

# The tests here check JuMP's model generation and communication with solvers.
# Model generation is checked by comparing the internal model with a serialized
# test model (in MOIU's lightweight text format).
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

        model = JuMP.JuMPMOIModel{Float64}()
        MOIU.loadfromstring!(model, modelstring)
        MOIU.test_models_equal(JuMP.caching_optimizer(m).model_cache, model, ["x","y"], ["c", "xub", "ylb"])

        JuMP.optimize(m, with_optimizer(MOIU.MockOptimizer, JuMP.JuMPMOIModel{Float64}(), evalobjective=false))

        mockoptimizer = JuMP.caching_optimizer(m).optimizer
        MOI.set!(mockoptimizer, MOI.TerminationStatus(), MOI.Success)
        MOI.set!(mockoptimizer, MOI.ObjectiveValue(), -1.0)
        MOI.set!(mockoptimizer, MOI.ResultCount(), 1)
        MOI.set!(mockoptimizer, MOI.PrimalStatus(), MOI.FeasiblePoint)
        MOI.set!(mockoptimizer, MOI.DualStatus(), MOI.FeasiblePoint)
        MOI.set!(mockoptimizer, MOI.VariablePrimal(), JuMP.optimizerindex(x), 1.0)
        MOI.set!(mockoptimizer, MOI.VariablePrimal(), JuMP.optimizerindex(y), 0.0)
        MOI.set!(mockoptimizer, MOI.ConstraintDual(), JuMP.optimizerindex(c), -1.0)
        MOI.set!(mockoptimizer, MOI.ConstraintDual(), JuMP.optimizerindex(JuMP.UpperBoundRef(x)), 0.0)
        MOI.set!(mockoptimizer, MOI.ConstraintDual(), JuMP.optimizerindex(JuMP.LowerBoundRef(y)), 1.0)

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
        mockoptimizer = MOIU.MockOptimizer(JuMP.JuMPMOIModel{Float64}(), evalobjective=false)

        m = JuMP.direct_model(mockoptimizer)
        @variable(m, x <= 2.0)
        @variable(m, y >= 0.0)
        @objective(m, Min, -x)

        c = @constraint(m, x + y <= 1)
        MOI.set!(mockoptimizer, MOI.TerminationStatus(), MOI.Success)
        MOI.set!(mockoptimizer, MOI.ObjectiveValue(), -1.0)
        MOI.set!(mockoptimizer, MOI.ResultCount(), 1)
        MOI.set!(mockoptimizer, MOI.PrimalStatus(), MOI.FeasiblePoint)
        MOI.set!(mockoptimizer, MOI.DualStatus(), MOI.FeasiblePoint)
        MOI.set!(mockoptimizer, MOI.VariablePrimal(), JuMP.optimizerindex(x), 1.0)
        MOI.set!(mockoptimizer, MOI.VariablePrimal(), JuMP.optimizerindex(y), 0.0)
        MOI.set!(mockoptimizer, MOI.ConstraintDual(), JuMP.optimizerindex(c), -1.0)
        MOI.set!(mockoptimizer, MOI.ConstraintDual(), JuMP.optimizerindex(JuMP.UpperBoundRef(x)), 0.0)
        MOI.set!(mockoptimizer, MOI.ConstraintDual(), JuMP.optimizerindex(JuMP.LowerBoundRef(y)), 1.0)

        @test_throws ErrorException JuMP.optimize(m, with_optimizer(MOIU.MockOptimizer, JuMP.JuMPMOIModel{Float64}()))
        JuMP.optimize(m)

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
        # Tests the solver= keyword.
        m = Model(with_optimizer(MOIU.MockOptimizer, JuMP.JuMPMOIModel{Float64}(), evalobjective=false), caching_mode = MOIU.Automatic)
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

        model = JuMP.JuMPMOIModel{Float64}()
        MOIU.loadfromstring!(model, modelstring)
        MOIU.test_models_equal(JuMP.caching_optimizer(m).model_cache, model, ["x","y"], ["xfix", "xint", "ybin"])

        MOIU.attachoptimizer!(m)

        mockoptimizer = JuMP.caching_optimizer(m).optimizer
        MOI.set!(mockoptimizer, MOI.TerminationStatus(), MOI.Success)
        MOI.set!(mockoptimizer, MOI.ObjectiveValue(), 1.0)
        MOI.set!(mockoptimizer, MOI.ResultCount(), 1)
        MOI.set!(mockoptimizer, MOI.PrimalStatus(), MOI.FeasiblePoint)
        MOI.set!(mockoptimizer, MOI.VariablePrimal(), JuMP.optimizerindex(x), 1.0)
        MOI.set!(mockoptimizer, MOI.VariablePrimal(), JuMP.optimizerindex(y), 0.0)

        @test_throws ErrorException JuMP.optimize(m, with_optimizer(MOIU.MockOptimizer, JuMP.JuMPMOIModel{Float64}()))
        JuMP.optimize(m)

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

        @constraint(m, c1, 2x*y <= 1)
        @constraint(m, c2, y^2 == x^2)
        @constraint(m, c3, 2x + 3y*x >= 2)

        modelstring = """
        variables: x, y
        minobjective: 1*x*x
        c1: 2*x*y <= 1.0
        c2: 1*y*y + -1*x*x == 0.0
        c3: 2x + 3*y*x >= 2.0
        """

        model = JuMP.JuMPMOIModel{Float64}()
        MOIU.loadfromstring!(model, modelstring)
        MOIU.test_models_equal(JuMP.caching_optimizer(m).model_cache, model, ["x","y"], ["c1", "c2", "c3"])

        JuMP.optimize(m, with_optimizer(MOIU.MockOptimizer, JuMP.JuMPMOIModel{Float64}(), evalobjective=false))

        mockoptimizer = JuMP.caching_optimizer(m).optimizer
        MOI.set!(mockoptimizer, MOI.TerminationStatus(), MOI.Success)
        MOI.set!(mockoptimizer, MOI.ObjectiveValue(), -1.0)
        MOI.set!(mockoptimizer, MOI.ResultCount(), 1)
        MOI.set!(mockoptimizer, MOI.PrimalStatus(), MOI.FeasiblePoint)
        MOI.set!(mockoptimizer, MOI.DualStatus(), MOI.FeasiblePoint)
        MOI.set!(mockoptimizer, MOI.VariablePrimal(), JuMP.optimizerindex(x), 1.0)
        MOI.set!(mockoptimizer, MOI.VariablePrimal(), JuMP.optimizerindex(y), 0.0)
        MOI.set!(mockoptimizer, MOI.ConstraintDual(), JuMP.optimizerindex(c1), -1.0)
        MOI.set!(mockoptimizer, MOI.ConstraintDual(), JuMP.optimizerindex(c2), 2.0)
        MOI.set!(mockoptimizer, MOI.ConstraintDual(), JuMP.optimizerindex(c3), 3.0)

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
        @constraint(m, varsoc, [x,y,z] in SecondOrderCone())
        # Equivalent to `[x+y,z,1.0] in MOI.SOCone()`
        @constraint(m, affsoc, [x+y,z,1.0] in MOI.SecondOrderCone(3))
        @constraint(m, rotsoc, [x+1,y,z] in RotatedSecondOrderCone())

        modelstring = """
        variables: x, y, z
        maxobjective: 1.0*x
        varsoc: [x,y,z] in SecondOrderCone(3)
        affsoc: [x+y,z,1.0] in SecondOrderCone(3)
        rotsoc: [x+1,y,z] in RotatedSecondOrderCone(3)
        """

        model = JuMP.JuMPMOIModel{Float64}()
        MOIU.loadfromstring!(model, modelstring)
        MOIU.test_models_equal(JuMP.caching_optimizer(m).model_cache, model, ["x","y","z"], ["varsoc", "affsoc", "rotsoc"])

        mockoptimizer = MOIU.MockOptimizer(JuMP.JuMPMOIModel{Float64}(), evalobjective=false)
        MOIU.resetoptimizer!(m, mockoptimizer)
        MOIU.attachoptimizer!(m)

        MOI.set!(mockoptimizer, MOI.TerminationStatus(), MOI.Success)
        MOI.set!(mockoptimizer, MOI.ResultCount(), 1)
        MOI.set!(mockoptimizer, MOI.PrimalStatus(), MOI.FeasiblePoint)
        MOI.set!(mockoptimizer, MOI.DualStatus(), MOI.FeasiblePoint)
        MOI.set!(mockoptimizer, MOI.VariablePrimal(), JuMP.optimizerindex(x), 1.0)
        MOI.set!(mockoptimizer, MOI.VariablePrimal(), JuMP.optimizerindex(y), 0.0)
        MOI.set!(mockoptimizer, MOI.VariablePrimal(), JuMP.optimizerindex(z), 0.0)
        MOI.set!(mockoptimizer, MOI.ConstraintDual(), JuMP.optimizerindex(varsoc), [-1.0,-2.0,-3.0])
        MOI.set!(mockoptimizer, MOI.ConstraintDual(), JuMP.optimizerindex(affsoc), [1.0,2.0,3.0])

        JuMP.optimize(m)

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
        @variable(m, x[1:2,1:2], Symmetric)
        setname(x[1,1], "x11")
        setname(x[1,2], "x12")
        setname(x[2,2], "x22")
        @objective(m, Max, trace(x))
        varpsd = @constraint(m, x in PSDCone())
        setname(varpsd, "varpsd")
        conpsd = @SDconstraint(m, x ⪰ [1.0 0.0; 0.0 1.0])
        setname(conpsd, "conpsd")

        modelstring = """
        variables: x11, x12, x22
        maxobjective: 1.0*x11 + 1.0*x22
        varpsd: [x11,x12,x22] in PositiveSemidefiniteConeTriangle(2)
        conpsd: [x11 + -1.0,x12,x22 + -1.0] in PositiveSemidefiniteConeTriangle(2)
        """

        model = JuMP.JuMPMOIModel{Float64}()
        MOIU.loadfromstring!(model, modelstring)
        MOIU.test_models_equal(JuMP.caching_optimizer(m).model_cache, model, ["x11","x12","x22"], ["varpsd", "conpsd"])

        mockoptimizer = MOIU.MockOptimizer(JuMP.JuMPMOIModel{Float64}(), evalobjective=false)
        MOIU.resetoptimizer!(m, mockoptimizer)
        MOIU.attachoptimizer!(m)

        MOI.set!(mockoptimizer, MOI.TerminationStatus(), MOI.Success)
        MOI.set!(mockoptimizer, MOI.ResultCount(), 1)
        MOI.set!(mockoptimizer, MOI.PrimalStatus(), MOI.FeasiblePoint)
        MOI.set!(mockoptimizer, MOI.DualStatus(), MOI.FeasiblePoint)
        MOI.set!(mockoptimizer, MOI.VariablePrimal(), JuMP.optimizerindex(x[1,1]), 1.0)
        MOI.set!(mockoptimizer, MOI.VariablePrimal(), JuMP.optimizerindex(x[1,2]), 2.0)
        MOI.set!(mockoptimizer, MOI.VariablePrimal(), JuMP.optimizerindex(x[2,2]), 4.0)
        MOI.set!(mockoptimizer, MOI.ConstraintDual(), JuMP.optimizerindex(varpsd), [1.0,2.0,3.0])
        MOI.set!(mockoptimizer, MOI.ConstraintDual(), JuMP.optimizerindex(conpsd), [4.0,5.0,6.0])

        JuMP.optimize(m)

        #@test JuMP.isattached(m)
        @test JuMP.hasresultvalues(m)

        @test JuMP.terminationstatus(m) == MOI.Success
        @test JuMP.primalstatus(m) == MOI.FeasiblePoint

        @test JuMP.resultvalue.(x) == [1.0 2.0; 2.0 4.0]
        @test JuMP.hasresultdual(m, typeof(varpsd))
        @test JuMP.resultdual(varpsd) == [1.0,2.0,3.0]
        @test JuMP.hasresultdual(m, typeof(conpsd))
        @test JuMP.resultdual(conpsd) == [4.0,5.0,6.0]

    end
end

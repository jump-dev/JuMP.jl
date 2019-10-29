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

using LinearAlgebra, Test
using JuMP

@testset "Generation and solve with fake solver" begin
    @testset "LP" begin
        m = Model()
        @variable(m, x <= 2.0)
        @variable(m, y >= 0.0)
        @objective(m, Min, -x)

        c = @constraint(m, x + y <= 1)

        JuMP.set_name(JuMP.UpperBoundRef(x), "xub")
        JuMP.set_name(JuMP.LowerBoundRef(y), "ylb")
        JuMP.set_name(c, "c")

        modelstring = """
        variables: x, y
        minobjective: -1.0*x
        xub: x <= 2.0
        ylb: y >= 0.0
        c: x + y <= 1.0
        """

        model = MOIU.Model{Float64}()
        MOIU.loadfromstring!(model, modelstring)
        MOIU.test_models_equal(JuMP.backend(m).model_cache, model, ["x","y"], ["c", "xub", "ylb"])

        set_optimizer(m, () -> MOIU.MockOptimizer(
                                 MOIU.Model{Float64}(),
                                 eval_objective_value=false))
        JuMP.optimize!(m)

        mockoptimizer = JuMP.backend(m).optimizer.model
        MOI.set(mockoptimizer, MOI.TerminationStatus(), MOI.OPTIMAL)
        MOI.set(mockoptimizer, MOI.RawStatusString(), "solver specific string")
        MOI.set(mockoptimizer, MOI.ObjectiveValue(), -1.0)
        MOI.set(mockoptimizer, MOI.ResultCount(), 1)
        MOI.set(mockoptimizer, MOI.PrimalStatus(), MOI.FEASIBLE_POINT)
        MOI.set(mockoptimizer, MOI.DualStatus(), MOI.FEASIBLE_POINT)
        MOI.set(mockoptimizer, MOI.VariablePrimal(), JuMP.optimizer_index(x), 1.0)
        MOI.set(mockoptimizer, MOI.VariablePrimal(), JuMP.optimizer_index(y), 0.0)
        MOI.set(mockoptimizer, MOI.ConstraintDual(), JuMP.optimizer_index(c), -1.0)
        MOI.set(mockoptimizer, MOI.ConstraintDual(), JuMP.optimizer_index(JuMP.UpperBoundRef(x)), 0.0)
        MOI.set(mockoptimizer, MOI.ConstraintDual(), JuMP.optimizer_index(JuMP.LowerBoundRef(y)), 1.0)

        #@test JuMP.isattached(m)
        @test JuMP.has_values(m)

        @test MOI.OPTIMAL == @inferred JuMP.termination_status(m)
        @test "solver specific string" == JuMP.raw_status(m)
        @test MOI.FEASIBLE_POINT == @inferred JuMP.primal_status(m)

        @test  1.0 == @inferred JuMP.value(x)
        @test  0.0 == @inferred JuMP.value(y)
        @test  1.0 == @inferred JuMP.value(x + y)
        @test  1.0 == @inferred JuMP.value(c)
        @test -1.0 == @inferred JuMP.objective_value(m)
        @test -1.0 == @inferred JuMP.dual_objective_value(m)

        @test MOI.FEASIBLE_POINT == @inferred JuMP.dual_status(m)
        @test -1.0 == @inferred JuMP.dual(c)
        @test  0.0 == @inferred JuMP.dual(JuMP.UpperBoundRef(x))
        @test  1.0 == @inferred JuMP.dual(JuMP.LowerBoundRef(y))
    end

    @testset "LP (Direct mode)" begin
        mockoptimizer = MOIU.MockOptimizer(MOIU.Model{Float64}(),
                                           eval_objective_value=false)

        m = JuMP.direct_model(mockoptimizer)
        @variable(m, x <= 2.0)
        @variable(m, y >= 0.0)
        @objective(m, Min, -x)

        c = @constraint(m, x + y <= 1)
        MOI.set(mockoptimizer, MOI.TerminationStatus(), MOI.OPTIMAL)
        MOI.set(mockoptimizer, MOI.RawStatusString(), "solver specific string")
        MOI.set(mockoptimizer, MOI.ObjectiveValue(), -1.0)
        MOI.set(mockoptimizer, MOI.ResultCount(), 1)
        MOI.set(mockoptimizer, MOI.PrimalStatus(), MOI.FEASIBLE_POINT)
        MOI.set(mockoptimizer, MOI.DualStatus(), MOI.FEASIBLE_POINT)
        MOI.set(mockoptimizer, MOI.VariablePrimal(), JuMP.optimizer_index(x), 1.0)
        MOI.set(mockoptimizer, MOI.VariablePrimal(), JuMP.optimizer_index(y), 0.0)
        MOI.set(mockoptimizer, MOI.ConstraintDual(), JuMP.optimizer_index(c), -1.0)
        MOI.set(mockoptimizer, MOI.ConstraintDual(), JuMP.optimizer_index(JuMP.UpperBoundRef(x)), 0.0)
        MOI.set(mockoptimizer, MOI.ConstraintDual(), JuMP.optimizer_index(JuMP.LowerBoundRef(y)), 1.0)

        JuMP.optimize!(m)

        #@test JuMP.isattached(m)
        @test JuMP.has_values(m)

        @test MOI.OPTIMAL == @inferred JuMP.termination_status(m)
        @test "solver specific string" == JuMP.raw_status(m)
        @test MOI.FEASIBLE_POINT == @inferred JuMP.primal_status(m)

        @test  1.0 == @inferred JuMP.value(x)
        @test  0.0 == @inferred JuMP.value(y)
        @test  1.0 == @inferred JuMP.value(x + y)
        @test -1.0 == @inferred JuMP.objective_value(m)

        @test MOI.FEASIBLE_POINT == @inferred JuMP.dual_status(m)
        @test -1.0 == @inferred JuMP.dual(c)
        @test  0.0 == @inferred JuMP.dual(JuMP.UpperBoundRef(x))
        @test  1.0 == @inferred JuMP.dual(JuMP.LowerBoundRef(y))
    end

    # TODO: test Manual mode

    @testset "IP" begin
        # Tests the solver= keyword.
        m = Model(() -> MOIU.MockOptimizer(
                            MOIU.Model{Float64}(),
                            eval_objective_value=false),
                  caching_mode = MOIU.AUTOMATIC)
        @variable(m, x == 1.0, Int)
        @variable(m, y, Bin)
        @objective(m, Max, x)

        JuMP.set_name(JuMP.FixRef(x), "xfix")
        JuMP.set_name(JuMP.IntegerRef(x), "xint")
        JuMP.set_name(JuMP.BinaryRef(y), "ybin")

        modelstring = """
        variables: x, y
        maxobjective: x
        xfix: x == 1.0
        xint: x in Integer()
        ybin: y in ZeroOne()
        """

        model = MOIU.Model{Float64}()
        MOIU.loadfromstring!(model, modelstring)
        MOIU.test_models_equal(JuMP.backend(m).model_cache, model, ["x","y"], ["xfix", "xint", "ybin"])

        MOIU.attach_optimizer(m)

        mockoptimizer = JuMP.backend(m).optimizer.model
        MOI.set(mockoptimizer, MOI.TerminationStatus(), MOI.OPTIMAL)
        MOI.set(mockoptimizer, MOI.RawStatusString(), "solver specific string")
        MOI.set(mockoptimizer, MOI.ObjectiveValue(), 1.0)
        MOI.set(mockoptimizer, MOI.ResultCount(), 1)
        MOI.set(mockoptimizer, MOI.PrimalStatus(), MOI.FEASIBLE_POINT)
        MOI.set(mockoptimizer, MOI.VariablePrimal(), JuMP.optimizer_index(x), 1.0)
        MOI.set(mockoptimizer, MOI.VariablePrimal(), JuMP.optimizer_index(y), 0.0)
        MOI.set(mockoptimizer, MOI.DualStatus(), MOI.NO_SOLUTION)

        JuMP.optimize!(m)

        #@test JuMP.isattached(m)
        @test JuMP.has_values(m)

        @test MOI.OPTIMAL == @inferred JuMP.termination_status(m)
        @test "solver specific string" == JuMP.raw_status(m)
        @test MOI.FEASIBLE_POINT == @inferred JuMP.primal_status(m)

        @test 1.0 == @inferred JuMP.value(x)
        @test 0.0 == @inferred JuMP.value(y)
        @test 1.0 == @inferred JuMP.objective_value(m)

        @test !JuMP.has_duals(m)
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

        model = MOIU.Model{Float64}()
        MOIU.loadfromstring!(model, modelstring)
        MOIU.test_models_equal(JuMP.backend(m).model_cache, model, ["x","y"], ["c1", "c2", "c3"])

        set_optimizer(m, () -> MOIU.MockOptimizer(
                                 MOIU.Model{Float64}(),
                                 eval_objective_value=false))
        JuMP.optimize!(m)

        mockoptimizer = JuMP.backend(m).optimizer.model
        MOI.set(mockoptimizer, MOI.TerminationStatus(), MOI.OPTIMAL)
        MOI.set(mockoptimizer, MOI.RawStatusString(), "solver specific string")
        MOI.set(mockoptimizer, MOI.ObjectiveValue(), -1.0)
        MOI.set(mockoptimizer, MOI.ResultCount(), 1)
        MOI.set(mockoptimizer, MOI.PrimalStatus(), MOI.FEASIBLE_POINT)
        MOI.set(mockoptimizer, MOI.DualStatus(), MOI.FEASIBLE_POINT)
        MOI.set(mockoptimizer, MOI.VariablePrimal(), JuMP.optimizer_index(x), 1.0)
        MOI.set(mockoptimizer, MOI.VariablePrimal(), JuMP.optimizer_index(y), 0.0)
        MOI.set(mockoptimizer, MOI.ConstraintDual(), JuMP.optimizer_index(c1), -1.0)
        MOI.set(mockoptimizer, MOI.ConstraintDual(), JuMP.optimizer_index(c2), 2.0)
        MOI.set(mockoptimizer, MOI.ConstraintDual(), JuMP.optimizer_index(c3), 3.0)

        #@test JuMP.isattached(m)
        @test JuMP.has_values(m)

        @test MOI.OPTIMAL == @inferred JuMP.termination_status(m)
        @test "solver specific string" == JuMP.raw_status(m)
        @test MOI.FEASIBLE_POINT == @inferred JuMP.primal_status(m)

        @test 1.0 == @inferred JuMP.value(x)
        @test 0.0 == @inferred JuMP.value(y)
        @test -1.0 == @inferred JuMP.objective_value(m)
        @test 5.0 == @inferred JuMP.dual_objective_value(m)

        @test MOI.FEASIBLE_POINT == @inferred JuMP.dual_status(m)
        @test -1.0 == @inferred JuMP.dual(c1)
        @test 2.0 == @inferred JuMP.dual(c2)
        @test 3.0 == @inferred JuMP.dual(c3)

        @test 2.0 == @inferred JuMP.value(2 * x + 3 * y * x)
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
        # Equivalent to `[x+y,z,1.0] in SecondOrderCone()`
        @constraint(m, affsoc, [x+y,z,1.0] in MOI.SecondOrderCone(3))
        @constraint(m, rotsoc, [x+1,y,z] in RotatedSecondOrderCone())

        modelstring = """
        variables: x, y, z
        maxobjective: 1.0*x
        varsoc: [x,y,z] in SecondOrderCone(3)
        affsoc: [x+y,z,1.0] in SecondOrderCone(3)
        rotsoc: [x+1,y,z] in RotatedSecondOrderCone(3)
        """

        model = MOIU.Model{Float64}()
        MOIU.loadfromstring!(model, modelstring)
        MOIU.test_models_equal(JuMP.backend(m).model_cache, model, ["x","y","z"], ["varsoc", "affsoc", "rotsoc"])

        mockoptimizer = MOIU.MockOptimizer(MOIU.Model{Float64}(),
                                           eval_objective_value=false,
                                           eval_variable_constraint_dual=false)
        MOIU.reset_optimizer(m, mockoptimizer)
        MOIU.attach_optimizer(m)

        MOI.set(mockoptimizer, MOI.TerminationStatus(), MOI.OPTIMAL)
        MOI.set(mockoptimizer, MOI.RawStatusString(), "solver specific string")
        MOI.set(mockoptimizer, MOI.ResultCount(), 1)
        MOI.set(mockoptimizer, MOI.PrimalStatus(), MOI.FEASIBLE_POINT)
        MOI.set(mockoptimizer, MOI.DualStatus(), MOI.FEASIBLE_POINT)
        MOI.set(mockoptimizer, MOI.VariablePrimal(), JuMP.optimizer_index(x), 1.0)
        MOI.set(mockoptimizer, MOI.VariablePrimal(), JuMP.optimizer_index(y), 0.0)
        MOI.set(mockoptimizer, MOI.VariablePrimal(), JuMP.optimizer_index(z), 0.0)
        MOI.set(mockoptimizer, MOI.ConstraintDual(), JuMP.optimizer_index(varsoc), [-1.0,-2.0,-3.0])
        MOI.set(mockoptimizer, MOI.ConstraintDual(), JuMP.optimizer_index(affsoc), [1.0,2.0,3.0])

        JuMP.optimize!(m)

        #@test JuMP.isattached(m)
        @test JuMP.has_values(m)

        @test MOI.OPTIMAL == @inferred JuMP.termination_status(m)
        @test "solver specific string" == JuMP.raw_status(m)
        @test MOI.FEASIBLE_POINT == @inferred JuMP.primal_status(m)

        @test 1.0 == @inferred JuMP.value(x)
        @test 0.0 == @inferred JuMP.value(y)
        @test 0.0 == @inferred JuMP.value(z)

        @test JuMP.has_duals(m)
        @test [-1.0, -2.0, -3.0] == @inferred JuMP.dual(varsoc)
        @test [ 1.0,  2.0,  3.0] == @inferred JuMP.dual(affsoc)
    end

    @testset "SDP" begin
        m = Model()
        @variable(m, x[1:2,1:2], Symmetric)
        set_name(x[1,1], "x11")
        set_name(x[1,2], "x12")
        set_name(x[2,2], "x22")
        @objective(m, Max, tr(x))
        var_psd = @constraint(m, x in PSDCone())
        set_name(var_psd, "var_psd")
        sym_psd = @constraint(m, Symmetric(x - [1.0 0.0; 0.0 1.0]) in PSDCone())
        set_name(sym_psd, "sym_psd")
        con_psd = @SDconstraint(m, x ⪰ [1.0 0.0; 0.0 1.0])
        set_name(con_psd, "con_psd")

        modelstring = """
        variables: x11, x12, x22
        maxobjective: 1.0*x11 + 1.0*x22
        var_psd: [x11,x12,x22] in PositiveSemidefiniteConeTriangle(2)
        sym_psd: [x11 + -1.0,x12,x22 + -1.0] in PositiveSemidefiniteConeTriangle(2)
        con_psd: [x11 + -1.0,x12,x12,x22 + -1.0] in PositiveSemidefiniteConeSquare(2)
        """

        model = MOIU.Model{Float64}()
        MOIU.loadfromstring!(model, modelstring)
        MOIU.test_models_equal(JuMP.backend(m).model_cache, model,
                               ["x11","x12","x22"],
                               ["var_psd", "sym_psd", "con_psd"])

        mockoptimizer = MOIU.MockOptimizer(MOIU.Model{Float64}(),
                                           eval_objective_value=false,
                                           eval_variable_constraint_dual=false)
        MOIU.reset_optimizer(m, mockoptimizer)
        MOIU.attach_optimizer(m)

        MOI.set(mockoptimizer, MOI.TerminationStatus(), MOI.OPTIMAL)
        MOI.set(mockoptimizer, MOI.RawStatusString(), "solver specific string")
        MOI.set(mockoptimizer, MOI.ResultCount(), 1)
        MOI.set(mockoptimizer, MOI.PrimalStatus(), MOI.FEASIBLE_POINT)
        MOI.set(mockoptimizer, MOI.DualStatus(), MOI.FEASIBLE_POINT)
        MOI.set(mockoptimizer, MOI.VariablePrimal(), JuMP.optimizer_index(x[1,1]), 1.0)
        MOI.set(mockoptimizer, MOI.VariablePrimal(), JuMP.optimizer_index(x[1,2]), 2.0)
        MOI.set(mockoptimizer, MOI.VariablePrimal(), JuMP.optimizer_index(x[2,2]), 4.0)
        MOI.set(mockoptimizer, MOI.ConstraintDual(),
                JuMP.optimizer_index(var_psd), [1.0, 2.0, 3.0])
        MOI.set(mockoptimizer, MOI.ConstraintDual(),
                JuMP.optimizer_index(sym_psd), [4.0, 5.0, 6.0])
        MOI.set(mockoptimizer, MOI.ConstraintDual(),
                JuMP.optimizer_index(con_psd), [7.0, 8.0, 9.0, 10.0])

        JuMP.optimize!(m)

        @test MOI.OPTIMAL == @inferred JuMP.termination_status(m)
        @test "solver specific string" == JuMP.raw_status(m)
        @test MOI.FEASIBLE_POINT == @inferred JuMP.primal_status(m)

        @test JuMP.has_values(m)
        @test [1.0 2.0; 2.0 4.0] == JuMP.value.(x)
        @test JuMP.value(var_psd) isa Symmetric
        @test [1.0 2.0; 2.0 4.0] == @inferred JuMP.value(var_psd)
        @test JuMP.value(sym_psd) isa Symmetric
        @test [0.0 2.0; 2.0 3.0] == @inferred JuMP.value(sym_psd)
        @test JuMP.value(con_psd) isa Matrix
        @test [0.0 2.0; 2.0 3.0] == @inferred JuMP.value(con_psd)

        @test JuMP.has_duals(m)
        @test JuMP.dual(var_psd) isa Symmetric
        @test [1.0 2.0; 2.0 3.0] == @inferred JuMP.dual(var_psd)
        @test JuMP.dual(sym_psd) isa Symmetric
        @test [4.0 5.0; 5.0 6.0] == @inferred JuMP.dual(sym_psd)
        @test JuMP.dual(con_psd) isa Matrix
        @test [7.0 9.0; 8.0 10.0] == @inferred JuMP.dual(con_psd)

    end

    @testset "Solver doesn't support nonlinear constraints" begin
        model = Model(() -> MOIU.MockOptimizer(MOIU.Model{Float64}()))
        @variable(model, x)
        @NLobjective(model, Min, sin(x))
        err = ErrorException("The solver does not support nonlinear problems " *
                             "(i.e., NLobjective and NLconstraint).")
        @test_throws err JuMP.optimize!(model)
    end
end

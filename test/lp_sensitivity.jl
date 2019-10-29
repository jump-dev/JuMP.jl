#  Copyright 2017, Iain Dunning, Joey Huchette, Miles Lubin, and contributors
#  This Source Code Form is subject to the terms of the Mozilla Public
#  License, v. 2.0. If a copy of the MPL was not distributed with this
#  file, You can obtain one at http://mozilla.org/MPL/2.0/.
#############################################################################
# JuMP
# An algebraic modeling language for Julia
# See http://github.com/JuliaOpt/JuMP.jl
#############################################################################
function test_lp_rhs_perturbation_range(model_string, primal_solution, basis_status, feasibility_ranges)
    model = JuMP.Model()
    MOIU.loadfromstring!(JuMP.backend(model), model_string)
    set_optimizer(model, () -> MOIU.MockOptimizer(
                                 MOIU.Model{Float64}(),
                                 eval_variable_constraint_dual=false))
    JuMP.optimize!(model)
    mock_optimizer = JuMP.backend(model).optimizer.model
    MOI.set(mock_optimizer, MOI.TerminationStatus(), MOI.OPTIMAL)

    JuMP.optimize!(model)
    for (variable_name, value) in primal_solution
        vi = MOI.get(JuMP.backend(model), MOI.VariableIndex, variable_name)
        variable_ref = JuMP.VariableRef(model, vi)
        MOI.set(mock_optimizer,  MOI.VariablePrimal(),
                JuMP.optimizer_index(variable_ref), value)
    end

    for (constraint_name, value) in basis_status
        ci = MOI.get(JuMP.backend(model), MOI.ConstraintIndex, constraint_name)
        constraint_ref = JuMP.ConstraintRef(model, ci, JuMP.ScalarShape())
        MOI.set(mock_optimizer,  MOI.ConstraintBasisStatus(),
                JuMP.optimizer_index(constraint_ref), value)
    end

    @testset "rhs range of $constraint_name" for constraint_name in keys(feasibility_ranges)
        ci = MOI.get(JuMP.backend(model), MOI.ConstraintIndex, constraint_name)
        constraint_ref = JuMP.ConstraintRef(model, ci, JuMP.ScalarShape())

        true_range = feasibility_ranges[constraint_name]
        range = lp_rhs_perturbation_range(constraint_ref)

        @test true_range[1] ≈ range[1]
        @test true_range[2] ≈ range[2]
    end
end

@testset "lp_rhs_perturbation_range" begin
    test_lp_rhs_perturbation_range("""
        variables: x, y, z, w
        xint: x in Interval(-1.0, 1.0)
        ylb: y >= 0.0
        zval: z == 1.0
        c1: x + y + z + w == 1.0
        c2: x + y <= 2.0
        """,
        # Optimal primals
        Dict("x" => 1.0, "y" => 1.0, "z" => 1.0, "w" => -2.0),
        # Basis status
        Dict("xint" => MOI.NONBASIC_AT_UPPER, "ylb" => MOI.BASIC, "zval" => MOI.NONBASIC, "c1" => MOI.NONBASIC, "c2" => MOI.NONBASIC),
        # Expected perturbation ranges
        Dict("ylb" => (-Inf, 1.0), "zval" => (-Inf, Inf), "c1" => (-Inf, Inf), "c2" => (-1.0, Inf)))

    test_lp_rhs_perturbation_range("""
        variables: x, y
        xub: x <= 1.0
        ylb: y >= 0.0
        c1: x + y <= 2.0
        c2: 2.0*x + y <= 3.5
        """,
        # Optimal primals
        Dict("x" => 1.0, "y" => 1.0),
        # Basis status
        Dict("xub" => MOI.NONBASIC, "ylb" => MOI.BASIC, "c1" => MOI.NONBASIC, "c2" => MOI.BASIC),
        # Expected perturbation ranges
        Dict("xub" => (-Inf, 0.5), "ylb" => (-Inf, 1.0), "c1" => (-1.0, 0.5), "c2" => (-0.5, Inf)))

    test_lp_rhs_perturbation_range("""
        variables: x, y
        xlb: x >= -1.0
        xub: x <= 1.0
        ylb: y >= 0.0
        c1: x + y <= 2.0
        c2: 2.0*x + y <= 3.0
        """,
        # Optimal primals
        Dict("x" => 1.0, "y" => 1.0),
        # Basis status
        Dict("xlb" => MOI.BASIC, "xub" => MOI.NONBASIC, "ylb" => MOI.BASIC, "c1" => MOI.BASIC, "c2" => MOI.NONBASIC),
        # Expected perturbation ranges
        Dict("xlb" => (-Inf, 2.0), "xub" => (0.0, 0.5), "ylb" => (-Inf, 1.0), "c1" => (0.0, Inf), "c2" => (-1.0, 0.0)))

    test_lp_rhs_perturbation_range("""
        variables: x, y
        xlb: x >= 0.0
        ylb: y >= 0.0
        c1: x + y >= 2.0
        c2: x + -1.0*y in Interval(-1.0, 1.0)
        """,
        # Optimal primals
        Dict("x" => 1.5, "y" => 0.5),
        # Basis status
        Dict("xlb" => MOI.BASIC, "ylb" => MOI.BASIC, "c1" => MOI.NONBASIC, "c2" => MOI.NONBASIC_AT_UPPER),
        # Expected perturbation ranges
        Dict("xlb" => (-Inf, 1.5), "ylb" => (-Inf, 0.5), "c1" => (-1.0, Inf)))
end

function test_lp_objective_perturbation_range(model_string, dual_solution, basis_status, optimality_ranges)
    model = JuMP.Model()
    MOIU.loadfromstring!(JuMP.backend(model), model_string)
    set_optimizer(model, () -> MOIU.MockOptimizer(
                                 MOIU.Model{Float64}(),
                                 eval_variable_constraint_dual=true))
    JuMP.optimize!(model)
    mock_optimizer = JuMP.backend(model).optimizer.model
    MOI.set(mock_optimizer, MOI.TerminationStatus(), MOI.OPTIMAL)
    MOI.set(mock_optimizer, MOI.DualStatus(), MOI.FEASIBLE_POINT)

    JuMP.optimize!(model)
    for (constraint_name, value) in dual_solution
        ci = MOI.get(JuMP.backend(model), MOI.ConstraintIndex, constraint_name)
        constraint_ref = JuMP.ConstraintRef(model, ci, JuMP.ScalarShape())
        MOI.set(mock_optimizer,  MOI.ConstraintDual(),
                JuMP.optimizer_index(constraint_ref), value)
    end

    for (constraint_name, value) in basis_status
        ci = MOI.get(JuMP.backend(model), MOI.ConstraintIndex, constraint_name)
        constraint_ref = JuMP.ConstraintRef(model, ci, JuMP.ScalarShape())
        MOI.set(mock_optimizer,  MOI.ConstraintBasisStatus(),
                JuMP.optimizer_index(constraint_ref), value)
    end

    @testset "optimality range of $variable_name" for variable_name in keys(optimality_ranges)
        vi = MOI.get(JuMP.backend(model), MOI.VariableIndex, variable_name)
        variable_ref = JuMP.VariableRef(model, vi)

        true_range = optimality_ranges[variable_name]
        range = lp_objective_perturbation_range(variable_ref)

        @test true_range[1] ≈ range[1]
        @test true_range[2] ≈ range[2]
    end
end

@testset "lp_objective_perturbation_range" begin
    test_lp_objective_perturbation_range("""
        variables: x, y, z, w
        minobjective: -1.0*x + -1.0*y
        xint: x in Interval(-1.0, 1.0)
        ylb: y >= 0.0
        zval: z == 1.0
        c1: x + y + z + w == 1.0
        c2: x + y <= 2.0
        """,
        # Optimal duals
        Dict("c1" => 0.0, "c2" => -1.0),
        # Basis status
        Dict("xint" => MOI.NONBASIC_AT_UPPER, "ylb" => MOI.BASIC, "zval" => MOI.NONBASIC, "c1" => MOI.NONBASIC, "c2" => MOI.NONBASIC),
        # Expected perturbation ranges
        Dict("x" => (-Inf, 0.0), "y" => (0.0, 1.0), "z" => (-Inf, Inf), "w" => (-1.0, Inf)))

    test_lp_objective_perturbation_range("""
        variables: x, y, z, w
        maxobjective: 1.0*x + 1.0*y
        xint: x in Interval(-1.0, 1.0)
        ylb: y >= 0.0
        zval: z == 1.0
        c1: x + y + z + w == 1.0
        c2: x + y <= 2.0
        """,
        # Optimal duals
        Dict("c1" => 0.0, "c2" => -1.0),
        # Basis status
        Dict("xint" => MOI.NONBASIC_AT_UPPER, "ylb" => MOI.BASIC, "zval" => MOI.NONBASIC, "c1" => MOI.NONBASIC, "c2" => MOI.NONBASIC),
        # Expected perturbation ranges
        Dict("x" => (0.0, Inf), "y" => (-1.0, 0.0), "z" => (-Inf, Inf), "w" => (-Inf, 1.0)))

    test_lp_objective_perturbation_range("""
        variables: x, y
        maxobjective: 1.0*x + 1.0*y
        xlb: x >= 0.0
        ylb: y >= 0.0
        c1: 1.0*x + 2.0*y in Interval(-1.0, 2.0)
        c2: x + y >= 0.5
        c3: 2.0*x + 1.0*y <= 2.0
        """,
        # Optimal duals
        Dict("c1" => -1/3, "c2" => 0.0, "c3" => -1/3),
        # Basis status
        Dict("xlb" => MOI.BASIC, "ylb" => MOI.BASIC, "c1" => MOI.NONBASIC_AT_UPPER, "c2" => MOI.BASIC, "c3" => MOI.NONBASIC),
        # Expected perturbation ranges
        Dict("x" => (-0.5, 1.0), "y" => (-0.5, 1.0)))

    test_lp_objective_perturbation_range("""
        variables: x, y
        maxobjective: 1.0*x + 1.0*y
        xlb: x >= 0.0
        ylb: y >= 0.0
        c1: 1.0*x + 2.0*y in Interval(-1.0, 2.0)
        c2: 1.0*x + 1.0*y >= 0.5
        c3: 2.0*x + 1.0*y <= 2.0
        """,
        # Optimal duals
        Dict("c1" => -1/3, "c2" => 0.0, "c3" => -1/3),
        # Basis status
        Dict("xlb" => MOI.BASIC, "ylb" => MOI.BASIC, "c1" => MOI.NONBASIC_AT_UPPER, "c2" => MOI.BASIC, "c3" => MOI.NONBASIC),
        # Expected perturbation ranges
        Dict("x" => (-0.5, 1.0), "y" => (-0.5, 1.0)))

    test_lp_objective_perturbation_range("""
        variables: x, y
        minobjective: 1.0*x + 1.0*y
        xlb: x >= 0.0
        ylb: y >= 0.0
        c1: 1.0*x + 2.0*y in Interval(-1.0, 2.0)
        c2: 1.0*x + 1.0*y >= 0.5
        c3: 2.0*x + 1.0*y <= 2.0
        """,
        # Optimal duals
        Dict("c1" => 0.0, "c2" => 1.0, "c3" => 0.0),
        # Basis status
        Dict("xlb" => MOI.BASIC, "ylb" => MOI.NONBASIC, "c1" => MOI.BASIC, "c2" => MOI.NONBASIC, "c3" => MOI.BASIC),
        # Expected perturbation ranges
        Dict("x" => (-1.0, 0.0), "y" => (0.0, Inf)))
end

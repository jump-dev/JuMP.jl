#  Copyright 2017, Iain Dunning, Joey Huchette, Miles Lubin, and contributors
#  This Source Code Form is subject to the terms of the Mozilla Public
#  License, v. 2.0. If a copy of the MPL was not distributed with this
#  file, You can obtain one at http://mozilla.org/MPL/2.0/.

module TestLPSensitivity

using Test
using JuMP
import LinearAlgebra

struct TestSensitivitySolution
    primal::Float64
    dual::Float64
    basis::MOI.BasisStatusCode
    range::Tuple{Float64,Float64}
end

function _test_sensitivity(model_string, solution)
    m = MOIU.MockOptimizer(
        MOIU.Model{Float64}();
        eval_variable_constraint_dual = false,
    )
    model = direct_model(m)
    MOI.Utilities.loadfromstring!(m, model_string)
    optimize!(model)
    MOI.set(m, MOI.TerminationStatus(), MOI.OPTIMAL)
    MOI.set(m, MOI.PrimalStatus(), MOI.FEASIBLE_POINT)
    MOI.set(m, MOI.DualStatus(), MOI.FEASIBLE_POINT)
    obj_map = Dict{Any,Any}()
    for (key, val) in solution
        if key isa String
            var = variable_by_name(model, key)
            if var !== nothing
                obj_map[key] = var
                MOI.set(m, MOI.VariablePrimal(), index(var), val.primal)
                MOI.set(m, MOI.VariableBasisStatus(), index(var), val.basis)
            else
                c = constraint_by_name(model, key)
                @assert c !== nothing
                obj_map[key] = c
                MOI.set(m, MOI.ConstraintDual(), index(c), val.dual)
                MOI.set(m, MOI.ConstraintBasisStatus(), index(c), val.basis)
            end
        else
            var = variable_by_name(model, key[1])
            c = if key[2] <: MOI.GreaterThan
                LowerBoundRef(var)
            elseif key[2] <: MOI.LessThan
                UpperBoundRef(var)
            else
                FixRef(var)
            end
            @assert c !== nothing
            obj_map[key] = c
            MOI.set(m, MOI.ConstraintDual(), index(c), val.dual)
            MOI.set(m, MOI.ConstraintBasisStatus(), index(c), val.basis)
        end
    end
    sens = lp_sensitivity_report(model)
    for (s_key, val) in solution
        key = obj_map[s_key]
        @test sens[key][1] ≈ val.range[1]
        @test sens[key][2] ≈ val.range[2]
    end
    return
end

function test_Error_handling()
    m = MOIU.MockOptimizer(
        MOIU.Model{Float64}();
        eval_variable_constraint_dual = false,
    )
    model = direct_model(m)
    optimize!(model)
    MOI.set(m, MOI.TerminationStatus(), MOI.OPTIMAL)
    MOI.set(m, MOI.PrimalStatus(), MOI.NO_SOLUTION)
    MOI.set(m, MOI.DualStatus(), MOI.NO_SOLUTION)
    err = ErrorException(
        "Unable to compute LP sensitivity: no primal solution available.",
    )
    @test_throws err lp_sensitivity_report(model)
    MOI.set(m, MOI.PrimalStatus(), MOI.FEASIBLE_POINT)
    err = ErrorException(
        "Unable to compute LP sensitivity: no dual solution available.",
    )
    @test_throws err lp_sensitivity_report(model)
    MOI.Utilities.loadfromstring!(
        m,
        """
variables: x
maxobjective: 1.0 * x
c: x in Interval(0.0, 1.0)
""",
    )
    err = ErrorException(
        "Unable to compute LP sensitivity because model is not a linear " *
        "program (or it contains interval constraints).",
    )
    @test_throws err lp_sensitivity_report(model)
    MOI.empty!(m)
    MOI.Utilities.loadfromstring!(
        m,
        """
variables: x, y
minobjective: 1.0 * x + 1.0 * y
c: [x, y] in Nonnegatives(2)
""",
    )
    err = ErrorException(
        "Unable to compute LP sensitivity because model is not a linear " *
        "program (or it contains interval constraints).",
    )
    @test_throws err lp_sensitivity_report(model)
    return
end

function test_Degeneracy()
    m = MOIU.MockOptimizer(
        MOIU.Model{Float64}();
        eval_variable_constraint_dual = false,
    )
    model = direct_model(m)
    MOI.Utilities.loadfromstring!(
        m,
        """
        variables: x, y
        maxobjective: 2.0 * x + 4.0 * y
        c1: x + 2 * y <= 4.0
        c2: x + 2 * y <= 4.0
        """,
    )
    optimize!(model)
    MOI.set(m, MOI.TerminationStatus(), MOI.OPTIMAL)
    MOI.set(m, MOI.PrimalStatus(), MOI.FEASIBLE_POINT)
    MOI.set(m, MOI.DualStatus(), MOI.FEASIBLE_POINT)
    for (key, val) in Dict("x" => 0, "y" => 2, "c1" => 0, "c2" => 2)
        var = variable_by_name(model, key)
        if var !== nothing
            MOI.set(m, MOI.VariablePrimal(), index(var), val)
            MOI.set(m, MOI.VariableBasisStatus(), index(var), MOI.BASIC)
        else
            c = constraint_by_name(model, key)
            MOI.set(m, MOI.ConstraintDual(), index(c), val)
            MOI.set(
                m,
                MOI.ConstraintBasisStatus(),
                index(c),
                MOI.NONBASIC_AT_UPPER,
            )
        end
    end
    @test_throws(LinearAlgebra.SingularException, lp_sensitivity_report(model))
    return
end

function test_Simple_bounds_Min()
    _test_sensitivity(
        """
        variables: x
        minobjective: x
        x >= -1.0
        x <= 1.0
        """,
        Dict(
            "x" => TestSensitivitySolution(
                -1.0,
                NaN,
                MOI.NONBASIC_AT_LOWER,
                (-1, Inf),
            ),
            ("x", MOI.GreaterThan{Float64}) =>
                TestSensitivitySolution(NaN, 1.0, MOI.NONBASIC, (-Inf, 2.0)),
            ("x", MOI.LessThan{Float64}) =>
                TestSensitivitySolution(NaN, 0.0, MOI.BASIC, (-2, Inf)),
        ),
    )
    return
end

function test_Simple_bounds_Max()
    _test_sensitivity(
        """
        variables: x
        maxobjective: x
        x >= -1.0
        x <= 1.0
        """,
        Dict(
            "x" => TestSensitivitySolution(
                1.0,
                NaN,
                MOI.NONBASIC_AT_UPPER,
                (-1, Inf),
            ),
            ("x", MOI.GreaterThan{Float64}) =>
                TestSensitivitySolution(NaN, 0.0, MOI.BASIC, (-Inf, 2.0)),
            ("x", MOI.LessThan{Float64}) =>
                TestSensitivitySolution(NaN, -1.0, MOI.NONBASIC, (-2, Inf)),
        ),
    )
    return
end

function test_Max_I()
    _test_sensitivity(
        """
        variables: x, y, z, w
        maxobjective: 1.1 * x + y
        x >= -1.0
        x <= 1.0
        y >= 0.0
        z == 1.0
        c1: x + y + z + w == 1.0
        c2: x + y         <= 2.0
        """,
        Dict(
            "x" => TestSensitivitySolution(
                1.0,
                NaN,
                MOI.NONBASIC_AT_UPPER,
                (-0.1, Inf),
            ),
            "y" => TestSensitivitySolution(1.0, NaN, MOI.BASIC, (-1, 0.1)),
            "z" => TestSensitivitySolution(1.0, NaN, MOI.NONBASIC, (-Inf, Inf)),
            "w" => TestSensitivitySolution(-2.0, NaN, MOI.BASIC, (-Inf, 1.0)),
            ("x", MOI.GreaterThan{Float64}) =>
                TestSensitivitySolution(NaN, 0.0, MOI.BASIC, (-Inf, 2.0)),
            ("x", MOI.LessThan{Float64}) =>
                TestSensitivitySolution(NaN, -0.1, MOI.NONBASIC, (-2, 1.0)),
            ("y", MOI.GreaterThan{Float64}) =>
                TestSensitivitySolution(NaN, 0.0, MOI.BASIC, (-Inf, 1.0)),
            "c1" =>
                TestSensitivitySolution(NaN, 0.0, MOI.NONBASIC, (-Inf, Inf)),
            "c2" =>
                TestSensitivitySolution(NaN, -1.0, MOI.NONBASIC, (-1.0, Inf)),
            ("z", MOI.EqualTo{Float64}) =>
                TestSensitivitySolution(NaN, 0.0, MOI.NONBASIC, (-Inf, Inf)),
        ),
    )
    return
end

function test_Min_I()
    _test_sensitivity(
        """
        variables: x, y, z, w
        minobjective: -1.1 * x + -1.0 * y
        x >= -1.0
        x <= 1.0
        y >= 0.0
        z == 1.0
        c1: x + y + z + w == 1.0
        c2: x + y         <= 2.0
        """,
        Dict(
            "x" => TestSensitivitySolution(
                1.0,
                NaN,
                MOI.NONBASIC_AT_UPPER,
                (-Inf, 0.1),
            ),
            "y" => TestSensitivitySolution(1.0, NaN, MOI.BASIC, (-0.1, 1)),
            "z" => TestSensitivitySolution(1.0, NaN, MOI.NONBASIC, (-Inf, Inf)),
            "w" => TestSensitivitySolution(-2.0, NaN, MOI.BASIC, (-1.0, Inf)),
            ("x", MOI.GreaterThan{Float64}) =>
                TestSensitivitySolution(NaN, 0.0, MOI.BASIC, (-Inf, 2.0)),
            ("x", MOI.LessThan{Float64}) =>
                TestSensitivitySolution(NaN, -0.1, MOI.NONBASIC, (-2, 1.0)),
            ("y", MOI.GreaterThan{Float64}) =>
                TestSensitivitySolution(NaN, 0.0, MOI.BASIC, (-Inf, 1.0)),
            "c1" =>
                TestSensitivitySolution(NaN, 0.0, MOI.NONBASIC, (-Inf, Inf)),
            "c2" =>
                TestSensitivitySolution(NaN, -1.0, MOI.NONBASIC, (-1.0, Inf)),
            ("z", MOI.EqualTo{Float64}) =>
                TestSensitivitySolution(NaN, 0.0, MOI.NONBASIC, (-Inf, Inf)),
        ),
    )
    return
end

function test_Max_II()
    _test_sensitivity(
        """
        variables: x, y
        maxobjective: -1.0 * x + -1.0 * y
        x >= 0.0
        y >= 0.0
        c1l: x + 2 * y >= -1.0
        c1u: x + 2 * y <= 2.0
        c2: x + y      >= 0.5
        c3: 2 * x + y  <= 2.0
        """,
        Dict(
            "x" => TestSensitivitySolution(0.5, NaN, MOI.BASIC, (0, 1)),
            "y" => TestSensitivitySolution(
                0.0,
                NaN,
                MOI.NONBASIC_AT_LOWER,
                (-Inf, 0),
            ),
            ("x", MOI.GreaterThan{Float64}) =>
                TestSensitivitySolution(NaN, 0.0, MOI.BASIC, (-Inf, 0.5)),
            ("y", MOI.GreaterThan{Float64}) =>
                TestSensitivitySolution(NaN, 0.0, MOI.NONBASIC, (-1, 0.5)),
            "c1l" => TestSensitivitySolution(NaN, 0.0, MOI.BASIC, (-Inf, 1.5)),
            "c1u" => TestSensitivitySolution(NaN, 0.0, MOI.BASIC, (-1.5, Inf)),
            "c2" =>
                TestSensitivitySolution(NaN, 1.0, MOI.NONBASIC, (-0.5, 0.5)),
            "c3" => TestSensitivitySolution(NaN, 1.0, MOI.BASIC, (-1.0, Inf)),
        ),
    )
    return
end

function test_Min_II()
    _test_sensitivity(
        """
        variables: x, y
        minobjective: 1.0 * x + 1.0 * y
        x >= 0.0
        y >= 0.0
        c1l: x + 2 * y >= -1.0
        c1u: x + 2 * y <= 2.0
        c2: x + y      >= 0.5
        c3: 2 * x + y  <= 2.0
        """,
        Dict(
            "x" => TestSensitivitySolution(0.5, NaN, MOI.BASIC, (-1, 0)),
            "y" => TestSensitivitySolution(
                0.0,
                NaN,
                MOI.NONBASIC_AT_LOWER,
                (0, Inf),
            ),
            ("x", MOI.GreaterThan{Float64}) =>
                TestSensitivitySolution(NaN, 0.0, MOI.BASIC, (-Inf, 0.5)),
            ("y", MOI.GreaterThan{Float64}) =>
                TestSensitivitySolution(NaN, 0.0, MOI.NONBASIC, (-1, 0.5)),
            "c1l" => TestSensitivitySolution(NaN, 0.0, MOI.BASIC, (-Inf, 1.5)),
            "c1u" => TestSensitivitySolution(NaN, 0.0, MOI.BASIC, (-1.5, Inf)),
            "c2" =>
                TestSensitivitySolution(NaN, 1.0, MOI.NONBASIC, (-0.5, 0.5)),
            "c3" => TestSensitivitySolution(NaN, 1.0, MOI.BASIC, (-1.0, Inf)),
        ),
    )
    return
end

function test_Max_III()
    _test_sensitivity(
        """
        variables: x, y
        maxobjective: 6.0 * x + 4.0 * y
        x >= 0.0
        y >= 0.0
        c1: 1 * x + 1 * y <=  6.0
        c2: 2 * x + 1 * y <=  9.0
        c3: 2 * x + 3 * y <= 16.0
        """,
        Dict(
            "x" => TestSensitivitySolution(3.0, NaN, MOI.BASIC, (-2, 2)),
            "y" => TestSensitivitySolution(3.0, NaN, MOI.BASIC, (-1, 2)),
            ("x", MOI.GreaterThan{Float64}) =>
                TestSensitivitySolution(NaN, 0.0, MOI.BASIC, (-Inf, 3)),
            ("y", MOI.GreaterThan{Float64}) =>
                TestSensitivitySolution(NaN, 0.0, MOI.BASIC, (-Inf, 3)),
            "c1" =>
                TestSensitivitySolution(NaN, -2.0, MOI.NONBASIC, (-1.5, 0.25)),
            "c2" => TestSensitivitySolution(NaN, -2.0, MOI.NONBASIC, (-1, 3)),
            "c3" => TestSensitivitySolution(NaN, 0.0, MOI.BASIC, (-1, Inf)),
        ),
    )
    return
end

function test_Min_III()
    _test_sensitivity(
        """
        variables: x, y
        minobjective: -6.0 * x + -4.0 * y
        x >= 0.0
        y >= 0.0
        c1: 1 * x + 1 * y <=  6.0
        c2: 2 * x + 1 * y <=  9.0
        c3: 2 * x + 3 * y <= 16.0
        """,
        Dict(
            "x" => TestSensitivitySolution(3.0, NaN, MOI.BASIC, (-2, 2)),
            "y" => TestSensitivitySolution(3.0, NaN, MOI.BASIC, (-2, 1)),
            ("x", MOI.GreaterThan{Float64}) =>
                TestSensitivitySolution(NaN, 0.0, MOI.BASIC, (-Inf, 3)),
            ("y", MOI.GreaterThan{Float64}) =>
                TestSensitivitySolution(NaN, 0.0, MOI.BASIC, (-Inf, 3)),
            "c1" =>
                TestSensitivitySolution(NaN, -2.0, MOI.NONBASIC, (-1.5, 0.25)),
            "c2" => TestSensitivitySolution(NaN, -2.0, MOI.NONBASIC, (-1, 3)),
            "c3" => TestSensitivitySolution(NaN, 0.0, MOI.BASIC, (-1, Inf)),
        ),
    )
    return
end

function test_Max_IV()
    _test_sensitivity(
        """
        variables: x, y
        maxobjective: 1.0 * x + 1.0 * y
        x >= 0.0
        y >= 0.0
        c1l: x + 2 * y >= -1.0
        c1u: x + 2 * y <= 2.0
        c2: x + y      >= 0.5
        c3: 2 * x + y  <= 2.0
        """,
        Dict(
            "x" => TestSensitivitySolution(2 / 3, NaN, MOI.BASIC, (-0.5, 1)),
            "y" => TestSensitivitySolution(2 / 3, NaN, MOI.BASIC, (-0.5, 1)),
            ("x", MOI.GreaterThan{Float64}) =>
                TestSensitivitySolution(NaN, 0.0, MOI.BASIC, (-Inf, 2 / 3)),
            ("y", MOI.GreaterThan{Float64}) =>
                TestSensitivitySolution(NaN, 0.0, MOI.BASIC, (-Inf, 2 / 3)),
            "c1l" => TestSensitivitySolution(NaN, 0.0, MOI.BASIC, (-Inf, 3)),
            "c1u" =>
                TestSensitivitySolution(NaN, -1 / 3, MOI.NONBASIC, (-1, 2)),
            "c2" => TestSensitivitySolution(NaN, 0.0, MOI.BASIC, (-Inf, 5 / 6)),
            "c3" =>
                TestSensitivitySolution(NaN, -1 / 3, MOI.NONBASIC, (-1.0, 2)),
        ),
    )
    return
end

function test_Min_IV()
    _test_sensitivity(
        """
        variables: x, y
        minobjective: -1.0 * x + -1.0 * y
        x >= 0.0
        y >= 0.0
        c1l: x + 2 * y >= -1.0
        c1u: x + 2 * y <= 2.0
        c2: x + y      >= 0.5
        c3: 2 * x + y  <= 2.0
        """,
        Dict(
            "x" => TestSensitivitySolution(2 / 3, NaN, MOI.BASIC, (-1, 0.5)),
            "y" => TestSensitivitySolution(2 / 3, NaN, MOI.BASIC, (-1, 0.5)),
            ("x", MOI.GreaterThan{Float64}) =>
                TestSensitivitySolution(NaN, 0.0, MOI.BASIC, (-Inf, 2 / 3)),
            ("y", MOI.GreaterThan{Float64}) =>
                TestSensitivitySolution(NaN, 0.0, MOI.BASIC, (-Inf, 2 / 3)),
            "c1l" => TestSensitivitySolution(NaN, 0.0, MOI.BASIC, (-Inf, 3)),
            "c1u" =>
                TestSensitivitySolution(NaN, -1 / 3, MOI.NONBASIC, (-1, 2)),
            "c2" => TestSensitivitySolution(NaN, 0.0, MOI.BASIC, (-Inf, 5 / 6)),
            "c3" =>
                TestSensitivitySolution(NaN, -1 / 3, MOI.NONBASIC, (-1.0, 2)),
        ),
    )
    return
end

function test_Free_variable()
    _test_sensitivity(
        """
        variables: x, y
        minobjective: 1.0 * x
        x >= 1.0
        """,
        Dict(
            "x" => TestSensitivitySolution(
                1.0,
                NaN,
                MOI.NONBASIC_AT_LOWER,
                (-1.0, Inf),
            ),
            "y" =>
                TestSensitivitySolution(0.0, NaN, MOI.SUPER_BASIC, (0.0, 0.0)),
            ("x", MOI.GreaterThan{Float64}) => TestSensitivitySolution(
                NaN,
                1.0,
                MOI.NONBASIC_AT_LOWER,
                (-Inf, Inf),
            ),
        ),
    )
    return
end

end  # module

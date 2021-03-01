#  Copyright 2017, Iain Dunning, Joey Huchette, Miles Lubin, and contributors
#  This Source Code Form is subject to the terms of the Mozilla Public
#  License, v. 2.0. If a copy of the MPL was not distributed with this
#  file, You can obtain one at http://mozilla.org/MPL/2.0/.

using Test

struct TestSensitivitySolution
    primal::Float64
    dual::Float64
    basis::Union{Nothing,MOI.BasisStatusCode}
    range::Tuple{Float64,Float64}
end

function _test_sensitivity(model_string, solution)
    m = MOIU.MockOptimizer(
        MOIU.Model{Float64}(),
        eval_variable_constraint_dual = false,
    )
    model = direct_model(m)
    MOI.Utilities.loadfromstring!(m, model_string)
    optimize!(model)
    MOI.set(m, MOI.TerminationStatus(), MOI.OPTIMAL)
    MOI.set(m, MOI.PrimalStatus(), MOI.FEASIBLE_POINT)
    MOI.set(m, MOI.DualStatus(), MOI.FEASIBLE_POINT)
    obj_map = Dict{String,Any}()
    for (key, val) in solution
        var = variable_by_name(model, key)
        if var !== nothing
            obj_map[key] = var
            MOI.set(model, MOI.VariablePrimal(), var, val.primal)
            continue
        end
        c = constraint_by_name(model, key)
        @assert c !== nothing
        obj_map[key] = c
        MOI.set(model, MOI.ConstraintDual(), c, val.dual)
        MOI.set(model, MOI.ConstraintBasisStatus(), c, val.basis)
    end
    sens = lp_sensitivity_report(model)
    @testset "$(s_key)" for (s_key, val) in solution
        key = obj_map[s_key]
        @test sens[key][1] ≈ val.range[1]
        @test sens[key][2] ≈ val.range[2]
    end
end

@testset "Error handling" begin
    m = MOIU.MockOptimizer(
        MOIU.Model{Float64}(),
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
end

@testset "Problems" begin
    @testset "Simple bounds: Min" begin
        _test_sensitivity(
            """
            variables: x
            minobjective: x
            xlb: x >= -1.0
            xub: x <= 1.0
            """,
            Dict(
                "x" => TestSensitivitySolution(-1.0, NaN, nothing, (-1, Inf)),
                "xlb" => TestSensitivitySolution(
                    NaN,
                    1.0,
                    MOI.NONBASIC,
                    (-Inf, 2.0),
                ),
                "xub" =>
                    TestSensitivitySolution(NaN, 0.0, MOI.BASIC, (-2, Inf)),
            ),
        )
    end
    @testset "Simple bounds: Max" begin
        _test_sensitivity(
            """
            variables: x
            maxobjective: x
            xlb: x >= -1.0
            xub: x <= 1.0
            """,
            Dict(
                "x" => TestSensitivitySolution(1.0, NaN, nothing, (-1, Inf)),
                "xlb" =>
                    TestSensitivitySolution(NaN, 0.0, MOI.BASIC, (-Inf, 2.0)),
                "xub" =>
                    TestSensitivitySolution(NaN, -1.0, MOI.NONBASIC, (-2, Inf)),
            ),
        )
    end
    @testset "Max I" begin
        _test_sensitivity(
            """
            variables: x, y, z, w
            maxobjective: 1.1 * x + y
            xlb: x >= -1.0
            xub: x <= 1.0
            ylb: y >= 0.0
            zfx: z == 1.0
            c1: x + y + z + w == 1.0
            c2: x + y         <= 2.0
            """,
            Dict(
                "x" => TestSensitivitySolution(1.0, NaN, nothing, (-0.1, Inf)),
                "y" => TestSensitivitySolution(1.0, NaN, nothing, (-1, 0.1)),
                "z" => TestSensitivitySolution(1.0, NaN, nothing, (-Inf, Inf)),
                "w" => TestSensitivitySolution(-2.0, NaN, nothing, (-Inf, 1.0)),
                "xlb" =>
                    TestSensitivitySolution(NaN, 0.0, MOI.BASIC, (-Inf, 2.0)),
                "xub" =>
                    TestSensitivitySolution(NaN, -0.1, MOI.NONBASIC, (-2, 1.0)),
                "ylb" =>
                    TestSensitivitySolution(NaN, 0.0, MOI.BASIC, (-Inf, 1.0)),
                "c1" => TestSensitivitySolution(
                    NaN,
                    0.0,
                    MOI.NONBASIC,
                    (-Inf, Inf),
                ),
                "c2" => TestSensitivitySolution(
                    NaN,
                    -1.0,
                    MOI.NONBASIC,
                    (-1.0, Inf),
                ),
                "zfx" => TestSensitivitySolution(
                    NaN,
                    0.0,
                    MOI.NONBASIC,
                    (-Inf, Inf),
                ),
            ),
        )
    end
    @testset "Min I" begin
        _test_sensitivity(
            """
            variables: x, y, z, w
            minobjective: -1.1 * x + -1.0 * y
            xlb: x >= -1.0
            xub: x <= 1.0
            ylb: y >= 0.0
            zfx: z == 1.0
            c1: x + y + z + w == 1.0
            c2: x + y         <= 2.0
            """,
            Dict(
                "x" => TestSensitivitySolution(1.0, NaN, nothing, (-Inf, 0.1)),
                "y" => TestSensitivitySolution(1.0, NaN, nothing, (-0.1, 1)),
                "z" => TestSensitivitySolution(1.0, NaN, nothing, (-Inf, Inf)),
                "w" => TestSensitivitySolution(-2.0, NaN, nothing, (-1.0, Inf)),
                "xlb" =>
                    TestSensitivitySolution(NaN, 0.0, MOI.BASIC, (-Inf, 2.0)),
                "xub" =>
                    TestSensitivitySolution(NaN, -0.1, MOI.NONBASIC, (-2, 1.0)),
                "ylb" =>
                    TestSensitivitySolution(NaN, 0.0, MOI.BASIC, (-Inf, 1.0)),
                "c1" => TestSensitivitySolution(
                    NaN,
                    0.0,
                    MOI.NONBASIC,
                    (-Inf, Inf),
                ),
                "c2" => TestSensitivitySolution(
                    NaN,
                    -1.0,
                    MOI.NONBASIC,
                    (-1.0, Inf),
                ),
                "zfx" => TestSensitivitySolution(
                    NaN,
                    0.0,
                    MOI.NONBASIC,
                    (-Inf, Inf),
                ),
            ),
        )
    end
    @testset "Max II" begin
        _test_sensitivity(
            """
            variables: x, y
            maxobjective: -1.0 * x + -1.0 * y
            xlb: x >= 0.0
            ylb: y >= 0.0
            c1l: x + 2 * y >= -1.0
            c1u: x + 2 * y <= 2.0
            c2: x + y      >= 0.5
            c3: 2 * x + y  <= 2.0
            """,
            Dict(
                "x" => TestSensitivitySolution(0.5, NaN, nothing, (0, 1)),
                "y" => TestSensitivitySolution(0.0, NaN, nothing, (-Inf, 0)),
                "xlb" =>
                    TestSensitivitySolution(NaN, 0.0, MOI.BASIC, (-Inf, 0.5)),
                "ylb" =>
                    TestSensitivitySolution(NaN, 0.0, MOI.NONBASIC, (-1, 0.5)),
                "c1l" =>
                    TestSensitivitySolution(NaN, 0.0, MOI.BASIC, (-Inf, 1.5)),
                "c1u" =>
                    TestSensitivitySolution(NaN, 0.0, MOI.BASIC, (-1.5, Inf)),
                "c2" => TestSensitivitySolution(
                    NaN,
                    1.0,
                    MOI.NONBASIC,
                    (-0.5, 0.5),
                ),
                "c3" =>
                    TestSensitivitySolution(NaN, 1.0, MOI.BASIC, (-1.0, Inf)),
            ),
        )
    end
    @testset "Min II" begin
        _test_sensitivity(
            """
            variables: x, y
            minobjective: 1.0 * x + 1.0 * y
            xlb: x >= 0.0
            ylb: y >= 0.0
            c1l: x + 2 * y >= -1.0
            c1u: x + 2 * y <= 2.0
            c2: x + y      >= 0.5
            c3: 2 * x + y  <= 2.0
            """,
            Dict(
                "x" => TestSensitivitySolution(0.5, NaN, nothing, (-1, 0)),
                "y" => TestSensitivitySolution(0.0, NaN, nothing, (0, Inf)),
                "xlb" =>
                    TestSensitivitySolution(NaN, 0.0, MOI.BASIC, (-Inf, 0.5)),
                "ylb" =>
                    TestSensitivitySolution(NaN, 0.0, MOI.NONBASIC, (-1, 0.5)),
                "c1l" =>
                    TestSensitivitySolution(NaN, 0.0, MOI.BASIC, (-Inf, 1.5)),
                "c1u" =>
                    TestSensitivitySolution(NaN, 0.0, MOI.BASIC, (-1.5, Inf)),
                "c2" => TestSensitivitySolution(
                    NaN,
                    1.0,
                    MOI.NONBASIC,
                    (-0.5, 0.5),
                ),
                "c3" =>
                    TestSensitivitySolution(NaN, 1.0, MOI.BASIC, (-1.0, Inf)),
            ),
        )
    end
    @testset "Max III" begin
        _test_sensitivity(
            """
            variables: x, y
            maxobjective: 6.0 * x + 4.0 * y
            xlb: x >= 0.0
            ylb: y >= 0.0
            c1: 1 * x + 1 * y <=  6.0
            c2: 2 * x + 1 * y <=  9.0
            c3: 2 * x + 3 * y <= 16.0
            """,
            Dict(
                "x" => TestSensitivitySolution(3.0, NaN, nothing, (-2, 2)),
                "y" => TestSensitivitySolution(3.0, NaN, nothing, (-1, 2)),
                "xlb" =>
                    TestSensitivitySolution(NaN, 0.0, MOI.BASIC, (-Inf, 3)),
                "ylb" =>
                    TestSensitivitySolution(NaN, 0.0, MOI.BASIC, (-Inf, 3)),
                "c1" => TestSensitivitySolution(
                    NaN,
                    -2.0,
                    MOI.NONBASIC,
                    (-1.5, 0.25),
                ),
                "c2" =>
                    TestSensitivitySolution(NaN, -2.0, MOI.NONBASIC, (-1, 3)),
                "c3" => TestSensitivitySolution(NaN, 0.0, MOI.BASIC, (-1, Inf)),
            ),
        )
    end
    @testset "Min III" begin
        _test_sensitivity(
            """
            variables: x, y
            minobjective: -6.0 * x + -4.0 * y
            xlb: x >= 0.0
            ylb: y >= 0.0
            c1: 1 * x + 1 * y <=  6.0
            c2: 2 * x + 1 * y <=  9.0
            c3: 2 * x + 3 * y <= 16.0
            """,
            Dict(
                "x" => TestSensitivitySolution(3.0, NaN, nothing, (-2, 2)),
                "y" => TestSensitivitySolution(3.0, NaN, nothing, (-2, 1)),
                "xlb" =>
                    TestSensitivitySolution(NaN, 0.0, MOI.BASIC, (-Inf, 3)),
                "ylb" =>
                    TestSensitivitySolution(NaN, 0.0, MOI.BASIC, (-Inf, 3)),
                "c1" => TestSensitivitySolution(
                    NaN,
                    -2.0,
                    MOI.NONBASIC,
                    (-1.5, 0.25),
                ),
                "c2" =>
                    TestSensitivitySolution(NaN, -2.0, MOI.NONBASIC, (-1, 3)),
                "c3" => TestSensitivitySolution(NaN, 0.0, MOI.BASIC, (-1, Inf)),
            ),
        )
    end
    @testset "Max IV" begin
        _test_sensitivity(
            """
            variables: x, y
            maxobjective: 1.0 * x + 1.0 * y
            xlb: x >= 0.0
            ylb: y >= 0.0
            c1l: x + 2 * y >= -1.0
            c1u: x + 2 * y <= 2.0
            c2: x + y      >= 0.5
            c3: 2 * x + y  <= 2.0
            """,
            Dict(
                "x" => TestSensitivitySolution(2 / 3, NaN, nothing, (-0.5, 1)),
                "y" => TestSensitivitySolution(2 / 3, NaN, nothing, (-0.5, 1)),
                "xlb" =>
                    TestSensitivitySolution(NaN, 0.0, MOI.BASIC, (-Inf, 2 / 3)),
                "ylb" =>
                    TestSensitivitySolution(NaN, 0.0, MOI.BASIC, (-Inf, 2 / 3)),
                "c1l" =>
                    TestSensitivitySolution(NaN, 0.0, MOI.BASIC, (-Inf, 3)),
                "c1u" =>
                    TestSensitivitySolution(NaN, -1 / 3, MOI.NONBASIC, (-1, 2)),
                "c2" =>
                    TestSensitivitySolution(NaN, 0.0, MOI.BASIC, (-Inf, 5 / 6)),
                "c3" => TestSensitivitySolution(
                    NaN,
                    -1 / 3,
                    MOI.NONBASIC,
                    (-1.0, 2),
                ),
            ),
        )
    end
    @testset "Min IV" begin
        _test_sensitivity(
            """
            variables: x, y
            minobjective: -1.0 * x + -1.0 * y
            xlb: x >= 0.0
            ylb: y >= 0.0
            c1l: x + 2 * y >= -1.0
            c1u: x + 2 * y <= 2.0
            c2: x + y      >= 0.5
            c3: 2 * x + y  <= 2.0
            """,
            Dict(
                "x" => TestSensitivitySolution(2 / 3, NaN, nothing, (-1, 0.5)),
                "y" => TestSensitivitySolution(2 / 3, NaN, nothing, (-1, 0.5)),
                "xlb" =>
                    TestSensitivitySolution(NaN, 0.0, MOI.BASIC, (-Inf, 2 / 3)),
                "ylb" =>
                    TestSensitivitySolution(NaN, 0.0, MOI.BASIC, (-Inf, 2 / 3)),
                "c1l" =>
                    TestSensitivitySolution(NaN, 0.0, MOI.BASIC, (-Inf, 3)),
                "c1u" =>
                    TestSensitivitySolution(NaN, -1 / 3, MOI.NONBASIC, (-1, 2)),
                "c2" =>
                    TestSensitivitySolution(NaN, 0.0, MOI.BASIC, (-Inf, 5 / 6)),
                "c3" => TestSensitivitySolution(
                    NaN,
                    -1 / 3,
                    MOI.NONBASIC,
                    (-1.0, 2),
                ),
            ),
        )
    end
end

#  Copyright 2017, Iain Dunning, Joey Huchette, Miles Lubin, and contributors
#  This Source Code Form is subject to the terms of the Mozilla Public
#  License, v. 2.0. If a copy of the MPL was not distributed with this
#  file, You can obtain one at https://mozilla.org/MPL/2.0/.
#############################################################################
# JuMP
# An algebraic modeling language for Julia
# See https://github.com/jump-dev/JuMP.jl
#############################################################################

module TestConstraint

using JuMP
using Test

import LinearAlgebra
import SparseArrays

include(joinpath(@__DIR__, "utilities.jl"))

function _test_constraint_name_util(constraint, s_name, F::Type, S::Type)
    @test s_name == @inferred name(constraint)
    model = constraint.model
    @test constraint.index == constraint_by_name(model, s_name).index
    if model isa Model
        @test constraint.index == constraint_by_name(model, s_name, F, S).index
    end
    return
end

function test_extension_VariableIndex_constraints(
    ModelType = Model,
    VariableRefType = VariableRef,
)
    T = value_type(ModelType)
    m = ModelType()
    @variable(m, x)
    # x <= 10.0 doesn't translate to a SingleVariable constraint because
    # the LHS is first subtracted to form x - 10.0 <= 0.
    @constraint(m, cref, x in MOI.LessThan{T}(10))
    c = constraint_object(cref)
    @test c.func == x
    @test c.set == MOI.LessThan(T(10))
    @variable(m, y[1:2])
    @constraint(m, cref2[i = 1:2], y[i] in MOI.LessThan{T}(i))
    c = constraint_object(cref2[1])
    @test c.func == y[1]
    @test c.set == MOI.LessThan(T(1))
    return
end

function test_extension_Container_constraints(
    ModelType = Model,
    VariableRefType = VariableRef,
)
    m = ModelType()
    S = ["a", "b"]
    @variable(m, x[S], Bin)
    cref = @constraint(m, x in SOS2())
    c = constraint_object(cref)
    @test c.func == [x["a"], x["b"]]
    return
end

function test_extension_VectorOfVariables_constraints(
    ModelType = Model,
    VariableRefType = VariableRef,
)
    m = ModelType()
    @variable(m, x[1:2])
    cref = @constraint(m, x in MOI.Zeros(2))
    c = constraint_object(cref)
    @test c.func == x
    @test c.set == MOI.Zeros(2)
    cref = @constraint(m, [x[2], x[1]] in MOI.Zeros(2))
    c = constraint_object(cref)
    @test c.func == [x[2], x[1]]
    @test c.set == MOI.Zeros(2)
    @test_throws DimensionMismatch @constraint(m, [x[2], x[1]] in MOI.Zeros(3))
    return
end

function test_extension_AffExpr_scalar_constraints(
    ModelType = Model,
    VariableRefType = VariableRef,
)
    AffExprType = JuMP.GenericAffExpr{value_type(ModelType),VariableRefType}
    model = ModelType()
    T = value_type(ModelType)
    @variable(model, x)
    cref = @constraint(model, 2x <= 10)
    @test "" == @inferred name(cref)
    set_name(cref, "c")
    _test_constraint_name_util(cref, "c", AffExprType, MOI.LessThan{Float64})
    c = constraint_object(cref)
    @test isequal_canonical(c.func, 2x)
    @test c.set == MOI.LessThan(T(10))
    cref = @constraint(model, 3x + 1 ≥ 10)
    c = constraint_object(cref)
    @test isequal_canonical(c.func, 3x)
    @test c.set == MOI.GreaterThan(T(9))
    cref = @constraint(model, 1 == -x)
    c = constraint_object(cref)
    @test isequal_canonical(c.func, one(T) * x)
    @test c.set == MOI.EqualTo(-one(T))
    cref = @constraint(model, 2 == 1)
    c = constraint_object(cref)
    @test isequal_canonical(c.func, zero(AffExpr))
    @test c.set == MOI.EqualTo(-one(T))
    return
end

function test_extension_AffExpr_vectorized_constraints(
    ModelType = Model,
    VariableRefType = VariableRef,
)
    T = value_type(ModelType)
    model = ModelType()
    @variable(model, x)
    err = ErrorException(
        "In `@constraint(model, [x, 2x] in MOI.EqualTo(1.0))`: " *
        "Unexpected vector in scalar constraint. The left- and right-hand " *
        "sides of the constraint must have the same dimension.",
    )
    @test_throws_strip err @constraint(model, [x, 2x] in MOI.EqualTo(1.0))
    VT = typeof([x, 2x])
    err = ErrorException(
        "Operation `sub_mul` between `$VT` and `$Int` is not " *
        "allowed. This most often happens when you write a constraint like " *
        "`x >= y` where `x` is an array and `y` is a constant. Use the " *
        "broadcast syntax `x .- y >= 0` instead.",
    )
    @test_throws err @constraint(model, [x, 2x] == 1)
    err = ErrorException(
        "Operation `sub_mul` between `$Int` and `$VT` is not " *
        "allowed. This most often happens when you write a constraint like " *
        "`x >= y` where `x` is a constant and `y` is an array. Use the " *
        "broadcast syntax `x .- y >= 0` instead.",
    )
    a = 1
    @test_throws err @constraint(model, a == [x, 2x])
    @test_macro_throws ErrorException begin
        @constraint(model, [x == 1 - x, 2x == 3])
    end
    cref = @constraint(model, [x, 2x] .== [1 - x, 3])
    c = constraint_object.(cref)
    @test isequal_canonical(c[1].func, 2x)
    @test c[1].set == MOI.EqualTo(T(1))
    @test isequal_canonical(c[2].func, 2x)
    @test c[2].set == MOI.EqualTo(T(3))
    return
end

function test_extension_AffExpr_vectorized_interval_constraints(
    ModelType = Model,
    VariableRefType = VariableRef,
)
    T = value_type(ModelType)
    model = ModelType()
    @variable(model, x[1:2])
    err = ErrorException(
        "In `@constraint(model, b <= x <= b)`: Unexpected vectors in " *
        "scalar constraint. Did you mean to use the dot comparison " *
        "operators `l .<= f(x) .<= u` instead?",
    )
    b = T[5, 6]
    @test_throws_strip err @constraint(model, b <= x <= b)
    cref = @constraint(model, b .<= x .<= b)
    c = constraint_object.(cref)
    for i in 1:2
        @test isequal_canonical(c[i].func, 1 * x[i])
        @test c[i].set == MOI.Interval(b[i], b[i])
    end
    return
end

function test_extension_AffExpr_vector_constraints(
    ModelType = Model,
    VariableRefType = VariableRef,
)
    AffExprType = JuMP.GenericAffExpr{value_type(ModelType),VariableRefType}
    model = ModelType()
    cref = @constraint(model, [1, 2] in MOI.Zeros(2))
    c = constraint_object(cref)
    @test isequal_canonical(c.func[1], zero(AffExpr) + 1)
    @test isequal_canonical(c.func[2], zero(AffExpr) + 2)
    @test c.set == MOI.Zeros(2)
    @test c.shape isa VectorShape
    @test_throws DimensionMismatch @constraint(model, [1, 2] in MOI.Zeros(3))
    return
end

function test_extension_delete_constraints(
    ModelType = Model,
    VariableRefType = VariableRef,
)
    model = ModelType()
    @variable(model, x)
    constraint_ref = @constraint(model, 2x <= 1)
    @test is_valid(model, constraint_ref)
    delete(model, constraint_ref)
    @test !is_valid(model, constraint_ref)
    second_model = ModelType()
    @test_throws Exception delete(second_model, constraint_ref)
    return
end

function test_extension_batch_delete_constraints(
    ModelType = Model,
    VariableRefType = VariableRef,
)
    model = ModelType()
    @variable(model, x[1:9])
    cons = [@constraint(model, sum(x[1:2:9]) <= 3)]
    push!(cons, @constraint(model, sum(x[2:2:8]) <= 2))
    push!(cons, @constraint(model, sum(x[1:3:9]) <= 1))
    @test all(is_valid.(model, cons))
    delete(model, cons[[1, 3]])
    @test all((!is_valid).(model, cons[[1, 3]]))
    @test is_valid(model, cons[2])
    second_model = ModelType()
    @test_throws Exception delete(second_model, cons[[1, 3]])
    @test_throws Exception delete(second_model, [cons[2]])
    return
end

function test_extension_two_sided_constraints(
    ModelType = Model,
    VariableRefType = VariableRef,
)
    T = value_type(ModelType)
    AffExprType = JuMP.GenericAffExpr{T,VariableRefType}
    m = ModelType()
    @variable(m, x)
    @variable(m, y)
    @constraint(m, cref, 1 <= x + y + 1 <= 2)
    _test_constraint_name_util(cref, "cref", AffExprType, MOI.Interval{T})
    c = constraint_object(cref)
    @test isequal_canonical(c.func, x + y)
    @test c.set ==
          MOI.Interval(zero(value_type(ModelType)), one(value_type(ModelType)))
    cref = @constraint(m, 2x - y + T(2) ∈ MOI.Interval(-one(T), one(T)))
    @test "" == @inferred name(cref)
    c = constraint_object(cref)
    @test isequal_canonical(c.func, 2x - y)
    @test c.set == MOI.Interval(-3one(T), -one(T))
    return
end

function test_extension_broadcasted_constraint_eq(
    ModelType = Model,
    VariableRefType = VariableRef,
)
    T = value_type(ModelType)
    m = ModelType()
    @variable(m, x[1:2])
    A = T[1 2; 3 4]
    b = T[4, 5]
    cref = @constraint(m, A * x .== b)
    @test (2,) == @inferred size(cref)
    c1 = constraint_object(cref[1])
    @test isequal_canonical(c1.func, x[1] + 2x[2])
    @test c1.set == MOI.EqualTo(T(4))
    c2 = constraint_object(cref[2])
    @test isequal_canonical(c2.func, 3x[1] + 4x[2])
    @test c2.set == MOI.EqualTo(T(5))
    return
end

function test_extension_broadcasted_constraint_leq(
    ModelType = Model,
    VariableRefType = VariableRef,
)
    T = value_type(ModelType)
    m = ModelType()
    @variable(m, x[1:2, 1:2])
    UB = T[1 2; 3 4]
    cref = @constraint(m, x .+ 1 .<= UB)
    @test (2, 2) == @inferred size(cref)
    for i in 1:2
        for j in 1:2
            c = constraint_object(cref[i, j])
            @test isequal_canonical(c.func, x[i, j] + 0)
            @test c.set == MOI.LessThan(UB[i, j] - 1)
        end
    end
    return
end

function test_extension_broadcasted_two_sided_constraint(
    ModelType = Model,
    VariableRefType = VariableRef,
)
    T = value_type(ModelType)
    m = ModelType()
    @variable(m, x[1:2])
    @variable(m, y[1:2])
    l = T[1, 2]
    u = T[3, 4]
    cref = @constraint(m, l .<= x + y .+ 1 .<= u)
    @test (2,) == @inferred size(cref)
    for i in 1:2
        c = constraint_object(cref[i])
        @test isequal_canonical(c.func, x[i] + y[i])
        @test c.set == MOI.Interval(l[i] - 1, u[i] - 1)
    end
    return
end

function test_extension_broadcasted_constraint_with_indices(
    ModelType = Model,
    VariableRefType = VariableRef,
)
    m = ModelType()
    @variable m x[1:2]
    @constraint m cref1[i = 2:4] x .== [i, i + 1]
    ConstraintRefType = eltype(cref1[2])
    @test cref1 isa Containers.DenseAxisArray{Vector{ConstraintRefType}}
    @constraint m cref2[i = 1:3, j = 1:4] x .≤ [i + j, i - j]
    ConstraintRefType = eltype(cref2[1])
    @test cref2 isa Matrix{Vector{ConstraintRefType}}
    @variable m y[1:2, 1:2]
    @constraint m cref3[i = 1:2] y[i, :] .== 1
    ConstraintRefType = eltype(cref3[1])
    @test cref3 isa Vector{Vector{ConstraintRefType}}
    return
end

function test_extension_quadexpr_constraints(
    ModelType = Model,
    VariableRefType = VariableRef,
)
    T = value_type(ModelType)
    model = ModelType()
    @variable(model, x)
    @variable(model, y)

    cref = @constraint(model, x^2 + x <= 1)
    c = constraint_object(cref)
    @test isequal_canonical(c.func, x^2 + x)
    @test c.set == MOI.LessThan(one(T))

    cref = @constraint(model, y * x - 1 == 0)
    c = constraint_object(cref)
    @test isequal_canonical(c.func, x * y)
    @test c.set == MOI.EqualTo(one(T))

    cref = @constraint(
        model,
        [2x - 4x * y + 3x^2 - 1, -3y + 2x * y - 2x^2 + 1] in SecondOrderCone()
    )
    c = constraint_object(cref)
    @test isequal_canonical(c.func[1], -1 + 3x^2 - 4x * y + 2x)
    @test isequal_canonical(c.func[2], 1 - 2x^2 + 2x * y - 3y)
    @test c.set == MOI.SecondOrderCone(2)
    return
end

function test_extension_syntax_error_constraint(
    ModelType = Model,
    VariableRefType = VariableRef,
)
    model = ModelType()
    @variable(model, x[1:2])
    err = ErrorException(
        "In `@constraint(model, [3, x] in SecondOrderCone())`: Unable to " *
        "add the constraint because we don't recognize $([3, x]) as a " *
        "valid JuMP function.",
    )
    @test_throws_strip err @constraint(model, [3, x] in SecondOrderCone())
    return
end

function test_extension_indicator_constraint(
    ModelType = Model,
    VariableRefType = VariableRef,
)
    T = value_type(ModelType)
    model = ModelType()
    @variable(model, a, Bin)
    @variable(model, b, Bin)
    @variable(model, x)
    @variable(model, y)
    for cref in [
        @constraint(model, a => {x + 2y <= 1})
        @constraint(model, a ⇒ {x + 2y ≤ 1})
        @constraint(model, a --> {x + 2y ≤ 1})
    ]
        c = constraint_object(cref)
        @test c.func == [a, x + 2y]
        @test c.set == MOI.Indicator{MOI.ACTIVATE_ON_ONE}(MOI.LessThan(one(T)))
    end
    for cref in [
        @constraint(model, !b => {2x + y <= 1})
        @constraint(model, ¬b ⇒ {2x + y ≤ 1})
        @constraint(model, ¬b --> {2x + y ≤ 1})
        # This returns a vector of constraints that is concatenated.
        @constraint(model, ![b, b] .=> {[2x + y, 2x + y] .≤ 1})
    ]
        c = constraint_object(cref)
        @test c.func == [b, 2x + y]
        @test c.set == MOI.Indicator{MOI.ACTIVATE_ON_ZERO}(MOI.LessThan(one(T)))
    end
    err = ErrorException(
        "In `@constraint(model, !(a, b) => {x <= 1})`: Invalid binary variable expression `!(a, b)` for indicator constraint.",
    )
    @test_macro_throws err @constraint(model, !(a, b) => {x <= 1})
    err = ErrorException(
        "In `@constraint(model, a => x)`: Invalid right-hand side `x` of indicator constraint. Expected constraint surrounded by `{` and `}`.",
    )
    @test_macro_throws err @constraint(model, a => x)
    err = ErrorException(
        "In `@constraint(model, a => x <= 1)`: Invalid right-hand side `x <= 1` of indicator constraint. Expected constraint surrounded by `{` and `}`.",
    )
    @test_macro_throws err @constraint(model, a => x <= 1)
    err = ErrorException(
        "In `@constraint(model, a => {x <= 1, x >= 0})`: Invalid right-hand side `{x <= 1, x >= 0}` of indicator constraint. Expected constraint surrounded by `{` and `}`.",
    )
    @test_macro_throws err @constraint(model, a => {x <= 1, x >= 0})
    err = ErrorException(
        "In `@constraint(model, [a, b] .=> {x + y <= 1})`: Inconsistent use of `.` in symbols to indicate vectorization.",
    )
    @test_macro_throws err @constraint(model, [a, b] .=> {x + y <= 1})
    return
end

function test_extension_SDP_constraint(
    ModelType = Model,
    VariableRefType = VariableRef,
)
    T = value_type(ModelType)
    AffExprType = JuMP.GenericAffExpr{T,VariableRefType}
    m = ModelType()
    @variable(m, x)
    @variable(m, y)
    @variable(m, z)
    @variable(m, w)

    cref = @constraint(m, [x y; z w] in PSDCone())
    c = constraint_object(cref)
    @test c.func == [x, z, y, w]
    @test c.set == MOI.PositiveSemidefiniteConeSquare(2)
    @test c.shape isa SquareMatrixShape

    @constraint(
        m,
        sym_ref,
        LinearAlgebra.Symmetric([x 1; 1 -y] - [1 x; x -2]) in PSDCone()
    )
    _test_constraint_name_util(
        sym_ref,
        "sym_ref",
        Vector{AffExprType},
        MOI.PositiveSemidefiniteConeTriangle,
    )
    c = constraint_object(sym_ref)
    @test isequal_canonical(c.func[1], x - 1)
    @test isequal_canonical(c.func[2], 1 - x)
    @test isequal_canonical(c.func[3], 2 - y)
    @test c.set == MOI.PositiveSemidefiniteConeTriangle(2)
    @test c.shape isa SymmetricMatrixShape

    @constraint(m, cref, [x 1; 1 -y] >= [1 x; x -2], PSDCone())
    _test_constraint_name_util(
        cref,
        "cref",
        Vector{AffExprType},
        MOI.PositiveSemidefiniteConeSquare,
    )
    c = constraint_object(cref)
    @test isequal_canonical(c.func[1], x - 1)
    @test isequal_canonical(c.func[2], 1 - x)
    @test isequal_canonical(c.func[3], 1 - x)
    @test isequal_canonical(c.func[4], 2 - y)
    @test c.set == MOI.PositiveSemidefiniteConeSquare(2)
    @test c.shape isa SquareMatrixShape

    @constraint(m, iref[i = 1:2], 0 <= [x+i x+y; x+y -y], PSDCone())
    for i in 1:2
        _test_constraint_name_util(
            iref[i],
            "iref[$i]",
            Vector{AffExprType},
            MOI.PositiveSemidefiniteConeSquare,
        )
        c = constraint_object(iref[i])
        @test isequal_canonical(c.func[1], x + i)
        @test isequal_canonical(c.func[2], x + y)
        @test isequal_canonical(c.func[3], x + y)
        @test isequal_canonical(c.func[4], -y)
        @test c.set == MOI.PositiveSemidefiniteConeSquare(2)
        @test c.shape isa SquareMatrixShape
    end

    @constraint(m, con_d, 0 <= LinearAlgebra.Diagonal([x, y]), PSDCone())
    c = constraint_object(con_d)
    @test c.func isa Vector{AffExprType}
    @test isequal_canonical(c.func[1], 1x)
    @test iszero(c.func[2])
    @test iszero(c.func[3])
    @test isequal_canonical(c.func[4], 1y)
    @test c.set == MOI.PositiveSemidefiniteConeSquare(2)

    @constraint(
        m,
        con_d_sym,
        0 <= LinearAlgebra.Symmetric(LinearAlgebra.Diagonal([x, y])),
        PSDCone()
    )
    c = constraint_object(con_d_sym)
    @test c.func isa Vector{AffExprType}
    @test isequal_canonical(c.func[1], 1x)
    @test iszero(c.func[2])
    @test isequal_canonical(c.func[3], 1y)
    @test c.set == MOI.PositiveSemidefiniteConeTriangle(2)

    @constraint(
        m,
        con_td,
        LinearAlgebra.Tridiagonal([z], [x, y], [w]) >= 0,
        PSDCone()
    )
    c = constraint_object(con_td)
    @test c.func isa Vector{AffExprType}
    @test isequal_canonical(c.func[1], 1x)
    @test isequal_canonical(c.func[2], 1z)
    @test isequal_canonical(c.func[3], 1w)
    @test isequal_canonical(c.func[4], 1y)
    @test c.set == MOI.PositiveSemidefiniteConeSquare(2)

    @constraint(
        m,
        con_td_sym,
        LinearAlgebra.Symmetric(LinearAlgebra.Tridiagonal([z], [x, y], [w])) >=
        0,
        PSDCone(),
    )
    c = constraint_object(con_td_sym)
    @test c.func isa Vector{AffExprType}
    @test isequal_canonical(c.func[1], 1x)
    @test isequal_canonical(c.func[2], 1w)
    @test isequal_canonical(c.func[3], 1y)
    @test c.set == MOI.PositiveSemidefiniteConeTriangle(2)

    @constraint(
        m,
        con_ut,
        LinearAlgebra.UpperTriangular([x y; z w]) >= 0,
        PSDCone()
    )
    c = constraint_object(con_ut)
    @test c.func isa Vector{AffExprType}
    @test isequal_canonical(c.func[1], 1x)
    @test iszero(c.func[2])
    @test isequal_canonical(c.func[3], 1y)
    @test isequal_canonical(c.func[4], 1w)
    @test c.set == MOI.PositiveSemidefiniteConeSquare(2)

    @constraint(
        m,
        con_lt,
        0 <= LinearAlgebra.LowerTriangular([x y; z w]),
        PSDCone()
    )
    c = constraint_object(con_lt)
    @test c.func isa Vector{AffExprType}
    @test isequal_canonical(c.func[1], 1x)
    @test isequal_canonical(c.func[2], 1z)
    @test iszero(c.func[3])
    @test isequal_canonical(c.func[4], 1w)
    @test c.set == MOI.PositiveSemidefiniteConeSquare(2)
    return
end

function test_extension_SDP_errors(
    ModelType = Model,
    VariableRefType = VariableRef,
)
    AffExprType = JuMP.GenericAffExpr{value_type(ModelType),VariableRefType}
    model = ModelType()
    @variable(model, x)
    @variable(model, y)
    @variable(model, z)
    @variable(model, w)
    aff_str = "$AffExprType"
    err = ErrorException(
        "In `@constraint(model, [x 1; 1 -y] >= [1 x; x -2], PSDCone(), unknown_kw = 1)`:" *
        " Unrecognized constraint building format. Tried to invoke " *
        "`build_constraint(error, $(aff_str)[x - " *
        "1 -x + 1; -x + 1 -y + 2], $(MOI.GreaterThan(0.0)), $(PSDCone()); unknown_kw = 1)`, but no " *
        "such method exists. This is due to specifying an unrecognized " *
        "function, constraint set, and/or extra positional/keyword " *
        "arguments.\n\nIf you're trying to create a JuMP extension, you " *
        "need to implement `build_constraint` to accomodate these arguments.",
    )
    @test_throws_strip(
        err,
        @constraint(
            model,
            [x 1; 1 -y] >= [1 x; x -2],
            PSDCone(),
            unknown_kw = 1,
        ),
    )
    # Invalid sense == in SDP constraint
    @test_throws(
        ErrorException,
        @constraint(model, [x 1; 1 -y] == [1 x; x -2], PSDCone()),
    )
    return
end

function _test_constraint_name_util(ModelType, VariableRefType)
    model = ModelType()
    @variable(model, x)
    @constraint(model, con, x^2 == 1)
    _test_constraint_name_util(
        con,
        "con",
        GenericQuadExpr{Float64,VariableRefType},
        MOI.EqualTo{Float64},
    )
    set_name(con, "kon")
    @test constraint_by_name(model, "con") isa Nothing
    _test_constraint_name_util(
        con,
        "kon",
        GenericQuadExpr{Float64,VariableRefType},
        MOI.EqualTo{Float64},
    )
    y = @constraint(model, kon, [x^2, x] in SecondOrderCone())
    err(name) = ErrorException("Multiple constraints have the name $name.")
    @test_throws err("kon") constraint_by_name(model, "kon")
    set_name(kon, "con")
    _test_constraint_name_util(
        con,
        "kon",
        GenericQuadExpr{Float64,VariableRefType},
        MOI.EqualTo{Float64},
    )
    _test_constraint_name_util(
        kon,
        "con",
        Vector{GenericQuadExpr{Float64,VariableRefType}},
        MOI.SecondOrderCone,
    )
    set_name(con, "con")
    @test_throws err("con") constraint_by_name(model, "con")
    @test constraint_by_name(model, "kon") isa Nothing
    set_name(kon, "kon")
    _test_constraint_name_util(
        con,
        "con",
        GenericQuadExpr{Float64,VariableRefType},
        MOI.EqualTo{Float64},
    )
    _test_constraint_name_util(
        kon,
        "kon",
        Vector{GenericQuadExpr{Float64,VariableRefType}},
        MOI.SecondOrderCone,
    )
    return
end

function test_extension_PSD_constraint_errors(
    ModelType = Model,
    VariableRefType = VariableRef,
)
    model = ModelType()
    @variable(model, X[1:2, 1:2])
    err = ErrorException(
        "In `@constraint(model, X in MOI.PositiveSemidefiniteConeSquare(2))`:" *
        " instead of `MathOptInterface.PositiveSemidefiniteConeSquare(2)`," *
        " use `JuMP.PSDCone()`.",
    )
    @test_throws_strip(
        err,
        @constraint(model, X in MOI.PositiveSemidefiniteConeSquare(2))
    )
    err = ErrorException(
        "In `@constraint(model, X in MOI.PositiveSemidefiniteConeTriangle(2))`:" *
        " instead of `MathOptInterface.PositiveSemidefiniteConeTriangle(2)`," *
        " use `JuMP.PSDCone()`.",
    )
    @test_throws_strip(
        err,
        @constraint(model, X in MOI.PositiveSemidefiniteConeTriangle(2))
    )
    return
end

function test_extension_matrix_constraint_errors(
    ModelType = Model,
    VariableRefType = VariableRef,
)
    model = ModelType()
    @variable(model, X[1:2, 1:2])
    err = ErrorException(
        "In `@constraint(model, X in MOI.SecondOrderCone(4))`: unexpected " *
        "matrix in vector constraint. Do you need to flatten the matrix " *
        "into a vector using `vec()`?",
    )
    # Note: this should apply to any MOI.AbstractVectorSet. We just pick
    # SecondOrderCone for convenience.
    @test_throws_strip(err, @constraint(model, X in MOI.SecondOrderCone(4)))
    return
end

function test_extension_nonsensical_SDP_constraint(
    ModelType = Model,
    VariableRefType = VariableRef,
)
    T = value_type(ModelType)
    m = ModelType()
    @test_throws_strip(
        ErrorException(
            "In `@variable(m, unequal[1:5, 1:6], PSD)`: Symmetric variables must be square. Got size (5, 6).",
        ),
        @variable(m, unequal[1:5, 1:6], PSD)
    )
    # Some of these errors happen at compile time, so we can't use @test_throws
    @test_throws MethodError @variable(m, notone[1:5, 2:6], PSD)
    @test_throws MethodError @variable(m, oneD[1:5], PSD)
    @test_throws MethodError @variable(m, threeD[1:5, 1:5, 1:5], PSD)
    Y = T[1 2; 21//10 3]
    function _ErrorException(m)
        return ErrorException(
            "In `$m`: Non-symmetric bounds, integrality or starting values " *
            "for symmetric variable.",
        )
    end
    @test_throws_strip(
        _ErrorException("@variable(m, foo[i = 1:2, j = 1:2] >= Y[i, j], PSD)"),
        @variable(m, foo[i = 1:2, j = 1:2] >= Y[i, j], PSD),
    )
    @test_throws_strip(
        _ErrorException("@variable(m, foo[i = 1:2, j = 1:2] <= Y[i, j], PSD)"),
        @variable(m, foo[i = 1:2, j = 1:2] <= Y[i, j], PSD),
    )
    @test_throws_strip(
        _ErrorException(
            "@variable(m, foo[i = 1:2, j = 1:2] >= Y[i, j], Symmetric)",
        ),
        @variable(m, foo[i = 1:2, j = 1:2] >= Y[i, j], Symmetric),
    )
    @test_throws_strip(
        _ErrorException(
            "@variable(m, foo[i = 1:2, j = 1:2] <= Y[i, j], Symmetric)",
        ),
        @variable(m, foo[i = 1:2, j = 1:2] <= Y[i, j], Symmetric),
    )
    return
end

function test_extension_sum_constraint(
    ModelType = Model,
    VariableRefType = VariableRef,
)
    model = ModelType()
    @variable(model, x[1:3, 1:3])
    @variable(model, y)
    C = [1 2 3; 4 5 6; 7 8 9]
    @test_expression sum(C[i, j] * x[i, j] for i in 1:2, j in 2:3)
    @test_expression sum(C[i, j] * x[i, j] for i in 1:3, j in 1:3 if i != j) - y
    @test isequal_canonical(
        @expression(model, sum(C[i, j] * x[i, j] for i in 1:3, j in 1:i)),
        sum(C[i, j] * x[i, j] for i in 1:3 for j in 1:i),
    )
    @test_expression sum(C[i, j] * x[i, j] for i in 1:3 for j in 1:i)
    @test_expression sum(C[i, j] * x[i, j] for i in 1:3 if true for j in 1:i)
    @test_expression sum(
        C[i, j] * x[i, j] for i in 1:3 if true for j in 1:i if true
    )
    @test_expression sum(0 * x[i, 1] for i in 1:3)
    @test_expression sum(0 * x[i, 1] + y for i in 1:3)
    @test_expression sum(0 * x[i, 1] + y for i in 1:3 for j in 1:3)
    return
end

function test_all_constraints_scalar()
    model = Model()
    @variable(model, x >= 0)
    @test 1 == @inferred num_constraints(
        model,
        VariableRef,
        MOI.GreaterThan{Float64},
    )
    ref =
        @inferred all_constraints(model, VariableRef, MOI.GreaterThan{Float64})
    @test ref == [LowerBoundRef(x)]
    @test 0 ==
          @inferred num_constraints(model, AffExpr, MOI.GreaterThan{Float64})
    aff_constraints = all_constraints(model, AffExpr, MOI.GreaterThan{Float64})
    @test isempty(aff_constraints)
    err = ErrorException(
        "`MathOptInterface.GreaterThan` is not a " *
        "concrete type. Did you miss a type parameter?",
    )
    @test_throws err num_constraints(model, AffExpr, MOI.GreaterThan)
    @test_throws err all_constraints(model, AffExpr, MOI.GreaterThan)
    err = ErrorException(
        "`JuMP.GenericAffExpr` is not a concrete type. Did you miss a type parameter?",
    )
    @test_throws err try
        num_constraints(model, GenericAffExpr, MOI.ZeroOne)
    catch e
        error(replace(e.msg, "" => ""))
    end
    @test_throws err try
        all_constraints(model, GenericAffExpr, MOI.ZeroOne)
    catch e
        error(replace(e.msg, "" => ""))
    end
    return
end

function test_all_constraints_vector()
    model = Model()
    @variable(model, x[1:2, 1:2], Symmetric)
    csdp = @constraint(model, x in PSDCone())
    csoc = @constraint(model, [x[1], 1] in SecondOrderCone())
    csos = @constraint(model, [x[2]^2, 1] in MOI.SOS1([1.0, 2.0]))
    @test_throws(
        DimensionMismatch,
        @constraint(model, [x[2]^2, 1] in MOI.SOS1([1.0, 2.0, 3.0]))
    )
    @test_throws(
        DimensionMismatch,
        @constraint(model, [x[2]^2, 1] in MOI.SOS2([1.0, 2.0, 3.0]))
    )
    @test 1 == @inferred num_constraints(
        model,
        Vector{VariableRef},
        MOI.PositiveSemidefiniteConeTriangle,
    )
    ref = all_constraints(
        model,
        Vector{VariableRef},
        MOI.PositiveSemidefiniteConeTriangle,
    )
    @test ref == [csdp]
    @test 1 ==
          @inferred num_constraints(model, Vector{AffExpr}, MOI.SecondOrderCone)
    ref = all_constraints(model, Vector{AffExpr}, MOI.SecondOrderCone)
    @test ref == [csoc]
    @test 1 ==
          @inferred num_constraints(model, Vector{QuadExpr}, MOI.SOS1{Float64})
    ref = all_constraints(model, Vector{QuadExpr}, MOI.SOS1{Float64})
    @test ref == [csos]
    @test 0 == @inferred num_constraints(
        model,
        Vector{AffExpr},
        MOI.PositiveSemidefiniteConeTriangle,
    )
    aff_constraints = all_constraints(
        model,
        Vector{AffExpr},
        MOI.PositiveSemidefiniteConeTriangle,
    )
    @test isempty(aff_constraints)
    err = ErrorException(
        replace(
            "`$(GenericAffExpr{Float64})` is not a concrete type. Did you " *
            "miss a type parameter?",
            "" => "",
        ),
    )
    @test_throws err try
        num_constraints(
            model,
            Vector{GenericAffExpr{Float64}},
            MOI.PositiveSemidefiniteConeTriangle,
        )
    catch e
        error(replace(e.msg, "" => ""))
    end
    @test_throws err try
        all_constraints(
            model,
            Vector{GenericAffExpr{Float64}},
            MOI.SecondOrderCone,
        )
    catch e
        error(replace(e.msg, "" => ""))
    end
    err = ErrorException(
        "`MathOptInterface.SOS1` is not a " *
        "concrete type. Did you miss a type parameter?",
    )
    @test_throws err all_constraints(
        model,
        Vector{GenericQuadExpr{Float64,VariableRef}},
        MOI.SOS1,
    )
    return
end

function test_list_of_constraint_types()
    model = Model()
    @variable(model, x >= 0, Bin)
    @constraint(model, 2x <= 1)
    @constraint(model, [x, x] in SecondOrderCone())
    @constraint(model, [2x 1; 1 x] in PSDCone())
    @constraint(model, [x^2, x] in RotatedSecondOrderCone())
    constraint_types = @inferred list_of_constraint_types(model)
    @test Set(constraint_types) == Set([
        (VariableRef, MOI.ZeroOne),
        (VariableRef, MOI.GreaterThan{Float64}),
        (AffExpr, MOI.LessThan{Float64}),
        (Vector{VariableRef}, MOI.SecondOrderCone),
        (Vector{AffExpr}, MOI.PositiveSemidefiniteConeSquare),
        (Vector{QuadExpr}, MOI.RotatedSecondOrderCone),
    ])
    return
end

function test_dual_start()
    model = Model()
    @variable(model, x)
    con = @constraint(model, 2x <= 1)
    @test dual_start_value(con) === nothing
    set_dual_start_value(con, 2)
    @test dual_start_value(con) == 2.0
    set_dual_start_value(con, nothing)
    @test dual_start_value(con) === nothing
    return
end

function test_dual_start_vector()
    model = Model()
    @variable(model, x)
    con_vec = @constraint(model, [x, x] in SecondOrderCone())
    @test dual_start_value(con_vec) === nothing
    set_dual_start_value(con_vec, [1.0, 3.0])
    @test dual_start_value(con_vec) == [1.0, 3.0]
    set_dual_start_value(con_vec, nothing)
    @test dual_start_value(con_vec) === nothing
    return
end

function test_primal_start()
    model = Model()
    @variable(model, x)
    con = @constraint(model, 2x <= 1)
    @test start_value(con) === nothing
    set_start_value(con, 2)
    @test start_value(con) == 2.0
    set_start_value(con, nothing)
    @test start_value(con) === nothing
    return
end

function test_primal_start_vector()
    model = Model()
    @variable(model, x)
    con_vec = @constraint(model, [x, x] in SecondOrderCone())
    @test start_value(con_vec) === nothing
    set_start_value(con_vec, [1.0, 3.0])
    @test start_value(con_vec) == [1.0, 3.0]
    set_start_value(con_vec, nothing)
    @test start_value(con_vec) === nothing
    return
end

function test_change_coefficient()
    model = Model()
    x = @variable(model)
    con_ref = @constraint(model, 2 * x == -1)
    @test normalized_coefficient(con_ref, x) == 2.0
    set_normalized_coefficient(con_ref, x, 1.0)
    @test normalized_coefficient(con_ref, x) == 1.0
    set_normalized_coefficient(con_ref, x, 3)  # Check type promotion.
    @test normalized_coefficient(con_ref, x) == 3.0
    quad_con = @constraint(model, x^2 == 0)
    @test normalized_coefficient(quad_con, x) == 0.0
    set_normalized_coefficient(quad_con, x, 2)
    @test normalized_coefficient(quad_con, x) == 2.0
    @test isequal_canonical(constraint_object(quad_con).func, x^2 + 2x)
    return
end

function test_change_coefficients()
    model = Model()
    @variable(model, x)
    @constraint(model, con, [2x + 3x, 4x] in MOI.Nonnegatives(2))
    @test isequal_canonical(constraint_object(con).func, [5.0x, 4.0x])
    set_normalized_coefficients(con, x, [(Int64(1), 3.0)])
    @test isequal_canonical(constraint_object(con).func, [3.0x, 4.0x])
    set_normalized_coefficients(con, x, [(Int64(1), 2.0), (Int64(2), 5.0)])
    @test isequal_canonical(constraint_object(con).func, [2.0x, 5.0x])
    return
end

function test_change_rhs()
    model = Model()
    x = @variable(model)
    con_ref = @constraint(model, 2 * x <= 1)
    @test normalized_rhs(con_ref) == 1.0
    set_normalized_rhs(con_ref, 2.0)
    @test normalized_rhs(con_ref) == 2.0
    con_ref = @constraint(model, 2 * x - 1 == 1)
    @test normalized_rhs(con_ref) == 2.0
    set_normalized_rhs(con_ref, 3)
    @test normalized_rhs(con_ref) == 3.0
    con_ref = @constraint(model, 0 <= 2 * x <= 1)
    @test_throws MethodError set_normalized_rhs(con_ref, 3)
    return
end

function test_add_to_function_constant_scalar()
    model = Model()
    x = @variable(model)
    con_ref = @constraint(model, 2 <= 2 * x <= 3)
    con = constraint_object(con_ref)
    @test isequal_canonical(jump_function(con), 2x)
    @test moi_set(con) == MOI.Interval(2.0, 3.0)
    add_to_function_constant(con_ref, 1.0)
    con = constraint_object(con_ref)
    @test isequal_canonical(jump_function(con), 2x)
    @test moi_set(con) == MOI.Interval(1.0, 2.0)
    return
end

function test_add_to_function_constant_vector()
    model = Model()
    x = @variable(model)
    con_ref = @constraint(model, [x + 1, x - 1] in MOI.Nonnegatives(2))
    con = constraint_object(con_ref)
    @test isequal_canonical(jump_function(con), [x + 1, x - 1])
    @test moi_set(con) == MOI.Nonnegatives(2)
    add_to_function_constant(con_ref, [2, 3])
    con = constraint_object(con_ref)
    @test isequal_canonical(jump_function(con), [x + 3, x + 2])
    @test moi_set(con) == MOI.Nonnegatives(2)
    return
end

function test_value_constraint_var(
    ModelType = Model,
    VariableRefType = VariableRef,
)
    model = ModelType()
    @variable(model, x[1:2])
    @constraint(model, c1, x[1] + x[2] <= 3.0)
    @constraint(model, c2, x[1]^2 + x[2]^2 <= 3.0)
    @constraint(model, c3, [1.0, x[1], x[2]] in SecondOrderCone())
    vals = Dict(x[1] => 1.0, x[2] => 2.0)
    f = vidx -> vals[vidx]
    @test value(f, c1) === 3.0 # Affine expression
    @test value(f, c2) === 5.0 # Quadratic expression
    @test value(f, c3) == [1.0, 1.0, 2.0] # Vector expression
    return
end

function _test_shadow_price_util(
    model_string,
    constraint_dual,
    constraint_shadow,
)
    model = Model()
    MOIU.loadfromstring!(backend(model), model_string)
    set_optimizer(
        model,
        () -> MOIU.MockOptimizer(
            MOIU.Model{Float64}();
            eval_objective_value = false,
            eval_variable_constraint_dual = false,
        ),
    )
    optimize!(model)
    mock_optimizer = unsafe_backend(model)
    MOI.set(mock_optimizer, MOI.TerminationStatus(), MOI.OPTIMAL)
    MOI.set(mock_optimizer, MOI.DualStatus(), MOI.FEASIBLE_POINT)
    optimize!(model)
    for (key, val) in constraint_dual
        ci = if key isa String
            MOI.get(backend(model), MOI.ConstraintIndex, key)
        else
            x = MOI.get(backend(model), MOI.VariableIndex, key[1])
            MOI.ConstraintIndex{MOI.VariableIndex,key[2]}(x.value)
        end
        constraint_ref = ConstraintRef(model, ci, ScalarShape())
        MOI.set(
            mock_optimizer,
            MOI.ConstraintDual(),
            optimizer_index(constraint_ref),
            val,
        )
        @test dual(constraint_ref) == val
        @test shadow_price(constraint_ref) == constraint_shadow[key]
    end
    return
end

function test_shadow_price()
    _test_shadow_price_util(
        """
        variables: x, y
        minobjective: -1.0*x
        x <= 2.0
        y >= 0.0
        c: x + y <= 1.0
        """,
        Dict(
            ("x", MOI.LessThan{Float64}) => 0.0,
            ("y", MOI.GreaterThan{Float64}) => 1.0,
            "c" => -1.0,
        ),
        Dict(
            ("x", MOI.LessThan{Float64}) => 0.0,
            ("y", MOI.GreaterThan{Float64}) => -1.0,
            "c" => -1.0,
        ),
    )
    _test_shadow_price_util(
        """
        variables: x, y
        maxobjective: 1.0*x
        x <= 2.0
        y >= 0.0
        c: x + y <= 1.0
        """,
        Dict(
            ("x", MOI.LessThan{Float64}) => 0.0,
            ("y", MOI.GreaterThan{Float64}) => 1.0,
            "c" => -1.0,
        ),
        Dict(
            ("x", MOI.LessThan{Float64}) => 0.0,
            ("y", MOI.GreaterThan{Float64}) => 1.0,
            "c" => 1.0,
        ),
    )
    _test_shadow_price_util(
        """
        variables: x, y
        maxobjective: 1.0*x
        x <= 2.0
        y >= 0.0
        """,
        Dict(
            ("x", MOI.LessThan{Float64}) => -1.0,
            ("y", MOI.GreaterThan{Float64}) => 0.0,
        ),
        Dict(
            ("x", MOI.LessThan{Float64}) => 1.0,
            ("y", MOI.GreaterThan{Float64}) => 0.0,
        ),
    )
    _test_shadow_price_util(
        """
        variables: x
        maxobjective: 1.0*x
        x == 2.0
        """,
        Dict(("x", MOI.EqualTo{Float64}) => -1.0),
        Dict(("x", MOI.EqualTo{Float64}) => 1.0),
    )
    _test_shadow_price_util(
        """
        variables: x
        minobjective: 1.0*x
        x == 2.0
        """,
        Dict(("x", MOI.EqualTo{Float64}) => 1.0),
        Dict(("x", MOI.EqualTo{Float64}) => -1.0),
    )
    return
end

function test_extension_abstractarray_vector_constraint(
    ModelType = Model,
    VariableRefType = VariableRef,
)
    T = value_type(ModelType)
    model = ModelType()
    @variable(model, x[1:2, 1:2])
    c = @constraint(model, view(x, 1:4) in SOS1())
    obj = constraint_object(c)
    @test obj.func == x[1:4]
    @test obj.set == MOI.SOS1(T[1, 2, 3, 4])
    return
end

function test_extension_constraint_inference(
    ModelType = Model,
    VariableRefType = VariableRef,
)
    T = value_type(ModelType)
    model = ModelType()
    @variable(model, x)
    foo(model, x) = @constraint(model, 2x <= 1)
    c = @inferred foo(model, x)
    obj = constraint_object(c)
    @test obj.func == 2x
    @test obj.set == MOI.LessThan(one(T))
    return
end

struct _UnsupportedConstraintName <: MOI.AbstractOptimizer end

MOI.add_variable(::_UnsupportedConstraintName) = MOI.VariableIndex(1)

function MOI.supports_constraint(
    ::_UnsupportedConstraintName,
    ::Type{MOI.VectorOfVariables},
    ::Type{MOI.SOS1{Float64}},
)
    return true
end

function MOI.add_constraint(
    ::_UnsupportedConstraintName,
    ::MOI.VectorOfVariables,
    ::MOI.SOS1{Float64},
)
    return MOI.ConstraintIndex{MOI.VectorOfVariables,MOI.SOS1{Float64}}(1)
end

MOI.is_empty(::_UnsupportedConstraintName) = true

function test_unsupported_ConstraintName()
    model = direct_model(_UnsupportedConstraintName())
    @variable(model, x)
    @constraint(model, c, [x, x] in SOS1())
    @test c isa ConstraintRef
    @test name(c) == ""
    return
end

function test_PSDCone_constraints()
    model = Model()
    @variable(model, X[1:2, 1:2])
    Y = [1 1; 1 1]
    f = reshape(X .- Y, 4)
    c1 = @constraint(model, X >= Y, PSDCone())
    obj = constraint_object(c1)
    @test obj.func == f
    @test obj.set == MOI.PositiveSemidefiniteConeSquare(2)
    c2 = @constraint(model, Y <= X, PSDCone())
    @test constraint_object(c2).func == f
    c3 = @constraint(model, X - Y >= 0, PSDCone())
    @test constraint_object(c3).func == f
    c4 = @constraint(model, Y - X <= 0, PSDCone())
    @test constraint_object(c4).func == f
    @test_throws ErrorException @constraint(model, X >= 1, PSDCone())
    return
end

function test_PSDCone_Symmetric_constraints()
    model = Model()
    @variable(model, X[1:2, 1:2], Symmetric)
    Y = LinearAlgebra.Symmetric([1 1; 1 1])
    f = reshape(X - Y, 4)
    # Julia 1.0 doesn't maintain symmetry for X - Y, so only test this on
    # Julia 1.6 and higher.
    c1 = @constraint(model, X >= Y, PSDCone())
    obj = constraint_object(c1)
    @test obj.func == f[[1, 3, 4]]
    @test obj.set == MOI.PositiveSemidefiniteConeTriangle(2)
    c2 = @constraint(model, Y <= X, PSDCone())
    @test constraint_object(c2).func == f[[1, 3, 4]]
    c3 = @constraint(model, LinearAlgebra.Symmetric(X - Y) >= 0, PSDCone())
    @test constraint_object(c3).func == f[[1, 3, 4]]
    c4 = @constraint(model, LinearAlgebra.Symmetric(Y - X) <= 0, PSDCone())
    @test constraint_object(c4).func == f[[1, 3, 4]]
    c5 = @constraint(model, X >= 0, PSDCone())
    @test constraint_object(c5).func == 1.0 .* X[[1, 3, 4]]
    return
end

"""
    test_set_inequalities()

Test the syntax `@constraint(model, X >= Y, Set())` is the same as
`@constraint(model, X - Y in Set())`.
"""
function test_set_inequalities()
    model = Model()
    @variable(model, X[1:2])
    Y = [3.0, 4.0]
    f1 = 1 .* X
    f2 = X .- Y
    f3 = X .- 1
    c1 = @constraint(model, X >= 0, MOI.Nonnegatives(2))
    c2 = @constraint(model, X >= Y, MOI.Nonnegatives(2))
    c3 = @constraint(model, 0 <= X, MOI.Nonnegatives(2))
    c4 = @constraint(model, Y <= X, MOI.Nonnegatives(2))
    c5 = @constraint(model, 0 <= X .- 1, MOI.Nonnegatives(2))
    c6 = @constraint(model, X .- 1 >= 0, MOI.Nonnegatives(2))
    @test constraint_object(c1).func == f1
    @test constraint_object(c2).func == f2
    @test constraint_object(c3).func == f1
    @test constraint_object(c4).func == f2
    @test constraint_object(c5).func == f3
    @test constraint_object(c6).func == f3
    @test constraint_object(c1).set == MOI.Nonnegatives(2)
    @test constraint_object(c2).set == MOI.Nonnegatives(2)
    @test constraint_object(c3).set == MOI.Nonnegatives(2)
    @test constraint_object(c4).set == MOI.Nonnegatives(2)
    @test constraint_object(c5).set == MOI.Nonnegatives(2)
    @test constraint_object(c6).set == MOI.Nonnegatives(2)
    @test_throws ErrorException @constraint(model, X >= 1, MOI.Nonnegatives(2))
    @test_throws ErrorException @constraint(model, 1 <= X, MOI.Nonnegatives(2))
    return
end

function test_num_constraints()
    model = Model()
    @variable(model, x >= 0, Int)
    @constraint(model, 2x <= 1)
    @test num_constraints(model; count_variable_in_set_constraints = true) == 3
    @test num_constraints(model; count_variable_in_set_constraints = false) == 1
    return
end

function test_num_constraints_UndefKeywordError()
    model = Model()
    @variable(model, x >= 0, Int)
    @constraint(model, 2x <= 1)
    @test_throws UndefKeywordError num_constraints(model)
    return
end

function test_num_constraints_nonlinear()
    model = Model()
    @variable(model, x >= 0, Int)
    @constraint(model, 2x <= 1)
    @NLconstraint(model, sqrt(x) <= 3)
    @test num_constraints(model; count_variable_in_set_constraints = true) == 4
    @test num_constraints(model; count_variable_in_set_constraints = false) == 2
    return
end

function test_all_constraints()
    model = Model()
    @variable(model, x >= 0, Int)
    c = @constraint(model, 2x <= 1)
    nl_c = @NLconstraint(model, x^2 <= 1)
    ret = all_constraints(model; include_variable_in_set_constraints = true)
    @test length(ret) == 4
    @test c in ret
    @test nl_c in ret
    @test IntegerRef(x) in ret
    @test LowerBoundRef(x) in ret
    ret = all_constraints(model; include_variable_in_set_constraints = false)
    @test length(ret) == 2
    @test c in ret
    @test nl_c in ret
    return
end

function test_relax_with_penalty!_default()
    model = Model()
    @variable(model, x >= 0)
    map = relax_with_penalty!(model)
    @test isempty(map)
    @constraint(model, c1, x <= 1)
    @constraint(model, c2, x == 0)
    map = relax_with_penalty!(model)
    @test length(map) == 2
    @test map[c1] isa AffExpr
    @test map[c2] isa AffExpr
    @test num_variables(model) == 4
    @test objective_sense(model) == MOI.MIN_SENSE
    @test objective_function(model) == map[c1] + map[c2]
    return
end

function test_relax_with_penalty!_max()
    model = Model()
    @variable(model, x >= 0)
    @constraint(model, c1, x <= 1)
    @constraint(model, c2, x == 0)
    @objective(model, Max, 1.0 * x + 2.5)
    map = relax_with_penalty!(model)
    @test length(map) == 2
    @test map[c1] isa AffExpr
    @test map[c2] isa AffExpr
    @test num_variables(model) == 4
    @test objective_sense(model) == MOI.MAX_SENSE
    @test objective_function(model) == x + 2.5 - map[c1] - map[c2]
    return
end

function test_relax_with_penalty!_constant()
    model = Model()
    @variable(model, x >= 0)
    map = relax_with_penalty!(model)
    @test isempty(map)
    @constraint(model, c1, x <= 1)
    @constraint(model, c2, x == 0)
    map = relax_with_penalty!(model; default = 2)
    @test length(map) == 2
    @test map[c1] isa AffExpr
    @test map[c2] isa AffExpr
    @test num_variables(model) == 4
    @test objective_sense(model) == MOI.MIN_SENSE
    @test objective_function(model) == 2.0 * map[c1] + 2.0 * map[c2]
    return
end

function test_relax_with_penalty!_specific()
    model = Model()
    @variable(model, x >= 0)
    map = relax_with_penalty!(model)
    @test isempty(map)
    @constraint(model, c1, x <= 1)
    @constraint(model, c2, x == 0)
    map = relax_with_penalty!(model, Dict(c1 => 3.0))
    @test length(map) == 1
    @test map[c1] isa AffExpr
    @test num_variables(model) == 2
    @test objective_sense(model) == MOI.MIN_SENSE
    @test objective_function(model) == 3.0 * map[c1]
    return
end

function test_relax_with_penalty!_specific_with_default()
    model = Model()
    @variable(model, x >= 0)
    map = relax_with_penalty!(model)
    @test isempty(map)
    @constraint(model, c1, x <= 1)
    @constraint(model, c2, x == 0)
    map = relax_with_penalty!(model, Dict(c1 => 3.0); default = 1)
    @test length(map) == 2
    @test map[c1] isa AffExpr
    @test num_variables(model) == 4
    @test objective_sense(model) == MOI.MIN_SENSE
    @test objective_function(model) == 3 * map[c1] + map[c2]
    return
end

function test_Hermitian_PSD_constraint()
    model = Model()
    set_optimizer(
        model,
        () -> MOIU.MockOptimizer(
            MOIU.Model{Float64}();
            eval_objective_value = false,
            eval_variable_constraint_dual = false,
        ),
    )
    @variable(model, x)
    @variable(model, y)
    @variable(model, w)
    A = [x 1im; -1im -y] - [1 (x+w*im); (x-w*im) -2]
    @constraint(model, href, LinearAlgebra.Hermitian(A) in HermitianPSDCone())
    _test_constraint_name_util(
        href,
        "href",
        Vector{GenericAffExpr{Float64,VariableRef}},
        MOI.HermitianPositiveSemidefiniteConeTriangle,
    )
    c = constraint_object(href)
    @test isequal_canonical(c.func[1], x - 1)
    @test isequal_canonical(c.func[2], -x)
    @test isequal_canonical(c.func[3], 2 - y)
    @test isequal_canonical(c.func[4], 1 - w)
    @test c.set == MOI.HermitianPositiveSemidefiniteConeTriangle(2)
    @test c.shape isa HermitianMatrixShape
    MOIU.attach_optimizer(model)
    model.is_model_dirty = false
    mock_optimizer = unsafe_backend(model)
    MOI.set(mock_optimizer, MOI.TerminationStatus(), MOI.OPTIMAL)
    MOI.set(mock_optimizer, MOI.DualStatus(), MOI.FEASIBLE_POINT)
    F = MOI.VectorAffineFunction{Float64}
    for (i, v) in enumerate(all_variables(model))
        MOI.set(mock_optimizer, MOI.VariablePrimal(), optimizer_index(v), i)
    end
    H = value(href)
    @test H isa LinearAlgebra.Hermitian
    @test parent(H) == [0 -1-2im; -1+2im 0]
    MOI.set(
        mock_optimizer,
        MOI.ConstraintDual(),
        optimizer_index(href),
        1:MOI.dimension(c.set),
    )
    H = dual(href)
    @test H isa LinearAlgebra.Hermitian
    @test parent(H) == [1 2+4im; 2-4im 3]
    return
end

function test_extension_HermitianPSDCone_errors(
    ModelType = Model,
    VariableRefType = VariableRef,
)
    AffExprType =
        JuMP.GenericAffExpr{Complex{value_type(ModelType)},VariableRefType}
    model = ModelType()
    T = value_type(ModelType)
    @variable(model, x)
    @variable(model, y)
    aff_str = "$AffExprType"
    err = ErrorException(
        "In `@constraint(model, H in HermitianPSDCone(), unknown_kw = 1)`:" *
        " Unrecognized constraint building format. Tried to invoke " *
        "`build_constraint(error, $aff_str[" *
        "x im; -im -y], $(HermitianPSDCone()); unknown_kw = 1)`, but no " *
        "such method exists. This is due to specifying an unrecognized " *
        "function, constraint set, and/or extra positional/keyword " *
        "arguments.\n\nIf you're trying to create a JuMP extension, you " *
        "need to implement `build_constraint` to accomodate these arguments.",
    )
    H = LinearAlgebra.Hermitian([x 1im; -1im -y])
    @test_throws_strip(
        err,
        @constraint(model, H in HermitianPSDCone(), unknown_kw = 1),
    )
    return
end

function test_extension_Nonnegatives(
    ModelType = Model,
    VariableRefType = VariableRef,
)
    model = ModelType()
    @variable(model, x[1:2])
    @constraint(model, c1, x in Nonnegatives())
    con = constraint_object(c1)
    @test con.func == x
    @test con.set == MOI.Nonnegatives(2)
    A = [1 2; 3 4]
    b = [5, 6]
    @constraint(model, c2, A * x >= b)
    con = constraint_object(c2)
    @test con.func == A * x - b
    @test con.set == MOI.Nonnegatives(2)
    return
end

function test_extension_Nonpositives(
    ModelType = Model,
    VariableRefType = VariableRef,
)
    model = ModelType()
    @variable(model, x[1:2])
    @constraint(model, c1, x in Nonpositives())
    con = constraint_object(c1)
    @test con.func == x
    @test con.set == MOI.Nonpositives(2)
    A = [1 2; 3 4]
    b = [5, 6]
    @constraint(model, c2, A * x <= b)
    con = constraint_object(c2)
    @test con.func == A * x - b
    @test con.set == MOI.Nonpositives(2)
    return
end

function test_extension_Zeros(ModelType = Model, VariableRefType = VariableRef)
    model = ModelType()
    @variable(model, x[1:2])
    @constraint(model, c1, x in Zeros())
    con = constraint_object(c1)
    @test con.func == x
    @test con.set == MOI.Zeros(2)
    A = [1 2; 3 4]
    b = [5, 6]
    @constraint(model, c2, A * x == b)
    con = constraint_object(c2)
    @test con.func == A * x - b
    @test con.set == MOI.Zeros(2)
    return
end

function test_hermitian_in_zeros()
    model = Model()
    @variable(model, H[1:2, 1:2] in HermitianPSDCone())
    c = @constraint(model, H == LinearAlgebra.I)
    con = constraint_object(c)
    @test isequal_canonical(con.func[1], real(H[1, 1]) - 1)
    @test isequal_canonical(con.func[2], real(H[1, 2]))
    @test isequal_canonical(con.func[3], real(H[2, 2]) - 1)
    @test isequal_canonical(con.func[4], imag(H[1, 2]))
    @test con.set == MOI.Zeros(4)
    @test occursin("Zeros()", sprint(show, c))
    @test value(x -> 1.0, c) isa LinearAlgebra.Hermitian
    return
end

function test_symmetric_in_zeros()
    model = Model()
    @variable(model, H[1:2, 1:2], Symmetric)
    c = @constraint(model, H == LinearAlgebra.I)
    con = constraint_object(c)
    @test isequal_canonical(con.func[1], H[1, 1] - 1)
    @test isequal_canonical(con.func[2], H[1, 2] - 0)
    @test isequal_canonical(con.func[3], H[2, 2] - 1)
    @test con.set == MOI.Zeros(3)
    @test occursin("Zeros()", sprint(show, c))
    @test value(x -> 1.0, c) isa LinearAlgebra.Symmetric
    return
end

function test_semicontinuous()
    model = Model()
    @variable(model, x in Semicontinuous(2, 3))
    @variable(model, y)
    c = @constraint(model, y + 1 in Semicontinuous(2.5, 3))
    c_obj = constraint_object(c)
    @test isequal_canonical(c_obj.func, y + 1)
    @test c_obj.set == MOI.Semicontinuous(2.5, 3.0)
    return
end

function test_semiinteger()
    model = Model()
    @variable(model, x in Semiinteger(2, 3))
    @variable(model, y)
    c = @constraint(model, y + 1 in Semiinteger(2.5, 3))
    c_obj = constraint_object(c)
    @test isequal_canonical(c_obj.func, y + 1)
    @test c_obj.set == MOI.Semiinteger(2.5, 3.0)
    return
end

function test_symmetric_vectorize_allocations()
    if VERSION < v"1.8"
        return
    end
    model = Model()
    @variable(model, x[1:2])
    C = SparseArrays.sparse([0 1; 0 0])
    X = LinearAlgebra.Symmetric(C - SparseArrays.spdiagm(x))
    @constraint(model, X in PSDCone())  # Once for compilation
    @test (@allocated @constraint(model, X in PSDCone())) <= 944
    return
end

function test_eltype_for_constraint_primal_float64()
    model = Model() do
        return MOI.Utilities.MockOptimizer(MOI.Utilities.Model{Float64}())
    end
    @variable(model, x >= 0)
    c1 = LowerBoundRef(x)
    @constraint(model, c2, 1.0 * x == 1.0)
    @constraint(model, c3, 1.0 * x^2 == 1.0)
    @constraint(model, c4, [x] in Zeros())
    @constraint(model, c5, [1.0 * x] == [1.0])
    @constraint(model, c6, [x^2] == [1.0])
    MOI.Utilities.attach_optimizer(model)
    mock = unsafe_backend(model)
    MOI.set(mock, MOI.TerminationStatus(), MOI.OPTIMAL)
    MOI.set(mock, MOI.DualStatus(), MOI.FEASIBLE_POINT)
    MOI.set(mock, MOI.VariablePrimal(), optimizer_index(x), 0.5)
    optimize!(model)
    @test JuMP._value_type(Model, Int) == Any  # For any user-defined functions.
    @test @inferred(value(c1)) == 0.5
    @test @inferred(value(c2)) == 0.5
    @test @inferred(value(c3)) == 0.25
    @test @inferred(value(c4)) == [0.5]
    @test @inferred(value(c5)) == [0.5 - 1]
    @test @inferred(value(c6)) == [0.25 - 1]
    return
end

function test_eltype_for_constraint_dual_float64()
    model = Model() do
        return MOI.Utilities.MockOptimizer(
            MOI.Utilities.Model{Float64}();
            eval_variable_constraint_dual = false,
        )
    end
    @variable(model, x >= 0)
    c1 = LowerBoundRef(x)
    @constraint(model, c2, 1.0 * x == 1.0)
    @constraint(model, c3, 1.0 * x^2 == 1.0)
    @constraint(model, c4, [x] in Zeros())
    @constraint(model, c5, [1.0 * x] == [1.0])
    @constraint(model, c6, [x^2] == [1.0])
    MOI.Utilities.attach_optimizer(model)
    mock = unsafe_backend(model)
    MOI.set(mock, MOI.TerminationStatus(), MOI.OPTIMAL)
    MOI.set(mock, MOI.DualStatus(), MOI.FEASIBLE_POINT)
    MOI.set(mock, MOI.ConstraintDual(), optimizer_index(c1), 0.5)
    MOI.set(mock, MOI.ConstraintDual(), optimizer_index(c2), 0.5)
    MOI.set(mock, MOI.ConstraintDual(), optimizer_index(c3), 0.5)
    MOI.set(mock, MOI.ConstraintDual(), optimizer_index(c4), [0.5])
    MOI.set(mock, MOI.ConstraintDual(), optimizer_index(c5), [0.5])
    MOI.set(mock, MOI.ConstraintDual(), optimizer_index(c6), [0.5])
    optimize!(model)
    @test @inferred(dual(c1)) == 0.5
    @test @inferred(dual(c2)) == 0.5
    @test @inferred(dual(c3)) == 0.5
    @test @inferred(dual(c4)) == [0.5]
    @test @inferred(dual(c5)) == [0.5]
    @test @inferred(dual(c6)) == [0.5]
    return
end

function test_eltype_for_constraint_primal_complex_float64()
    model = Model() do
        return MOI.Utilities.MockOptimizer(MOI.Utilities.Model{Float64}())
    end
    @variable(model, x in ComplexPlane())
    @constraint(model, c2, 1.0 * x == 1.0 + 2im)
    @constraint(model, c3, 1.0 * x^2 == 1.0 + 2im)
    @constraint(model, c4, [x] in Zeros())
    @constraint(model, c5, [1.0 * x] == [1.0 + 2im])
    @constraint(model, c6, [x^2] == [1.0 + 2im])
    MOI.Utilities.attach_optimizer(model)
    mock = unsafe_backend(model)
    MOI.set(mock, MOI.TerminationStatus(), MOI.OPTIMAL)
    MOI.set(mock, MOI.DualStatus(), MOI.FEASIBLE_POINT)
    x_re, x_im = first(keys(real(x).terms)), first(keys(imag(x).terms))
    MOI.set(mock, MOI.VariablePrimal(), optimizer_index(x_re), 0.5)
    MOI.set(mock, MOI.VariablePrimal(), optimizer_index(x_im), 0.25)
    optimize!(model)
    @test @inferred(value(c2)) == (0.5 + 0.25im)
    @test @inferred(value(c3)) == (0.5 + 0.25im)^2
    @test @inferred(value(c4)) == [(0.5 + 0.25im)]
    @test @inferred(value(c5)) == [(0.5 + 0.25im) - (1 + 2im)]
    @test @inferred(value(c6)) == [(0.5 + 0.25im)^2 - (1 + 2im)]
    return
end

end

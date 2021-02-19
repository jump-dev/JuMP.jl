module TestConstraint

using JuMP
using LinearAlgebra
using Test

include(joinpath(@__DIR__, "utilities.jl"))
include(joinpath(@__DIR__, "JuMPExtension.jl"))

function _test_constraint_name_util(constraint, name, F::Type, S::Type)
    @test name == @inferred JuMP.name(constraint)
    model = constraint.model
    @test constraint.index == JuMP.constraint_by_name(model, name).index
    if !(model isa JuMPExtension.MyModel)
        @test constraint.index ==
              JuMP.constraint_by_name(model, name, F, S).index
    end
end

function test_SingleVariable_constraints(ModelType, ::Any)
    m = ModelType()
    @variable(m, x)

    # x <= 10.0 doesn't translate to a SingleVariable constraint because
    # the LHS is first subtracted to form x - 10.0 <= 0.
    @constraint(m, cref, x in MOI.LessThan(10.0))
    _test_constraint_name_util(
        cref,
        "cref",
        JuMP.VariableRef,
        MOI.LessThan{Float64},
    )
    c = JuMP.constraint_object(cref)
    @test c.func == x
    @test c.set == MOI.LessThan(10.0)

    @variable(m, y[1:2])
    @constraint(m, cref2[i = 1:2], y[i] in MOI.LessThan(float(i)))
    _test_constraint_name_util(
        cref2[1],
        "cref2[1]",
        JuMP.VariableRef,
        MOI.LessThan{Float64},
    )
    c = JuMP.constraint_object(cref2[1])
    @test c.func == y[1]
    @test c.set == MOI.LessThan(1.0)
end

function test_VectorOfVariables_constraints(ModelType, ::Any)
    m = ModelType()
    @variable(m, x[1:2])

    cref = @constraint(m, x in MOI.Zeros(2))
    c = JuMP.constraint_object(cref)
    @test c.func == x
    @test c.set == MOI.Zeros(2)

    cref = @constraint(m, [x[2], x[1]] in MOI.Zeros(2))
    c = JuMP.constraint_object(cref)
    @test c.func == [x[2], x[1]]
    @test c.set == MOI.Zeros(2)

    @test_throws DimensionMismatch @constraint(m, [x[2], x[1]] in MOI.Zeros(3))
end

function test_AffExpr_scalar_constraints(ModelType, ::Any)
    model = ModelType()
    @variable(model, x)

    cref = @constraint(model, 2x <= 10)
    @test "" == @inferred JuMP.name(cref)
    JuMP.set_name(cref, "c")
    _test_constraint_name_util(cref, "c", JuMP.AffExpr, MOI.LessThan{Float64})

    c = JuMP.constraint_object(cref)
    @test JuMP.isequal_canonical(c.func, 2x)
    @test c.set == MOI.LessThan(10.0)

    cref = @constraint(model, 3x + 1 ≥ 10)
    c = JuMP.constraint_object(cref)
    @test JuMP.isequal_canonical(c.func, 3x)
    @test c.set == MOI.GreaterThan(9.0)

    cref = @constraint(model, 1 == -x)
    c = JuMP.constraint_object(cref)
    @test JuMP.isequal_canonical(c.func, 1.0x)
    @test c.set == MOI.EqualTo(-1.0)

    cref = @constraint(model, 2 == 1)
    c = JuMP.constraint_object(cref)
    @test JuMP.isequal_canonical(c.func, zero(JuMP.AffExpr))
    @test c.set == MOI.EqualTo(-1.0)
end

function test_AffExpr_vectorized_constraints(ModelType, ::Any)
    model = ModelType()
    @variable(model, x)
    err =
        ErrorException("In `@constraint(model, [x, 2x] == [1 - x, 3])`: Unexpected vector in scalar constraint. Did you mean to use the dot comparison operators like .==, .<=, and .>= instead?")
    @test_throws_strip err @constraint(model, [x, 2x] == [1 - x, 3])
    @test_macro_throws ErrorException begin
        @constraint(model, [x == 1 - x, 2x == 3])
    end
    cref = @constraint(model, [x, 2x] .== [1 - x, 3])
    c = JuMP.constraint_object.(cref)
    @test JuMP.isequal_canonical(c[1].func, 2.0x)
    @test c[1].set == MOI.EqualTo(1.0)
    @test JuMP.isequal_canonical(c[2].func, 2.0x)
    @test c[2].set == MOI.EqualTo(3.0)
end

function test_AffExpr_vector_constraints(ModelType, ::Any)
    model = ModelType()
    cref = @constraint(model, [1, 2] in MOI.Zeros(2))
    c = JuMP.constraint_object(cref)
    @test JuMP.isequal_canonical(c.func[1], zero(JuMP.AffExpr) + 1)
    @test JuMP.isequal_canonical(c.func[2], zero(JuMP.AffExpr) + 2)
    @test c.set == MOI.Zeros(2)
    @test c.shape isa JuMP.VectorShape
    @test_throws DimensionMismatch @constraint(model, [1, 2] in MOI.Zeros(3))
end

function test_delete_constraints(ModelType, ::Any)
    model = ModelType()
    @variable(model, x)
    constraint_ref = @constraint(model, 2x <= 1)
    @test JuMP.is_valid(model, constraint_ref)
    JuMP.delete(model, constraint_ref)
    @test !JuMP.is_valid(model, constraint_ref)
    second_model = ModelType()
    @test_throws Exception JuMP.delete(second_model, constraint_ref)
end

function test_batch_delete_constraints(ModelType, ::Any)
    model = ModelType()
    @variable(model, x[1:9])
    cons = [@constraint(model, sum(x[1:2:9]) <= 3)]
    push!(cons, @constraint(model, sum(x[2:2:8]) <= 2))
    push!(cons, @constraint(model, sum(x[1:3:9]) <= 1))
    @test all(JuMP.is_valid.(model, cons))
    JuMP.delete(model, cons[[1, 3]])
    @test all((!JuMP.is_valid).(model, cons[[1, 3]]))
    @test JuMP.is_valid(model, cons[2])
    second_model = ModelType()
    @test_throws Exception JuMP.delete(second_model, cons[[1, 3]])
    @test_throws Exception JuMP.delete(second_model, [cons[2]])
end

function test_two_sided_constraints(ModelType, ::Any)
    m = ModelType()
    @variable(m, x)
    @variable(m, y)

    @constraint(m, cref, 1.0 <= x + y + 1.0 <= 2.0)
    _test_constraint_name_util(
        cref,
        "cref",
        JuMP.AffExpr,
        MOI.Interval{Float64},
    )

    c = JuMP.constraint_object(cref)
    @test JuMP.isequal_canonical(c.func, x + y)
    @test c.set == MOI.Interval(0.0, 1.0)

    cref = @constraint(m, 2x - y + 2.0 ∈ MOI.Interval(-1.0, 1.0))
    @test "" == @inferred JuMP.name(cref)

    c = JuMP.constraint_object(cref)
    @test JuMP.isequal_canonical(c.func, 2x - y)
    @test c.set == MOI.Interval(-3.0, -1.0)
end

function test_broadcasted_constraint_eq(ModelType, ::Any)
    m = ModelType()
    @variable(m, x[1:2])

    A = [1.0 2.0; 3.0 4.0]
    b = [4.0, 5.0]

    cref = @constraint(m, A * x .== b)
    @test (2,) == @inferred size(cref)

    c1 = JuMP.constraint_object(cref[1])
    @test JuMP.isequal_canonical(c1.func, x[1] + 2x[2])
    @test c1.set == MOI.EqualTo(4.0)
    c2 = JuMP.constraint_object(cref[2])
    @test JuMP.isequal_canonical(c2.func, 3x[1] + 4x[2])
    @test c2.set == MOI.EqualTo(5.0)
end

function test_broadcasted_constraint_leq(ModelType, ::Any)
    m = ModelType()
    @variable(m, x[1:2, 1:2])

    UB = [1.0 2.0; 3.0 4.0]

    cref = @constraint(m, x .+ 1 .<= UB)
    @test (2, 2) == @inferred size(cref)
    for i in 1:2
        for j in 1:2
            c = JuMP.constraint_object(cref[i, j])
            @test JuMP.isequal_canonical(c.func, x[i, j] + 0)
            @test c.set == MOI.LessThan(UB[i, j] - 1)
        end
    end
end

function test_broadcasted_two_sided_constraint(ModelType, ::Any)
    m = ModelType()
    @variable(m, x[1:2])
    @variable(m, y[1:2])
    l = [1.0, 2.0]
    u = [3.0, 4.0]

    cref = @constraint(m, l .<= x + y .+ 1 .<= u)
    @test (2,) == @inferred size(cref)

    for i in 1:2
        c = JuMP.constraint_object(cref[i])
        @test JuMP.isequal_canonical(c.func, x[i] + y[i])
        @test c.set == MOI.Interval(l[i] - 1, u[i] - 1)
    end
end

function test_broadcasted_constraint_with_indices(ModelType, ::Any)
    m = ModelType()
    @variable m x[1:2]
    @constraint m cref1[i = 2:4] x .== [i, i + 1]
    ConstraintRefType = eltype(cref1[2])
    @test cref1 isa JuMP.Containers.DenseAxisArray{Vector{ConstraintRefType}}
    @constraint m cref2[i = 1:3, j = 1:4] x .≤ [i + j, i - j]
    ConstraintRefType = eltype(cref2[1])
    @test cref2 isa Matrix{Vector{ConstraintRefType}}
    @variable m y[1:2, 1:2]
    @constraint m cref3[i = 1:2] y[i, :] .== 1
    ConstraintRefType = eltype(cref3[1])
    @test cref3 isa Vector{Vector{ConstraintRefType}}
end

function test_quadexpr_constraints(ModelType, ::Any)
    model = ModelType()
    @variable(model, x)
    @variable(model, y)

    cref = @constraint(model, x^2 + x <= 1)
    c = JuMP.constraint_object(cref)
    @test JuMP.isequal_canonical(c.func, x^2 + x)
    @test c.set == MOI.LessThan(1.0)

    cref = @constraint(model, y * x - 1.0 == 0.0)
    c = JuMP.constraint_object(cref)
    @test JuMP.isequal_canonical(c.func, x * y)
    @test c.set == MOI.EqualTo(1.0)

    cref = @constraint(
        model,
        [2x - 4x * y + 3x^2 - 1, -3y + 2x * y - 2x^2 + 1] in SecondOrderCone()
    )
    c = JuMP.constraint_object(cref)
    @test JuMP.isequal_canonical(c.func[1], -1 + 3x^2 - 4x * y + 2x)
    @test JuMP.isequal_canonical(c.func[2], 1 - 2x^2 + 2x * y - 3y)
    @test c.set == MOI.SecondOrderCone(2)
end

function test_syntax_error_constraint(ModelType, ::Any)
    model = ModelType()
    @variable(model, x[1:2])
    err = ErrorException(
        "In `@constraint(model, [3, x] in SecondOrderCone())`: unable to " *
        "add the constraint because we don't recognize $([3, x]) as a " *
        "valid JuMP function.",
    )
    @test_throws_strip err @constraint(model, [3, x] in SecondOrderCone())
end

function test_indicator_constraint(ModelType, ::Any)
    model = ModelType()
    @variable(model, a, Bin)
    @variable(model, b, Bin)
    @variable(model, x)
    @variable(model, y)
    for cref in [
        @constraint(model, a => {x + 2y <= 1}),
        @constraint(model, a ⇒ {x + 2y ≤ 1})
    ]
        c = JuMP.constraint_object(cref)
        @test c.func == [a, x + 2y]
        @test c.set == MOI.IndicatorSet{MOI.ACTIVATE_ON_ONE}(MOI.LessThan(1.0))
    end
    for cref in [
        @constraint(model, !b => {2x + y <= 1})
        @constraint(model, ¬b ⇒ {2x + y ≤ 1})
        # This returns a vector of constraints that is concatenated.
        @constraint(model, ![b, b] .=> {[2x + y, 2x + y] .≤ 1})
    ]
        c = JuMP.constraint_object(cref)
        @test c.func == [b, 2x + y]
        @test c.set == MOI.IndicatorSet{MOI.ACTIVATE_ON_ZERO}(MOI.LessThan(1.0))
    end
    err =
        ErrorException("In `@constraint(model, !(a, b) => {x <= 1})`: Invalid binary variable expression `!(a, b)` for indicator constraint.")
    @test_macro_throws err @constraint(model, !(a, b) => {x <= 1})
    err =
        ErrorException("In `@constraint(model, a => x)`: Invalid right-hand side `x` of indicator constraint. Expected constraint surrounded by `{` and `}`.")
    @test_macro_throws err @constraint(model, a => x)
    err =
        ErrorException("In `@constraint(model, a => x <= 1)`: Invalid right-hand side `x <= 1` of indicator constraint. Expected constraint surrounded by `{` and `}`.")
    @test_macro_throws err @constraint(model, a => x <= 1)
    err =
        ErrorException("In `@constraint(model, a => {x <= 1, x >= 0})`: Invalid right-hand side `{x <= 1, x >= 0}` of indicator constraint. Expected constraint surrounded by `{` and `}`.")
    @test_macro_throws err @constraint(model, a => {x <= 1, x >= 0})
    err =
        ErrorException("In `@constraint(model, [a, b] .=> {x + y <= 1})`: Inconsistent use of `.` in symbols to indicate vectorization.")
    @test_macro_throws err @constraint(model, [a, b] .=> {x + y <= 1})
end

function test_SDP_constraint(ModelType, VariableRefType)
    m = ModelType()
    @variable(m, x)
    @variable(m, y)
    @variable(m, z)
    @variable(m, w)

    cref = @constraint(m, [x y; z w] in PSDCone())
    c = JuMP.constraint_object(cref)
    @test c.func == [x, z, y, w]
    @test c.set == MOI.PositiveSemidefiniteConeSquare(2)
    @test c.shape isa JuMP.SquareMatrixShape

    @constraint(m, sym_ref, Symmetric([x 1; 1 -y] - [1 x; x -2]) in PSDCone())
    _test_constraint_name_util(
        sym_ref,
        "sym_ref",
        Vector{AffExpr},
        MOI.PositiveSemidefiniteConeTriangle,
    )
    c = JuMP.constraint_object(sym_ref)
    @test JuMP.isequal_canonical(c.func[1], x - 1)
    @test JuMP.isequal_canonical(c.func[2], 1 - x)
    @test JuMP.isequal_canonical(c.func[3], 2 - y)
    @test c.set == MOI.PositiveSemidefiniteConeTriangle(2)
    @test c.shape isa JuMP.SymmetricMatrixShape

    @SDconstraint(m, cref, [x 1; 1 -y] ⪰ [1 x; x -2])
    _test_constraint_name_util(
        cref,
        "cref",
        Vector{AffExpr},
        MOI.PositiveSemidefiniteConeSquare,
    )
    c = JuMP.constraint_object(cref)
    @test JuMP.isequal_canonical(c.func[1], x - 1)
    @test JuMP.isequal_canonical(c.func[2], 1 - x)
    @test JuMP.isequal_canonical(c.func[3], 1 - x)
    @test JuMP.isequal_canonical(c.func[4], 2 - y)
    @test c.set == MOI.PositiveSemidefiniteConeSquare(2)
    @test c.shape isa JuMP.SquareMatrixShape

    @SDconstraint(m, iref[i = 1:2], 0 ⪯ [x+i x+y; x+y -y])
    for i in 1:2
        _test_constraint_name_util(
            iref[i],
            "iref[$i]",
            Vector{AffExpr},
            MOI.PositiveSemidefiniteConeSquare,
        )
        c = JuMP.constraint_object(iref[i])
        @test JuMP.isequal_canonical(c.func[1], x + i)
        @test JuMP.isequal_canonical(c.func[2], x + y)
        @test JuMP.isequal_canonical(c.func[3], x + y)
        @test JuMP.isequal_canonical(c.func[4], -y)
        @test c.set == MOI.PositiveSemidefiniteConeSquare(2)
        @test c.shape isa JuMP.SquareMatrixShape
    end

    @SDconstraint(m, con_d, 0 ⪯ Diagonal([x, y]))
    c = JuMP.constraint_object(con_d)
    @test c.func isa Vector{GenericAffExpr{Float64,VariableRefType}}
    @test JuMP.isequal_canonical(c.func[1], 1x)
    @test iszero(c.func[2])
    @test iszero(c.func[3])
    @test JuMP.isequal_canonical(c.func[4], 1y)
    @test c.set == MOI.PositiveSemidefiniteConeSquare(2)

    @SDconstraint(m, con_d_sym, 0 ⪯ Symmetric(Diagonal([x, y])))
    c = JuMP.constraint_object(con_d_sym)
    @test c.func isa Vector{GenericAffExpr{Float64,VariableRefType}}
    @test JuMP.isequal_canonical(c.func[1], 1x)
    @test iszero(c.func[2])
    @test JuMP.isequal_canonical(c.func[3], 1y)
    @test c.set == MOI.PositiveSemidefiniteConeTriangle(2)

    @SDconstraint(m, con_td, Tridiagonal([z], [x, y], [w]) ⪰ 0)
    c = JuMP.constraint_object(con_td)
    @test c.func isa Vector{GenericAffExpr{Float64,VariableRefType}}
    @test JuMP.isequal_canonical(c.func[1], 1x)
    @test JuMP.isequal_canonical(c.func[2], 1z)
    @test JuMP.isequal_canonical(c.func[3], 1w)
    @test JuMP.isequal_canonical(c.func[4], 1y)
    @test c.set == MOI.PositiveSemidefiniteConeSquare(2)

    @SDconstraint(m, con_td_sym, Symmetric(Tridiagonal([z], [x, y], [w])) ⪰ 0)
    c = JuMP.constraint_object(con_td_sym)
    @test c.func isa Vector{GenericAffExpr{Float64,VariableRefType}}
    @test JuMP.isequal_canonical(c.func[1], 1x)
    @test JuMP.isequal_canonical(c.func[2], 1w)
    @test JuMP.isequal_canonical(c.func[3], 1y)
    @test c.set == MOI.PositiveSemidefiniteConeTriangle(2)

    @SDconstraint(m, con_ut, UpperTriangular([x y; z w]) ⪰ 0)
    c = JuMP.constraint_object(con_ut)
    @test c.func isa Vector{GenericAffExpr{Float64,VariableRefType}}
    @test JuMP.isequal_canonical(c.func[1], 1x)
    @test iszero(c.func[2])
    @test JuMP.isequal_canonical(c.func[3], 1y)
    @test JuMP.isequal_canonical(c.func[4], 1w)
    @test c.set == MOI.PositiveSemidefiniteConeSquare(2)

    @SDconstraint(m, con_lt, 0 ⪯ LowerTriangular([x y; z w]))
    c = JuMP.constraint_object(con_lt)
    @test c.func isa Vector{GenericAffExpr{Float64,VariableRefType}}
    @test JuMP.isequal_canonical(c.func[1], 1x)
    @test JuMP.isequal_canonical(c.func[2], 1z)
    @test iszero(c.func[3])
    @test JuMP.isequal_canonical(c.func[4], 1w)
    @test c.set == MOI.PositiveSemidefiniteConeSquare(2)

    # Julia changed how it reports keyword arguments between 1.3 and 1.4!
    err = if VERSION < v"1.4"
        ErrorException("function build_constraint does not accept keyword arguments")
    else
        MethodError
    end
    @test_throws(
        err,
        @SDconstraint(m, [x 1; 1 -y] ⪰ [1 x; x -2], unknown_kw = 1)
    )
    # Invalid sense == in SDP constraint
    @test_macro_throws ErrorException @SDconstraint(
        m,
        [x 1; 1 -y] == [1 x; x -2]
    )
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
    JuMP.set_name(con, "kon")
    @test JuMP.constraint_by_name(model, "con") isa Nothing
    _test_constraint_name_util(
        con,
        "kon",
        GenericQuadExpr{Float64,VariableRefType},
        MOI.EqualTo{Float64},
    )
    y = @constraint(model, kon, [x^2, x] in SecondOrderCone())
    err(name) = ErrorException("Multiple constraints have the name $name.")
    @test_throws err("kon") JuMP.constraint_by_name(model, "kon")
    JuMP.set_name(kon, "con")
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
    JuMP.set_name(con, "con")
    @test_throws err("con") JuMP.constraint_by_name(model, "con")
    @test JuMP.constraint_by_name(model, "kon") isa Nothing
    JuMP.set_name(kon, "kon")
    _test_constraint_name_util(
        con,
        "con",
        GenericQuadExpr{Float64,VariableRefType},
        MOI.EqualTo{Float64},
    )
    return _test_constraint_name_util(
        kon,
        "kon",
        Vector{GenericQuadExpr{Float64,VariableRefType}},
        MOI.SecondOrderCone,
    )
end

function test_PSD_constraint_errors(ModelType, ::Any)
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
end

function test_matrix_constraint_errors(ModelType, VariableRefType)
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
end

function test_nonsensical_SDP_constraint(ModelType, ::Any)
    m = ModelType()
    @test_throws_strip(
        ErrorException("In `@variable(m, unequal[1:5, 1:6], PSD)`: Symmetric variables must be square. Got size (5, 6)."),
        @variable(m, unequal[1:5, 1:6], PSD)
    )
    # Some of these errors happen at compile time, so we can't use @test_throws
    @test_throws MethodError @variable(m, notone[1:5, 2:6], PSD)
    @test_throws MethodError @variable(m, oneD[1:5], PSD)
    @test_throws MethodError @variable(m, threeD[1:5, 1:5, 1:5], PSD)
    @test_throws MethodError @variable(m, psd[2] <= rand(2, 2), PSD)

    @test_throws_strip(
        ErrorException("In `@variable(m, -(ones(3, 4)) <= foo[1:4, 1:4] <= ones(4, 4), PSD)`: Non-symmetric bounds, integrality or starting values for symmetric variable."),
        @variable(m, -ones(3, 4) <= foo[1:4, 1:4] <= ones(4, 4), PSD)
    )
    @test_throws_strip(
        ErrorException("In `@variable(m, -(ones(3, 4)) <= foo[1:4, 1:4] <= ones(4, 4), Symmetric)`: Non-symmetric bounds, integrality or starting values for symmetric variable."),
        @variable(m, -ones(3, 4) <= foo[1:4, 1:4] <= ones(4, 4), Symmetric)
    )
    @test_throws_strip(
        ErrorException("In `@variable(m, -(ones(4, 4)) <= foo[1:4, 1:4] <= ones(4, 5), Symmetric)`: Non-symmetric bounds, integrality or starting values for symmetric variable."),
        @variable(m, -ones(4, 4) <= foo[1:4, 1:4] <= ones(4, 5), Symmetric)
    )
    @test_throws_strip(
        ErrorException("In `@variable(m, -(rand(5, 5)) <= nonsymmetric[1:5, 1:5] <= rand(5, 5), Symmetric)`: Non-symmetric bounds, integrality or starting values for symmetric variable."),
        @variable(
            m,
            -rand(5, 5) <= nonsymmetric[1:5, 1:5] <= rand(5, 5),
            Symmetric
        )
    )
end

function test_sum_constraint(ModelType, ::Any)
    model = ModelType()
    @variable(model, x[1:3, 1:3])
    @variable(model, y)
    C = [1 2 3; 4 5 6; 7 8 9]

    @test_expression sum(C[i, j] * x[i, j] for i in 1:2, j in 2:3)
    @test_expression sum(C[i, j] * x[i, j] for i = 1:3, j in 1:3 if i != j) - y
    @test JuMP.isequal_canonical(
        @expression(model, sum(C[i, j] * x[i, j] for i in 1:3, j in 1:i)),
        sum(C[i, j] * x[i, j] for i in 1:3 for j in 1:i),
    )
    @test_expression sum(C[i, j] * x[i, j] for i in 1:3 for j in 1:i)
    @test_expression sum(C[i, j] * x[i, j] for i = 1:3 if true for j in 1:i)
    @test_expression sum(
        C[i, j] * x[i, j] for i = 1:3 if true for j = 1:i if true
    )
    @test_expression sum(0 * x[i, 1] for i in 1:3)
    @test_expression sum(0 * x[i, 1] + y for i in 1:3)
    @test_expression sum(0 * x[i, 1] + y for i in 1:3 for j in 1:3)
end

function test_Model_all_constraints_scalar(::Any, ::Any)
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
    err =
        ErrorException("`GenericAffExpr` is not a concrete type. Did you miss a type parameter?")
    @test_throws err try
        num_constraints(model, JuMP.GenericAffExpr, MOI.ZeroOne)
    catch e
        error(replace(e.msg, "JuMP." => ""))
    end
    @test_throws err try
        all_constraints(model, JuMP.GenericAffExpr, MOI.ZeroOne)
    catch e
        error(replace(e.msg, "JuMP." => ""))
    end
end

function test_Model_all_constraints_vector(::Any, ::Any)
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
        "`$(GenericAffExpr{Float64})` is not a concrete type. Did you miss a " *
        "type parameter?",
    )
    @test_throws err try
        num_constraints(
            model,
            Vector{GenericAffExpr{Float64}},
            MOI.PositiveSemidefiniteConeTriangle,
        )
    catch e
        error(replace(e.msg, "JuMP." => ""))
    end
    @test_throws err try
        all_constraints(
            model,
            Vector{GenericAffExpr{Float64}},
            MOI.SecondOrderCone,
        )
    catch e
        error(replace(e.msg, "JuMP." => ""))
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
end

function test_Model_list_of_constraint_types(::Any, ::Any)
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
end

function test_Model_dual_start(::Any, ::Any)
    model = Model()
    @variable(model, x)
    con = @constraint(model, 2x <= 1)
    @test dual_start_value(con) === nothing
    set_dual_start_value(con, 2)
    @test dual_start_value(con) == 2.0
    set_dual_start_value(con, nothing)
    @test dual_start_value(con) === nothing
end

function test_Model_dual_start_vector(::Any, ::Any)
    model = Model()
    @variable(model, x)
    con_vec = @constraint(model, [x, x] in SecondOrderCone())
    @test dual_start_value(con_vec) === nothing
    set_dual_start_value(con_vec, [1.0, 3.0])
    @test dual_start_value(con_vec) == [1.0, 3.0]
    set_dual_start_value(con_vec, nothing)
    @test dual_start_value(con_vec) === nothing
end

function test_Model_change_coefficient(::Any, ::Any)
    model = JuMP.Model()
    x = @variable(model)
    con_ref = @constraint(model, 2 * x == -1)
    @test JuMP.normalized_coefficient(con_ref, x) == 2.0
    JuMP.set_normalized_coefficient(con_ref, x, 1.0)
    @test JuMP.normalized_coefficient(con_ref, x) == 1.0
    JuMP.set_normalized_coefficient(con_ref, x, 3)  # Check type promotion.
    @test JuMP.normalized_coefficient(con_ref, x) == 3.0
    quad_con = @constraint(model, x^2 == 0)
    @test JuMP.normalized_coefficient(quad_con, x) == 0.0
    JuMP.set_normalized_coefficient(quad_con, x, 2)
    @test JuMP.normalized_coefficient(quad_con, x) == 2.0
    @test JuMP.isequal_canonical(
        JuMP.constraint_object(quad_con).func,
        x^2 + 2x,
    )
end

function test_Model_change_rhs(::Any, ::Any)
    model = JuMP.Model()
    x = @variable(model)
    con_ref = @constraint(model, 2 * x <= 1)
    @test JuMP.normalized_rhs(con_ref) == 1.0
    JuMP.set_normalized_rhs(con_ref, 2.0)
    @test JuMP.normalized_rhs(con_ref) == 2.0
    con_ref = @constraint(model, 2 * x - 1 == 1)
    @test JuMP.normalized_rhs(con_ref) == 2.0
    JuMP.set_normalized_rhs(con_ref, 3)
    @test JuMP.normalized_rhs(con_ref) == 3.0
    con_ref = @constraint(model, 0 <= 2 * x <= 1)
    @test_throws MethodError JuMP.set_normalized_rhs(con_ref, 3)
end

function test_Model_add_to_function_constant_scalar(::Any, ::Any)
    model = JuMP.Model()
    x = @variable(model)
    con_ref = @constraint(model, 2 <= 2 * x <= 3)
    con = constraint_object(con_ref)
    @test JuMP.isequal_canonical(JuMP.jump_function(con), 2x)
    @test JuMP.moi_set(con) == MOI.Interval(2.0, 3.0)
    JuMP.add_to_function_constant(con_ref, 1.0)
    con = constraint_object(con_ref)
    @test JuMP.isequal_canonical(JuMP.jump_function(con), 2x)
    @test JuMP.moi_set(con) == MOI.Interval(1.0, 2.0)
end

function test_Model_add_to_function_constant_vector(::Any, ::Any)
    model = JuMP.Model()
    x = @variable(model)
    con_ref = @constraint(model, [x + 1, x - 1] in MOI.Nonnegatives(2))
    con = constraint_object(con_ref)
    @test JuMP.isequal_canonical(JuMP.jump_function(con), [x + 1, x - 1])
    @test JuMP.moi_set(con) == MOI.Nonnegatives(2)
    JuMP.add_to_function_constant(con_ref, [2, 3])
    con = constraint_object(con_ref)
    @test JuMP.isequal_canonical(JuMP.jump_function(con), [x + 3, x + 2])
    @test JuMP.moi_set(con) == MOI.Nonnegatives(2)
end

function test_Model_value_constraint_var(ModelType, ::Any)
    model = ModelType()
    @variable(model, x[1:2])
    @constraint(model, c1, x[1] + x[2] <= 3.0)
    @constraint(model, c2, x[1]^2 + x[2]^2 <= 3.0)
    @constraint(model, c3, [1.0, x[1], x[2]] in SecondOrderCone())

    vals = Dict(x[1] => 1.0, x[2] => 2.0)
    f = vidx -> vals[vidx]

    @test value(c1, f) === 3.0 # Affine expression
    @test value(c2, f) === 5.0 # Quadratic expression
    @test value(c3, f) == [1.0, 1.0, 2.0] # Vector expression
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
    JuMP.optimize!(model)
    mock_optimizer = JuMP.backend(model).optimizer.model
    MOI.set(mock_optimizer, MOI.TerminationStatus(), MOI.OPTIMAL)
    MOI.set(mock_optimizer, MOI.DualStatus(), MOI.FEASIBLE_POINT)
    JuMP.optimize!(model)
    for constraint_name in keys(constraint_dual)
        ci = MOI.get(backend(model), MOI.ConstraintIndex, constraint_name)
        constraint_ref = JuMP.ConstraintRef(model, ci, JuMP.ScalarShape())
        MOI.set(
            mock_optimizer,
            MOI.ConstraintDual(),
            JuMP.optimizer_index(constraint_ref),
            constraint_dual[constraint_name],
        )
        @test dual(constraint_ref) == constraint_dual[constraint_name]
        @test shadow_price(constraint_ref) == constraint_shadow[constraint_name]
    end
end

function test_Model_shadow_price(::Any, ::Any)
    _test_shadow_price_util(
        """
        variables: x, y
        minobjective: -1.0*x
        xub: x <= 2.0
        ylb: y >= 0.0
        c: x + y <= 1.0
        """,
        Dict("xub" => 0.0, "ylb" => 1.0, "c" => -1.0),
        Dict("xub" => 0.0, "ylb" => -1.0, "c" => -1.0),
    )

    _test_shadow_price_util(
        """
        variables: x, y
        maxobjective: 1.0*x
        xub: x <= 2.0
        ylb: y >= 0.0
        c: x + y <= 1.0
        """,
        Dict("xub" => 0.0, "ylb" => 1.0, "c" => -1.0),
        Dict("xub" => 0.0, "ylb" => 1.0, "c" => 1.0),
    )

    _test_shadow_price_util(
        """
        variables: x, y
        maxobjective: 1.0*x
        xub: x <= 2.0
        ylb: y >= 0.0
        """,
        Dict("xub" => -1.0, "ylb" => 0.0),
        Dict("xub" => 1.0, "ylb" => 0.0),
    )

    _test_shadow_price_util(
        """
        variables: x
        maxobjective: 1.0*x
        xeq: x == 2.0
        """,
        Dict("xeq" => -1.0),
        Dict("xeq" => 1.0),
    )

    return _test_shadow_price_util(
        """
        variables: x
        minobjective: 1.0*x
        xeq: x == 2.0
        """,
        Dict("xeq" => 1.0),
        Dict("xeq" => -1.0),
    )
end

function test_abstractarray_vector_constraint(ModelType, ::Any)
    model = ModelType()
    @variable(model, x[1:2, 1:2])
    c = @constraint(model, view(x, 1:4) in SOS1())
    obj = constraint_object(c)
    @test obj.func == x[1:4]
    @test obj.set == MOI.SOS1([1.0, 2.0, 3.0, 4.0])
end

function runtests()
    for name in names(@__MODULE__; all = true)
        if !startswith("$(name)", "test_")
            continue
        end
        f = getfield(@__MODULE__, name)
        @testset "$(name)" begin
            f(Model, VariableRef)
        end
        if !startswith("$(name)", "test_Model_")
            @testset "$(name)-JuMPExtension" begin
                f(JuMPExtension.MyModel, JuMPExtension.MyVariableRef)
            end
        end
    end
end

end

TestConstraint.runtests()

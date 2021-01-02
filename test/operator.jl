module TestOperators

using JuMP
using LinearAlgebra
using SparseArrays
using Test

const MA = JuMP._MA

include(joinpath(@__DIR__, "utilities.jl"))
include(joinpath(@__DIR__, "JuMPExtension.jl"))

# For "DimensionMismatch when performing vector-matrix multiplication with custom types #988"
import Base: +, *
struct MyType{T}
    a::T
end
struct MySumType{T}
    a::T
end
Base.copy(t::MyType) = t
Base.zero(::Type{MyType{T}}) where {T} = MyType(zero(T))
Base.one(::Type{MyType{T}}) where {T} = MyType(one(T))
Base.zero(::Type{MySumType{T}}) where {T} = MySumType(zero(T))
Base.zero(::MySumType{T}) where {T} = MySumType(zero(T))
Base.transpose(t::MyType) = MyType(t.a)
Base.transpose(t::MySumType) = MySumType(t.a)
LinearAlgebra.adjoint(t::Union{MyType,MySumType}) = t
function +(
    t1::MyT,
    t2::MyS,
) where {MyT<:Union{MyType,MySumType},MyS<:Union{MyType,MySumType}}
    return MySumType(t1.a + t2.a)
end
*(t1::MyType{S}, t2::T) where {S,T} = MyType(t1.a * t2)
*(t1::S, t2::MyType{T}) where {S,T} = MyType(t1 * t2.a)
*(t1::MyType{S}, t2::MyType{T}) where {S,T} = MyType(t1.a * t2.a)

function JuMP.isequal_canonical(t::MySumType, s::MySumType)
    return JuMP.isequal_canonical(t.a, s.a)
end
function JuMP.isequal_canonical(
    x::AbstractArray{<:MySumType},
    y::AbstractArray{<:MySumType},
)
    return size(x) == size(y) && all(JuMP.isequal_canonical.(x, y))
end

function test_promotion(ModelType, VariableRefType)
    m = ModelType()
    I = Int
    V = VariableRefType
    A = JuMP.GenericAffExpr{Float64,VariableRefType}
    Q = JuMP.GenericQuadExpr{Float64,VariableRefType}
    @test promote_type(V, I) == A
    @test promote_type(I, V) == A
    @test promote_type(A, I) == A
    @test promote_type(I, A) == A
    @test promote_type(A, V) == A
    @test promote_type(V, A) == A
    @test promote_type(Q, I) == Q
    @test promote_type(I, Q) == Q
    @test promote_type(Q, A) == Q
    @test promote_type(A, Q) == Q
    @test promote_type(Q, V) == Q
    @test promote_type(V, Q) == Q
end

function test_broadcast_division_error(ModelType, ::Any)
    model = ModelType()
    @variable(model, x[1:2, 1:2])
    A = [1 2; 3 4]
    B = sparse(A)
    y = SparseMatrixCSC(2, 2, copy(B.colptr), copy(B.rowval), vec(x))
    @test_throws ErrorException A ./ x
    @test_throws ErrorException B ./ x
    @test_throws ErrorException A ./ y
    @test_throws ErrorException B ./ y

    # TODO: Refactor to avoid calling the internal JuMP function
    # `_densify_with_jump_eltype`.
    #z = JuMP._densify_with_jump_eltype((2 .* y) ./ 3)
    #@test JuMP.isequal_canonical((2 .* x) ./ 3, z)
    #z = JuMP._densify_with_jump_eltype(2 * (y ./ 3))
    #@test JuMP.isequal_canonical(2 .* (x ./ 3), z)
    #z = JuMP._densify_with_jump_eltype((x[1,1],) .* B)
    #@test JuMP.isequal_canonical((x[1,1],) .* A, z)
end

function test_vectorized_comparisons(ModelType, ::Any)
    m = ModelType()
    @variable(m, x[1:3])
    A = [
        1 2 3
        0 4 5
        6 0 7
    ]
    B = sparse(A)
    # force vector output
    cref1 = @constraint(m, reshape(x, (1, 3)) * A * x .>= 1)
    c1 = JuMP.constraint_object.(cref1)
    f1 = map(c -> c.func, c1)
    @test JuMP.isequal_canonical(
        f1,
        [
            x[1] * x[1] +
            2x[1] * x[2] +
            4x[2] * x[2] +
            9x[1] * x[3] +
            5x[2] * x[3] +
            7x[3] * x[3],
        ],
    )
    @test all(c -> c.set.lower == 1, c1)

    cref2 = @constraint(m, x' * A * x >= 1)
    c2 = JuMP.constraint_object.(cref2)
    @test JuMP.isequal_canonical(f1[1], c2.func)

    mat = [
        3x[1] + 12x[3] + 4x[2]
        2x[1] + 12x[2] + 10x[3]
        15x[1] + 5x[2] + 21x[3]
    ]

    cref3 = @constraint(m, (x'A)' + 2A * x .<= 1)
    c3 = JuMP.constraint_object.(cref3)
    f3 = map(c -> c.func, c3)
    @test JuMP.isequal_canonical(f3, mat)
    @test all(c -> c.set.upper == 1, c3)
    @test JuMP.isequal_canonical((x'A)' + 2A * x, (x'A)' + 2B * x)
    @test JuMP.isequal_canonical((x'A)' + 2A * x, (x'B)' + 2A * x)
    @test JuMP.isequal_canonical((x'A)' + 2A * x, (x'B)' + 2B * x)
    @test JuMP.isequal_canonical((x'A)' + 2A * x, MA.@rewrite((x'A)' + 2A * x))
    @test JuMP.isequal_canonical((x'A)' + 2A * x, MA.@rewrite((x'B)' + 2A * x))
    @test JuMP.isequal_canonical((x'A)' + 2A * x, MA.@rewrite((x'A)' + 2B * x))
    @test JuMP.isequal_canonical((x'A)' + 2A * x, MA.@rewrite((x'B)' + 2B * x))

    cref4 = @constraint(m, -1 .<= (x'A)' + 2A * x .<= 1)
    c4 = JuMP.constraint_object.(cref4)
    f4 = map(c -> c.func, c4)
    @test JuMP.isequal_canonical(f4, mat)
    @test all(c -> c.set.lower == -1, c4)
    @test all(c -> c.set.upper == 1, c4)

    cref5 = @constraint(m, -[1:3;] .<= (x'A)' + 2A * x .<= 1)
    c5 = JuMP.constraint_object.(cref5)
    f5 = map(c -> c.func, c5)
    @test JuMP.isequal_canonical(f5, mat)
    @test map(c -> c.set.lower, c5) == -[1:3;]
    @test all(c -> c.set.upper == 1, c4)

    cref6 = @constraint(m, -[1:3;] .<= (x'A)' + 2A * x .<= [3:-1:1;])
    c6 = JuMP.constraint_object.(cref6)
    f6 = map(c -> c.func, c6)
    @test JuMP.isequal_canonical(f6, mat)
    @test map(c -> c.set.lower, c6) == -[1:3;]
    @test map(c -> c.set.upper, c6) == [3:-1:1;]

    cref7 = @constraint(m, -[1:3;] .<= (x'A)' + 2A * x .<= 3)
    c7 = JuMP.constraint_object.(cref7)
    f7 = map(c -> c.func, c7)
    @test JuMP.isequal_canonical(f7, mat)
    @test map(c -> c.set.lower, c7) == -[1:3;]
    @test all(c -> c.set.upper == 3, c7)
end

function test_custom_dimension_mismatch(ModelType, VariableRefType)
    ElemT = MySumType{JuMP.GenericAffExpr{Float64,VariableRefType}}
    model = ModelType()
    @variable(model, Q[1:3, 1:3], PSD)
    x = [MyType(1), MyType(2), MyType(3)]
    y = Q * x
    z = x' * Q
    @test y isa Vector{ElemT}
    @test size(y) == (3,)
    @test z isa Adjoint{ElemT,<:Any}
    @test size(z) == (1, 3)
    for i in 1:3
        # Q is symmetric
        a = zero(JuMP.GenericAffExpr{Float64,VariableRefType})
        a += Q[1, i]
        a += 2Q[2, i]
        a += 3Q[3, i]
        # Q[1,i] + 2Q[2,i] + 3Q[3,i] is rearranged as 2 Q[2,3] + Q[1,3] + 3 Q[3,3]
        @test z[i].a == y[i].a == a
    end
end

function test_matrix_multiplication(ModelType, VariableRefType)
    ElemT = MySumType{JuMP.GenericAffExpr{Float64,VariableRefType}}
    model = ModelType()
    @variable(model, Q[1:3, 1:3], PSD)
    X = MyType.((1:3)' .+ (1:3))
    @test_expression Q * X
    Y = Q * X
    @test_expression X' * Q
    Z = X' * Q
    @test Y isa Matrix{ElemT}
    @test size(Y) == (3, 3)
    @test Z isa Matrix{ElemT}
    @test size(Z) == (3, 3)
    @test JuMP.isequal_canonical(Z', Y)
    @test_expression Q * X'
    Y = Q * X'
    @test_expression X * Q
    Z = X * Q
    @test Y isa Matrix{ElemT}
    @test size(Y) == (3, 3)
    @test Z isa Matrix{ElemT}
    @test size(Z) == (3, 3)
    @test JuMP.isequal_canonical(Z', Y)
end

function test_operator_warn(ModelType, ::Any)
    model = ModelType()
    @variable model x[1:51]
    # JuMPExtension does not have the `operator_counter` field
    if ModelType <: Model
        @test model.operator_counter == 0
    end
    # Triggers the increment of operator_counter since sum(x) has more than 50 terms
    @test_expression(sum(x) + 2x[1])
    if ModelType <: Model
        # The following check verifies that this test covers the code incrementing `operator_counter`
        @test model.operator_counter == 1
    end
end

function test_uniform_scaling(ModelType, ::Any)
    model = ModelType()
    @variable(model, x)
    @test_expression_with_string x + 2I "x + 2"
    @test_expression_with_string (x + 1) + I "x + 2"
    @test_expression_with_string x - 2I "x - 2"
    @test_expression_with_string (x - 1) - I "x - 2"
    @test_expression_with_string 2I + x "x + 2"
    @test_expression_with_string I + (x + 1) "x + 2"
    @test_expression_with_string 2I - x "-x + 2"
    @test_expression_with_string I - (x - 1) "-x + 2"
    @test_expression_with_string I * x "x"
    @test_expression_with_string I * (x + 1) "x + 1"
    @test_expression_with_string (x + 1) * I "x + 1"
end

function test_basic_operators(ModelType, ::Any)
    model = ModelType()
    @variable(model, w)
    @variable(model, x)
    @variable(model, y)
    @variable(model, z)

    aff = @inferred 7.1 * x + 2.5
    @test_expression_with_string 7.1 * x + 2.5 "7.1 x + 2.5"
    aff2 = @inferred 1.2 * y + 1.2
    @test_expression_with_string 1.2 * y + 1.2 "1.2 y + 1.2"
    q = @inferred 2.5 * y * z + aff
    @test_expression_with_string 2.5 * y * z + aff "2.5 y*z + 7.1 x + 2.5"
    q2 = @inferred 8 * x * z + aff2
    @test_expression_with_string 8 * x * z + aff2 "8 x*z + 1.2 y + 1.2"
    @test_expression_with_string 2 * x * x + 1 * y * y + z + 3 "2 x² + y² + z + 3"
end

function test_basic_operators_number(ModelType, ::Any)
    model = ModelType()
    @variable(model, w)
    @variable(model, x)
    @variable(model, y)
    @variable(model, z)

    aff = @inferred 7.1 * x + 2.5
    q = @inferred 2.5 * y * z + aff

    # 1-1 Number--Number
    # ???
    # 1-2 Number--Variable
    @test_expression_with_string 4.13 + w "w + 4.13"
    @test_expression_with_string 3.16 - w "-w + 3.16"
    @test_expression_with_string 5.23 * w "5.23 w"
    @test_throws ErrorException 2.94 / w
    # 1-3 Number--AffExpr
    @test_expression_with_string 1.5 + aff "7.1 x + 4"
    @test_expression_with_string 1.5 - aff "-7.1 x - 1"
    @test_expression_with_string 2 * aff "14.2 x + 5"
    @test_throws ErrorException 2 / aff
    # 1-4 Number--QuadExpr
    @test_expression_with_string 1.5 + q "2.5 y*z + 7.1 x + 4"
    @test_expression_with_string 1.5 - q "-2.5 y*z - 7.1 x - 1"
    @test_expression_with_string 2 * q "5 y*z + 14.2 x + 5"
    @test_throws ErrorException 2 / q
end

function test_basic_operators_variable(ModelType, ::Any)
    model = ModelType()
    @variable(model, w)
    @variable(model, x)
    @variable(model, y)
    @variable(model, z)

    aff = @inferred 7.1 * x + 2.5
    q = @inferred 2.5 * y * z + aff

    # 2-0 Variable unary
    @test (+x) === x
    @test_expression_with_string -x "-x"
    # 2-1 Variable--Number
    @test_expression_with_string w + 4.13 "w + 4.13"
    @test_expression_with_string w - 4.13 "w - 4.13"
    @test_expression_with_string w * 4.13 "4.13 w"
    @test_expression_with_string w / 2.00 "0.5 w"
    @test w == w
    @test_expression_with_string x * y - 1 "x*y - 1"
    @test_expression_with_string x^2 "x²"
    @test_expression_with_string x^1 "x"
    @test_expression_with_string x^0 "1"
    @test_throws ErrorException x^3
    @test_throws ErrorException x^1.5
    # 2-2 Variable--Variable
    @test_expression_with_string w + x "w + x"
    @test_expression_with_string w - x "w - x"
    @test_expression_with_string w * x "w*x"
    @test_expression_with_string x - x "0"
    @test_throws ErrorException w / x
    @test_expression_with_string y * z - x "y*z - x"
    # 2-3 Variable--AffExpr
    @test_expression_with_string z + aff "z + 7.1 x + 2.5"
    @test_expression_with_string z - aff "z - 7.1 x - 2.5"
    @test_expression_with_string z * aff "7.1 z*x + 2.5 z"
    @test_throws ErrorException z / aff
    @test_throws MethodError z ≤ aff
    @test_expression_with_string 7.1 * x - aff "0 x - 2.5"
    # 2-4 Variable--QuadExpr
    @test_expression_with_string w + q "2.5 y*z + w + 7.1 x + 2.5"
    @test_expression_with_string w - q "-2.5 y*z + w - 7.1 x - 2.5"
    @test_throws ErrorException w * q
    @test_throws ErrorException w / q
    @test transpose(x) === x
    @test conj(x) === x
end

function test_basic_operators_affexpr(ModelType, ::Any)
    model = ModelType()
    @variable(model, w)
    @variable(model, x)
    @variable(model, y)
    @variable(model, z)

    aff = @inferred 7.1 * x + 2.5
    aff2 = @inferred 1.2 * y + 1.2
    q = @inferred 2.5 * y * z + aff
    # 3-0 AffExpr unary
    @test_expression_with_string +aff "7.1 x + 2.5"
    @test_expression_with_string -aff "-7.1 x - 2.5"
    # 3-1 AffExpr--Number
    @test_expression_with_string aff + 1.5 "7.1 x + 4"
    @test_expression_with_string aff - 1.5 "7.1 x + 1"
    @test_expression_with_string aff * 2 "14.2 x + 5"
    @test_expression_with_string aff / 2 "3.55 x + 1.25"
    @test_throws MethodError aff ≤ 1
    @test aff == aff
    @test_throws MethodError aff ≥ 1
    @test_expression_with_string aff - 1 "7.1 x + 1.5"
    @test_expression_with_string aff^2 "50.41 x² + 35.5 x + 6.25"
    @test_expression_with_string (7.1 * x + 2.5)^2 "50.41 x² + 35.5 x + 6.25"
    @test_expression_with_string aff^1 "7.1 x + 2.5"
    @test_expression_with_string (7.1 * x + 2.5)^1 "7.1 x + 2.5"
    @test_expression_with_string aff^0 "1"
    @test_expression_with_string (7.1 * x + 2.5)^0 "1"
    @test_throws ErrorException aff^3
    @test_throws ErrorException (7.1 * x + 2.5)^3
    @test_throws ErrorException aff^1.5
    @test_throws ErrorException (7.1 * x + 2.5)^1.5
    # 3-2 AffExpr--Variable
    @test_expression_with_string aff + z "7.1 x + z + 2.5"
    @test_expression_with_string aff - z "7.1 x - z + 2.5"
    @test_expression_with_string aff * z "7.1 x*z + 2.5 z"
    @test_throws ErrorException aff / z
    @test_expression_with_string aff - 7.1 * x "0 x + 2.5"
    # 3-3 AffExpr--AffExpr
    @test_expression_with_string aff + aff2 "7.1 x + 1.2 y + 3.7"
    @test_expression_with_string aff - aff2 "7.1 x - 1.2 y + 1.3"
    @test_expression_with_string aff * aff2 "8.52 x*y + 3 y + 8.52 x + 3"
    @test string((x + x) * (x + 3)) == string((x + 3) * (x + x))  # Issue #288
    @test_throws ErrorException aff / aff2
    @test_expression_with_string aff - aff "0 x"
    # 4-4 AffExpr--QuadExpr
    @test_expression_with_string aff2 + q "2.5 y*z + 1.2 y + 7.1 x + 3.7"
    @test_expression_with_string aff2 - q "-2.5 y*z + 1.2 y - 7.1 x - 1.3"
    @test_throws ErrorException aff2 * q
    @test_throws ErrorException aff2 / q
    @test transpose(aff) === aff
    @test conj(aff) === aff
end

function test_basic_operators_quadexpr(ModelType, ::Any)
    model = ModelType()
    @variable(model, w)
    @variable(model, x)
    @variable(model, y)
    @variable(model, z)

    aff = @inferred 7.1 * x + 2.5
    aff2 = @inferred 1.2 * y + 1.2
    q = @inferred 2.5 * y * z + aff
    q2 = @inferred 8 * x * z + aff2
    # 4-0 QuadExpr unary
    @test_expression_with_string +q "2.5 y*z + 7.1 x + 2.5"
    @test_expression_with_string -q "-2.5 y*z - 7.1 x - 2.5"
    # 4-1 QuadExpr--Number
    @test_expression_with_string q + 1.5 "2.5 y*z + 7.1 x + 4"
    @test_expression_with_string q - 1.5 "2.5 y*z + 7.1 x + 1"
    @test_expression_with_string q * 2 "5 y*z + 14.2 x + 5"
    @test_expression_with_string q / 2 "1.25 y*z + 3.55 x + 1.25"
    @test q == q
    @test_expression_with_string aff2 - q "-2.5 y*z + 1.2 y - 7.1 x - 1.3"
    # 4-2 QuadExpr--Variable
    @test_expression_with_string q + w "2.5 y*z + 7.1 x + w + 2.5"
    @test_expression_with_string q - w "2.5 y*z + 7.1 x - w + 2.5"
    @test_throws ErrorException q * w
    @test_throws ErrorException q / w
    # 4-3 QuadExpr--AffExpr
    @test_expression_with_string q + aff2 "2.5 y*z + 7.1 x + 1.2 y + 3.7"
    @test_expression_with_string q - aff2 "2.5 y*z + 7.1 x - 1.2 y + 1.3"
    @test_throws ErrorException q * aff2
    @test_throws ErrorException q / aff2
    # 4-4 QuadExpr--QuadExpr
    @test_expression_with_string q + q2 "2.5 y*z + 8 x*z + 7.1 x + 1.2 y + 3.7"
    @test_expression_with_string q - q2 "2.5 y*z - 8 x*z + 7.1 x - 1.2 y + 1.3"
    @test_throws ErrorException q * q2
    @test_throws ErrorException q / q2
    @test transpose(q) === q
    @test conj(q) === q
end

function test_dot(ModelType, ::Any)
    model = ModelType()
    @variable(model, 0 ≤ x[1:3] ≤ 1)

    @test_expression_with_string dot(x[1], x[1]) "x[1]²"
    @test_expression_with_string dot(2, x[1]) "2 x[1]"
    @test_expression_with_string dot(x[1], 2) "2 x[1]"

    c = vcat(1:3)
    @test_expression_with_string dot(c, x) "x[1] + 2 x[2] + 3 x[3]"
    @test_expression_with_string dot(x, c) "x[1] + 2 x[2] + 3 x[3]"

    A = [1 3; 2 4]
    @variable(model, 1 ≤ y[1:2, 1:2] ≤ 1)
    @test_expression_with_string dot(A, y) "y[1,1] + 2 y[2,1] + 3 y[1,2] + 4 y[2,2]"
    @test_expression_with_string dot(y, A) "y[1,1] + 2 y[2,1] + 3 y[1,2] + 4 y[2,2]"

    B = ones(2, 2, 2)
    @variable(model, 0 ≤ z[1:2, 1:2, 1:2] ≤ 1)
    @test_expression_with_string dot(B, z) "z[1,1,1] + z[2,1,1] + z[1,2,1] + z[2,2,1] + z[1,1,2] + z[2,1,2] + z[1,2,2] + z[2,2,2]"
    @test_expression_with_string dot(z, B) "z[1,1,1] + z[2,1,1] + z[1,2,1] + z[2,2,1] + z[1,1,2] + z[2,1,2] + z[1,2,2] + z[2,2,2]"

    @objective(model, Max, dot(x, ones(3)) - dot(y, ones(2, 2)))
    for i in 1:3
        JuMP.set_start_value(x[i], 1)
    end
    for i in 1:2, j in 1:2
        JuMP.set_start_value(y[i, j], 1)
    end
    for i in 1:2, j in 1:2, k in 1:2
        JuMP.set_start_value(z[i, j, k], 1)
    end
    @test dot(c, JuMP.start_value.(x)) ≈ 6
    @test dot(A, JuMP.start_value.(y)) ≈ 10
    @test dot(B, JuMP.start_value.(z)) ≈ 8

    if ModelType <: Model
        # Only `Model` is guaranteed to have `operator_counter`, so
        # only test for that case.
        @testset "dot doesn't trigger operator_counter" begin
            # Check that dot is not falling back to default, inefficient
            # addition (JuMP PR #943).
            model = ModelType()
            @test model.operator_counter == 0
            @variable(model, x[1:100])
            JuMP.set_start_value.(x, 1:100)
            @expression(model, test_sum, sum(x[i] * i for i in 1:100))
            @expression(model, test_dot1, dot(x, 1:100))
            @expression(model, test_dot2, dot(1:100, x))
            @test model.operator_counter == 0
            test_add = test_dot1 + test_dot2
            @test model.operator_counter == 1  # Check triggerable.
            test_sum_value = JuMP.value(test_sum, JuMP.start_value)
            @test test_sum_value ≈ JuMP.value(test_dot1, JuMP.start_value)
            @test test_sum_value ≈ JuMP.value(test_dot2, JuMP.start_value)
        end
    end
end

function test_higher_level(ModelType, ::Any)
    model = ModelType()
    @variable(model, 0 ≤ matrix[1:3, 1:3] ≤ 1, start = 1)
    @testset "sum(j::DenseAxisArray{Variable})" begin
        @test_expression_with_string sum(matrix) "matrix[1,1] + matrix[2,1] + matrix[3,1] + matrix[1,2] + matrix[2,2] + matrix[3,2] + matrix[1,3] + matrix[2,3] + matrix[3,3]"
    end

    @testset "sum(j::DenseAxisArray{T}) where T<:Real" begin
        @test sum(JuMP.start_value.(matrix)) ≈ 9
    end
    @testset "sum(j::Array{VariableRef})" begin
        @test string(sum(matrix[1:3, 1:3])) == string(sum(matrix))
    end
    @testset "sum(affs::Array{AffExpr})" begin
        @test_expression_with_string sum([
            2 * matrix[i, j] for i in 1:3, j in 1:3
        ]) "2 matrix[1,1] + 2 matrix[2,1] + 2 matrix[3,1] + 2 matrix[1,2] + 2 matrix[2,2] + 2 matrix[3,2] + 2 matrix[1,3] + 2 matrix[2,3] + 2 matrix[3,3]"
    end
    @testset "sum(quads::Array{QuadExpr})" begin
        @test_expression_with_string sum([
            2 * matrix[i, j]^2 for i in 1:3, j in 1:3
        ]) "2 matrix[1,1]² + 2 matrix[2,1]² + 2 matrix[3,1]² + 2 matrix[1,2]² + 2 matrix[2,2]² + 2 matrix[3,2]² + 2 matrix[1,3]² + 2 matrix[2,3]² + 2 matrix[3,3]²"
    end
    S = [1, 3]
    @variable(model, x[S], start = 1)
    @testset "sum(j::JuMPDict{VariableRef})" begin
        @test_expression sum(x)
        @test length(string(sum(x))) == 11 # order depends on hashing
        @test occursin("x[1]", string(sum(x)))
        @test occursin("x[3]", string(sum(x)))
    end
    @testset "sum(j::JuMPDict{T}) where T<:Real" begin
        @test sum(JuMP.start_value.(x)) == 2
    end
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
        @testset "$(name)-JuMPExtension" begin
            f(JuMPExtension.MyModel, JuMPExtension.MyVariableRef)
        end
    end
end

end

TestOperators.runtests()

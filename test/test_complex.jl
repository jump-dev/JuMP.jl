#  Copyright 2017, Iain Dunning, Joey Huchette, Miles Lubin, and contributors
#  This Source Code Form is subject to the terms of the Mozilla Public
#  License, v. 2.0. If a copy of the MPL was not distributed with this
#  file, You can obtain one at https://mozilla.org/MPL/2.0/.
#############################################################################
# JuMP
# An algebraic modeling language for Julia
# See https://github.com/jump-dev/JuMP.jl
#############################################################################

module TestComplexNumberSupport

using JuMP
using Test

import LinearAlgebra
import MutableArithmetics as MA
import SparseArrays

include(joinpath(@__DIR__, "utilities.jl"))

function _test_dot(a, b)
    @test LinearAlgebra.dot(a, b) == conj(a) * b
    @test LinearAlgebra.dot(b, a) == conj(b) * a
end

function test_complex_aff_expr()
    model = Model()
    @variable(model, x)
    y = (1 + 2im) * x + 1
    @test typeof(y) == GenericAffExpr{Complex{Float64},VariableRef}
    @test 1im + x == 1im + 1 * x
    @test 1im + x == 1im + (1 + 0im) * x
    @test 1im - x == -(1 + 0im)x + 1im
    @test 1im - 1 * x == -1 * x + 1im
    _test_dot(4 - 5im, y)
    _test_dot((3 - 7im) * x + 2 - im, y)
    return
end

function test_complex_quad_expr()
    model = Model()
    @variable(model, x)
    y = im * x^2
    @test typeof(y) == GenericQuadExpr{Complex{Float64},VariableRef}
    @test y == x^2 * im
    @test x^2 + y == y + x^2
    @test x^2 + y == MA.@rewrite(x^2 + y)
    @test x^2 + y == MA.@rewrite(x^2 + 1 * y)
    @test x^2 + y == MA.@rewrite(y + x^2)
    @test x^2 + y == MA.@rewrite(y + 1 * x^2)
    @test x^2 + y == MA.@rewrite(x^2 + im * x^2)

    @test x^2 - y == -(y - x^2)
    @test x^2 - y == MA.@rewrite(x^2 - y)
    @test x^2 - y == MA.@rewrite(x^2 + im * y * im)
    @test x^2 - y == MA.@rewrite(-(y - x^2))

    z = x^2
    @test y + z * im == im * z + y
    @test y + z * im == MA.@rewrite(y + z * im)
    @test y + z * im == MA.@rewrite(im * z + y)
    @test y - x^2 == -x^2 + y
    @test y - x^2 == MA.@rewrite(y - x^2)
    @test y - x^2 == MA.@rewrite(-x^2 + y)
    return
end

function test_complex_plus_variable()
    model = Model()
    @variable(model, x)
    y = x + im
    @test typeof(y) == GenericAffExpr{Complex{Float64},VariableRef}
    @test y == im + x
    return
end

function test_complex_minus_variable()
    model = Model()
    @variable(model, x)
    y = im - x
    @test typeof(y) == GenericAffExpr{Complex{Float64},VariableRef}
    @test -y == x - im
    return
end

function test_complex_aff_expr_convert()
    model = Model()
    @variable(model, x)
    y = (1 + 2im) * x + 1
    y_int = convert(GenericAffExpr{Complex{Int},VariableRef}, y)
    @test typeof(y_int) == GenericAffExpr{Complex{Int},VariableRef}
    @test y_int == y
    @test_throws InexactError convert(AffExpr, y)
    return
end

function test_complex_add_aff()
    model = Model()
    @variable(model, x)
    real_aff = 3x - 1
    complex_aff = (1 + 2im) * x + 1
    @test complex_aff == MA.@rewrite((1 + 2im) * x + 1)
    @test complex_aff == MA.@rewrite(1 + (1 + 2im) * x)
    @test complex_aff == MA.@rewrite(1 + (2im) * x + x)
    @test real_aff + complex_aff == complex_aff + real_aff
    @test real_aff - complex_aff == -(complex_aff - real_aff)
    return
end

function test_complex_vector_constraint()
    model = Model()
    @variable(model, x)
    con_ref = @constraint(model, [(1 + 2im) * x + 1] in MOI.Zeros(1))
    @test list_of_constraint_types(model) ==
          [(Vector{GenericAffExpr{Complex{Float64},VariableRef}}, MOI.Zeros)]
    con_obj = constraint_object(con_ref)
    @test jump_function(con_obj) == [(1 + 2im) * x + 1]
    @test moi_set(con_obj) == MOI.Zeros(1)
end

function test_complex_vector_constraint_quadratic()
    model = Model()
    @variable(model, x)
    con_ref = @constraint(model, [(1 + 2im) * x^2 + 1] in MOI.Zeros(1))
    @test list_of_constraint_types(model) ==
          [(Vector{GenericQuadExpr{Complex{Float64},VariableRef}}, MOI.Zeros)]
    con_obj = constraint_object(con_ref)
    @test jump_function(con_obj) == [(1 + 2im) * x^2 + 1]
    @test moi_set(con_obj) == MOI.Zeros(1)
end

function test_complex_scalar_affine_constraint()
    model = Model()
    @variable(model, x)
    con_ref = @constraint(model, (1 + 2im) * x == 1.0)
    @test list_of_constraint_types(model) ==
          [(GenericAffExpr{ComplexF64,VariableRef}, MOI.EqualTo{ComplexF64})]
    con_obj = constraint_object(con_ref)
    @test jump_function(con_obj) == (1 + 2im) * x
    @test moi_set(con_obj) == MOI.EqualTo(1.0 + 0.0im)
    return
end

function test_complex_scalar_quadratic_constraint()
    model = Model()
    @variable(model, x)
    con_ref = @constraint(model, (1 + 2im) * x^2 == 1.0)
    @test list_of_constraint_types(model) ==
          [(GenericQuadExpr{ComplexF64,VariableRef}, MOI.EqualTo{ComplexF64})]
    con_obj = constraint_object(con_ref)
    @test jump_function(con_obj) == (1 + 2im) * x^2
    @test moi_set(con_obj) == MOI.EqualTo(1.0 + 0.0im)
    return
end

function test_complex_print()
    model = Model()
    @variable(model, x)
    y = (1 + 2im) * x + 1
    @test sprint(show, y) == "(1 + 2im) x + 1"
    y = im * x
    @test sprint(show, y) == "x im"
    return
end

function test_complex_print_zeros()
    model = Model()
    @variable(model, x in ComplexPlane())
    @test sprint(show, real(x)) == "real(x)"
    @test sprint(show, imag(x)) == "imag(x)"
    return
end

function test_complex_conj()
    model = Model()
    @variable(model, x)
    @test conj(x) === x
    @test real(x) === x
    @test imag(x) == 0
    real_aff = 2 * x + 3
    @test conj(real_aff) === real_aff
    @test real(real_aff) === real_aff
    @test imag(real_aff) == 0
    complex_aff = (2 + im) * x + 3 - im
    @test conj(complex_aff) == (2 - im) * x + 3 + im
    @test real(complex_aff) == 2x + 3
    @test imag(complex_aff) == x - 1
    real_quad = 4x^2 + 2 * x + 3
    @test conj(real_quad) === real_quad
    @test real(real_quad) === real_quad
    @test imag(real_quad) === real_quad
    complex_quad = (4 - 5im) * x^2 + (2 + im) * x + 3 - im
    @test conj(complex_quad) == (4 + 5im) * x^2 + (2 - im) * x + 3 + im
    @test real(complex_quad) == 4x^2 + 2x + 3
    @test imag(complex_quad) == -5x^2 + x - 1
end

function test_complex_abs2()
    model = Model()
    @variable(model, x)
    @test abs2(x) == x^2
    @test abs2(x + 2im) == x^2 + 4
    @test abs2(x + 2) == x^2 + 4x + 4
    @test abs2(x * im + 2) == x^2 + 4
end

function test_hermitian()
    model = Model()
    @variable(model, x)
    A = [3 1im; -1im 2x]
    @test A isa Matrix{GenericAffExpr{ComplexF64,VariableRef}}
    A = [3x^2 1im; -1im 2x]
    @test A isa Matrix{GenericQuadExpr{ComplexF64,VariableRef}}
    A = [3x 1im; -1im 2x^2]
    @test A isa Matrix{GenericQuadExpr{ComplexF64,VariableRef}}
    A = [3x 1im; -1im 2x]
    @test A isa Matrix{GenericAffExpr{ComplexF64,VariableRef}}
    @test isequal_canonical(A', A)
    H = LinearAlgebra.Hermitian(A)
    T = GenericAffExpr{ComplexF64,VariableRef}
    @test H isa LinearAlgebra.Hermitian{T,Matrix{T}}
    @test isequal_canonical(A[1, 2], LinearAlgebra.adjoint(A[2, 1]))
    @test isequal_canonical(H[1, 2], LinearAlgebra.adjoint(H[2, 1]))
    for i in 1:2, j in 1:2
        @test isequal_canonical(A[i, j], H[i, j])
    end
    return
end

function test_complex_sparse_arrays_dropzeros()
    model = Model()
    @variable(model, x)
    a = 2.0 + 1.0im
    for rhs in (0.0 + 0.0im, 0.0 - 0.0im, -0.0 + 0.0im, -0.0 + -0.0im)
        # We need to explicitly set the .constant field to avoid a conversion to
        # 0.0 + 0.0im
        expr = a * x
        expr.constant = rhs
        @test isequal(SparseArrays.dropzeros(expr), a * x)
        expr.constant = 1.0 + rhs
        @test isequal(SparseArrays.dropzeros(expr), a * x + 1.0)
    end
    return
end

function test_complex_hermitian_constraint()
    model = Model()
    @variable(model, x[1:2, 1:2])
    H = LinearAlgebra.Hermitian(x)
    @test vectorize(H, HermitianMatrixShape(2)) ==
          [x[1, 1], x[1, 2], x[2, 2], 0.0]
    @constraint(model, c, H in HermitianPSDCone())
    @test constraint_object(c).func == [x[1, 1], x[1, 2], x[2, 2], 0.0]
    @test function_string(MIME("text/plain"), constraint_object(c)) ==
          "[x[1,1]  x[1,2]\n x[1,2]  x[2,2]]"
    return
end

function test_complex_hermitian_inequality_constraint()
    model = Model()
    @variable(model, x[1:2, 1:2])
    H = LinearAlgebra.Hermitian(x)
    @test vectorize(H, HermitianMatrixShape(2)) ==
          [x[1, 1], x[1, 2], x[2, 2], 0.0]
    @constraint(model, c, H >= 0, HermitianPSDCone())
    @test constraint_object(c).func == [x[1, 1], x[1, 2], x[2, 2], 0.0]
    @test function_string(MIME("text/plain"), constraint_object(c)) ==
          "[x[1,1]  x[1,2]\n x[1,2]  x[2,2]]"
    @constraint(model, c2, 0 <= H, HermitianPSDCone())
    @test constraint_object(c2).func == [x[1, 1], x[1, 2], x[2, 2], 0.0]
    @test function_string(MIME("text/plain"), constraint_object(c2)) ==
          "[x[1,1]  x[1,2]\n x[1,2]  x[2,2]]"
    return
end

function test_isreal()
    model = Model()
    @variable(model, x[1:2])
    @test isreal(x[1])
    @test isreal(2 * x[1] + 3) == true
    @test isreal(2 * x[1] + 3 * x[2]) == true
    @test isreal(2 * x[1] + 3 * x[2] + 4) == true
    @test isreal(2 * x[1] + 0im) == true
    @test isreal(2 * x[1] + 1im) == false
    @test isreal(2im * x[1] + 1) == false
    @test isreal(x[1] * x[2] + 1) == true
    @test isreal(x[1] * x[2] + 1im) == false
    @test isreal(x[1] * x[2] * 1im + 2) == false
    return
end

function test_HermitianPSDCone_general_matrix_error()
    model = Model()
    @variable(model, X[1:2, 1:2] in HermitianPSDCone())
    @variable(model, t)
    Y = X + LinearAlgebra.I(2) * t
    @test LinearAlgebra.ishermitian(Y)
    @test !(Y isa LinearAlgebra.Hermitian)
    @test_throws_runtime(
        ErrorException(
            "In `@constraint(model, Y in HermitianPSDCone())`: " *
            "Unable to add matrix in HermitianPSDCone because the matrix is " *
            "not a subtype of `LinearAlgebra.Hermitian`. To fix, wrap the " *
            "matrix `H` in `LinearAlgebra.Hermitian(H)`.",
        ),
        @constraint(model, Y in HermitianPSDCone()),
    )
    return
end

function test_complex_generic_number_type()
    for T in (Float32, Float64, Rational{BigInt}, BigFloat)
        model = GenericModel{T}()
        @variable(model, x in ComplexPlane())
        @test x isa GenericAffExpr{Complex{T},GenericVariableRef{T}}
    end
    return
end

function test_hermitian_generic_number_type()
    for T in (Float32, Float64, Rational{BigInt}, BigFloat)
        model = GenericModel{T}()
        @variable(model, x[1:2, 1:2], Hermitian)
        @test x[1, 1] isa GenericAffExpr{Complex{T},GenericVariableRef{T}}
    end
    return
end

function test_hermitian_jump_scalar()
    model = Model()
    @variable(model, x[1:2, 1:2] in HermitianPSDCone())
    Y = im * x.data
    Z = [0 real(x[1, 2])*im-imag(x[1, 2]); -real(x[1, 2])*im-imag(x[1, 2]) 0]
    @test !isequal_canonical(Y, Z)
    @test isequal_canonical(
        LinearAlgebra.Hermitian(Y),
        LinearAlgebra.Hermitian(Z),
    )
    return
end

function test_mul_real_hermitian()
    A = [1+1im 1+1im; 1-2im 3+1im]
    model = Model()
    @variable(model, x)
    for s in (:L, :U), f in (x, x + 1, x^2)
        B = LinearAlgebra.Hermitian(A, s)
        @test f * B isa LinearAlgebra.Hermitian
        @test isequal_canonical(f * B, LinearAlgebra.Hermitian(f * A, s))
        @test B * f isa LinearAlgebra.Hermitian
        @test isequal_canonical(B * f, LinearAlgebra.Hermitian(A * f, s))
        @test isequal_canonical(f * B + f * B, (2 * f) * B)
    end
    return
end

end  # TestComplexNumberSupport

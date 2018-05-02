using JuMP
using Compat
using Compat.Test

@testset "Testing basic operator overloads" begin
    m = Model()
    @variable(m, w)
    @variable(m, x)
    @variable(m, y)
    @variable(m, z)

    aff = @inferred 7.1 * x + 2.5
    @test_expr_str 7.1 * x + 2.5 "7.1 x + 2.5"
    aff2 = @inferred 1.2 * y + 1.2
    @test_expr_str 1.2 * y + 1.2 "1.2 y + 1.2"
    q = @inferred 2.5 * y * z + aff
    @test_expr_str 2.5 * y * z + aff "2.5 y*z + 7.1 x + 2.5"
    q2 = @inferred 8 * x * z + aff2
    @test_expr_str 8 * x * z + aff2 "8 x*z + 1.2 y + 1.2"
    @test_expr_str 2 * x * x + 1 * y * y + z + 3 "2 x² + y² + z + 3"

    @test JuMP.isequal_canonical((@inferred 3w + 2y), @inferred 2y + 3w)
    @test !JuMP.isequal_canonical((@inferred 3w + 2y + 1), @inferred 3w + 2y)
    @test !JuMP.isequal_canonical((@inferred 3w + 2y), @inferred 3y + 2w)
    @test !JuMP.isequal_canonical((@inferred 3w + 2y), @inferred 3w + y)

    @test !JuMP.isequal_canonical(aff, aff2)
    @test !JuMP.isequal_canonical(aff2, aff)

    @test  JuMP.isequal_canonical(q, @inferred 2.5z*y + aff)
    @test !JuMP.isequal_canonical(q, @inferred 2.5y*z + aff2)
    @test !JuMP.isequal_canonical(q, @inferred 2.5x*z + aff)
    @test !JuMP.isequal_canonical(q, @inferred 2.5y*x + aff)
    @test !JuMP.isequal_canonical(q, @inferred 1.5y*z + aff)
    @test  JuMP.isequal_canonical(q2, @inferred 8z*x + aff2)
    @test !JuMP.isequal_canonical(q2, @inferred 8x*z + aff)
    @test !JuMP.isequal_canonical(q2, @inferred 7x*z + aff2)
    @test !JuMP.isequal_canonical(q2, @inferred 8x*y + aff2)
    @test !JuMP.isequal_canonical(q2, @inferred 8y*z + aff2)

    # Different objects that must all interact:
    # 1. Number
    # 2. Variable
    # 3. AffExpr
    # 4. QuadExpr

    # 1. Number tests
    @testset "Number--???" begin
        # 1-1 Number--Number - nope!
        # 1-2 Number--Variable
        @test_expr_str 4.13 + w "w + 4.13"
        @test_expr_str 3.16 - w "-w + 3.16"
        @test_expr_str 5.23 * w "5.23 w"
        @test_throws ErrorException 2.94 / w
        # 1-3 Number--AffExpr
        @test_expr_str 1.5 + aff "7.1 x + 4"
        @test_expr_str 1.5 - aff "-7.1 x - 1"
        @test_expr_str 2 * aff "14.2 x + 5"
        @test_throws ErrorException 2 / aff
        # 1-4 Number--QuadExpr
        @test_expr_str 1.5 + q "2.5 y*z + 7.1 x + 4"
        @test_expr_str 1.5 - q "-2.5 y*z - 7.1 x - 1"
        @test_expr_str 2 * q "5 y*z + 14.2 x + 5"
        @test_throws ErrorException 2 / q
    end

    # 2. Variable tests
    @testset "Variable--???" begin
        # 2-0 Variable unary
        @test (+x) === x
        @test_expr_str -x "-x"
        # 2-1 Variable--Number
        @test_expr_str w + 4.13 "w + 4.13"
        @test_expr_str w - 4.13 "w - 4.13"
        @test_expr_str w * 4.13 "4.13 w"
        @test_expr_str w / 2.00 "0.5 w"
        @test w == w
        @test_expr_str x*y - 1 "x*y - 1"
        # 2-2 Variable--Variable
        @test_expr_str w + x "w + x"
        @test_expr_str w - x "w - x"
        @test_expr_str w * x "w*x"
        @test_expr_str x - x "0"
        @test_throws ErrorException w / x
        @test_expr_str y*z - x "y*z - x"
        @test_expr_str x - x "0"
        # 2-3 Variable--AffExpr
        @test_expr_str z + aff "z + 7.1 x + 2.5"
        @test_expr_str z - aff "z - 7.1 x - 2.5"
        @test_expr_str z * aff "7.1 z*x + 2.5 z"
        @test_throws ErrorException z / aff
        @test_throws MethodError z ≤ aff
        @test_expr_str 7.1 * x - aff "-2.5"
        # 2-4 Variable--QuadExpr
        @test_expr_str w + q "2.5 y*z + w + 7.1 x + 2.5"
        @test_expr_str w - q "-2.5 y*z + w - 7.1 x - 2.5"
        @test_throws ErrorException w*q
        @test_throws ErrorException w/q
    end

    # 3. AffExpr tests
    @testset "AffExpr--???" begin
        # 3-0 AffExpr unary
        @test_expr_str +aff "7.1 x + 2.5"
        @test_expr_str -aff "-7.1 x - 2.5"
        # 3-1 AffExpr--Number
        @test_expr_str aff + 1.5 "7.1 x + 4"
        @test_expr_str aff - 1.5 "7.1 x + 1"
        @test_expr_str aff * 2 "14.2 x + 5"
        @test_expr_str aff / 2 "3.55 x + 1.25"
        @test_throws MethodError aff ≤ 1
        @test aff == aff
        @test_throws MethodError aff ≥ 1
        @test_expr_str aff - 1 "7.1 x + 1.5"
        # 3-2 AffExpr--Variable
        @test_expr_str aff + z "7.1 x + z + 2.5"
        @test_expr_str aff - z "7.1 x - z + 2.5"
        @test_expr_str aff * z "7.1 x*z + 2.5 z"
        @test_throws ErrorException aff/z
        @test_expr_str aff - 7.1 * x "2.5"
        # 3-3 AffExpr--AffExpr
        @test_expr_str aff + aff2 "7.1 x + 1.2 y + 3.7"
        @test_expr_str aff - aff2 "7.1 x - 1.2 y + 1.3"
        @test_expr_str aff * aff2 "8.52 x*y + 3 y + 8.52 x + 3"
        @test string((x+x)*(x+3)) == string((x+3)*(x+x))  # Issue #288
        @test_throws ErrorException aff/aff2
        @test_expr_str aff-aff "0"
        # 4-4 AffExpr--QuadExpr
        @test_expr_str aff2 + q "2.5 y*z + 1.2 y + 7.1 x + 3.7"
        @test_expr_str aff2 - q "-2.5 y*z + 1.2 y - 7.1 x - 1.3"
        @test_throws ErrorException aff2 * q
        @test_throws ErrorException aff2 / q
    end

    # 4. QuadExpr
    @testset "QuadExpr--???" begin
        # 4-0 QuadExpr unary
        @test_expr_str +q "2.5 y*z + 7.1 x + 2.5"
        @test_expr_str -q "-2.5 y*z - 7.1 x - 2.5"
        # 4-1 QuadExpr--Number
        @test_expr_str q + 1.5 "2.5 y*z + 7.1 x + 4"
        @test_expr_str q - 1.5 "2.5 y*z + 7.1 x + 1"
        @test_expr_str q * 2 "5 y*z + 14.2 x + 5"
        @test_expr_str q / 2 "1.25 y*z + 3.55 x + 1.25"
        @test q == q
        @test_expr_str aff2 - q "-2.5 y*z + 1.2 y - 7.1 x - 1.3"
        # 4-2 QuadExpr--Variable
        @test_expr_str q + w "2.5 y*z + w + 7.1 x + 2.5"
        @test_expr_str q - w "2.5 y*z + 7.1 x - w + 2.5"
        @test_throws ErrorException q*w
        @test_throws ErrorException q/w
        # 4-3 QuadExpr--AffExpr
        @test_expr_str q + aff2 "2.5 y*z + 7.1 x + 1.2 y + 3.7"
        @test_expr_str q - aff2 "2.5 y*z + 7.1 x - 1.2 y + 1.3"
        @test_throws ErrorException q * aff2
        @test_throws ErrorException q / aff2
        # 4-4 QuadExpr--QuadExpr
        @test_expr_str q + q2 "2.5 y*z + 8 x*z + 7.1 x + 1.2 y + 3.7"
        @test_expr_str q - q2 "2.5 y*z - 8 x*z + 7.1 x - 1.2 y + 1.3"
        @test_throws ErrorException q * q2
        @test_throws ErrorException q / q2
    end
end

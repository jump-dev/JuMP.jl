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
    @test_expression_with_string 7.1 * x + 2.5 "7.1 x + 2.5"
    aff2 = @inferred 1.2 * y + 1.2
    @test_expression_with_string 1.2 * y + 1.2 "1.2 y + 1.2"
    q = @inferred 2.5 * y * z + aff
    @test_expression_with_string 2.5 * y * z + aff "2.5 y*z + 7.1 x + 2.5"
    q2 = @inferred 8 * x * z + aff2
    @test_expression_with_string 8 * x * z + aff2 "8 x*z + 1.2 y + 1.2"
    @test_expression_with_string 2 * x * x + 1 * y * y + z + 3 "2 x² + y² + z + 3"

    @testset "isequal_canonical" begin
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
    end

    # Different objects that must all interact:
    # 1. Number
    # 2. Variable
    # 3. AffExpr
    # 4. QuadExpr

    # 1. Number tests
    @testset "Number--???" begin
        # 1-1 Number--Number - nope!
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

    # 2. Variable tests
    @testset "Variable--???" begin
        # 2-0 Variable unary
        @test (+x) === x
        @test_expression_with_string -x "-x"
        # 2-1 Variable--Number
        @test_expression_with_string w + 4.13 "w + 4.13"
        @test_expression_with_string w - 4.13 "w - 4.13"
        @test_expression_with_string w * 4.13 "4.13 w"
        @test_expression_with_string w / 2.00 "0.5 w"
        @test w == w
        @test_expression_with_string x*y - 1 "x*y - 1"
        # 2-2 Variable--Variable
        @test_expression_with_string w + x "w + x"
        @test_expression_with_string w - x "w - x"
        @test_expression_with_string w * x "w*x"
        @test_expression_with_string x - x "0"
        @test_throws ErrorException w / x
        @test_expression_with_string y*z - x "y*z - x"
        @test_expression_with_string x - x "0"
        # 2-3 Variable--AffExpr
        @test_expression_with_string z + aff "z + 7.1 x + 2.5"
        @test_expression_with_string z - aff "z - 7.1 x - 2.5"
        @test_expression_with_string z * aff "7.1 z*x + 2.5 z"
        @test_throws ErrorException z / aff
        @test_throws MethodError z ≤ aff
        @test_expression_with_string 7.1 * x - aff "-2.5"
        # 2-4 Variable--QuadExpr
        @test_expression_with_string w + q "2.5 y*z + w + 7.1 x + 2.5"
        @test_expression_with_string w - q "-2.5 y*z + w - 7.1 x - 2.5"
        @test_throws ErrorException w*q
        @test_throws ErrorException w/q
    end

    # 3. AffExpr tests
    @testset "AffExpr--???" begin
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
        # 3-2 AffExpr--Variable
        @test_expression_with_string aff + z "7.1 x + z + 2.5"
        @test_expression_with_string aff - z "7.1 x - z + 2.5"
        @test_expression_with_string aff * z "7.1 x*z + 2.5 z"
        @test_throws ErrorException aff/z
        @test_expression_with_string aff - 7.1 * x "2.5"
        # 3-3 AffExpr--AffExpr
        @test_expression_with_string aff + aff2 "7.1 x + 1.2 y + 3.7"
        @test_expression_with_string aff - aff2 "7.1 x - 1.2 y + 1.3"
        @test_expression_with_string aff * aff2 "8.52 x*y + 3 y + 8.52 x + 3"
        @test string((x+x)*(x+3)) == string((x+3)*(x+x))  # Issue #288
        @test_throws ErrorException aff/aff2
        @test_expression_with_string aff-aff "0"
        # 4-4 AffExpr--QuadExpr
        @test_expression_with_string aff2 + q "2.5 y*z + 1.2 y + 7.1 x + 3.7"
        @test_expression_with_string aff2 - q "-2.5 y*z + 1.2 y - 7.1 x - 1.3"
        @test_throws ErrorException aff2 * q
        @test_throws ErrorException aff2 / q
    end

    # 4. QuadExpr
    @testset "QuadExpr--???" begin
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
        @test_expression_with_string q + w "2.5 y*z + w + 7.1 x + 2.5"
        @test_expression_with_string q - w "2.5 y*z + 7.1 x - w + 2.5"
        @test_throws ErrorException q*w
        @test_throws ErrorException q/w
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
    end

    @testset "Higher-level operators" begin
        @testset "sum" begin
            sum_m = Model()
            @variable(sum_m, 0 ≤ matrix[1:3,1:3] ≤ 1, start = 1)
            @testset "sum(j::JuMPArray{Variable})" begin
                @test_expression_with_string sum(matrix) "matrix[1,1] + matrix[2,1] + matrix[3,1] + matrix[1,2] + matrix[2,2] + matrix[3,2] + matrix[1,3] + matrix[2,3] + matrix[3,3]"
            end

            @testset "sum(j::JuMPArray{T}) where T<:Real" begin
                @test sum(JuMP.startvalue.(matrix)) ≈ 9
            end
            @testset "sum(j::Array{VariableRef})" begin
                @test string(sum(matrix[1:3,1:3])) == string(sum(matrix))
            end
            @testset "sum(affs::Array{AffExpr})" begin
                @test_expression_with_string sum([2*matrix[i,j] for i in 1:3, j in 1:3]) "2 matrix[1,1] + 2 matrix[2,1] + 2 matrix[3,1] + 2 matrix[1,2] + 2 matrix[2,2] + 2 matrix[3,2] + 2 matrix[1,3] + 2 matrix[2,3] + 2 matrix[3,3]"
            end

            S = [1,3]
            @variable(sum_m, x[S], start=1)
            @testset "sum(j::JuMPDict{VariableRef})" begin
                @test_expression sum(x)
                @test length(string(sum(x))) == 11 # order depends on hashing
                @test contains(string(sum(x)),"x[1]")
                @test contains(string(sum(x)),"x[3]")
            end
            @testset "sum(j::JuMPDict{T}) where T<:Real" begin
                @test sum(JuMP.startvalue.(x)) == 2
            end
        end

        @testset "dot" begin
            dot_m = Model()
            @variable(dot_m, 0 ≤ x[1:3] ≤ 1)
            c = vcat(1:3)
            @test_expression_with_string dot(c,x) "x[1] + 2 x[2] + 3 x[3]"
            @test_expression_with_string dot(x,c) "x[1] + 2 x[2] + 3 x[3]"

            A = [1 3 ; 2 4]
            @variable(dot_m, 1 ≤ y[1:2,1:2] ≤ 1)
            @test_expression_with_string vecdot(A,y) "y[1,1] + 2 y[2,1] + 3 y[1,2] + 4 y[2,2]"
            @test_expression_with_string vecdot(y,A) "y[1,1] + 2 y[2,1] + 3 y[1,2] + 4 y[2,2]"

            B = ones(2,2,2)
            @variable(dot_m, 0 ≤ z[1:2,1:2,1:2] ≤ 1)
            @test_expression_with_string vecdot(B,z) "z[1,1,1] + z[2,1,1] + z[1,2,1] + z[2,2,1] + z[1,1,2] + z[2,1,2] + z[1,2,2] + z[2,2,2]"
            @test_expression_with_string vecdot(z,B) "z[1,1,1] + z[2,1,1] + z[1,2,1] + z[2,2,1] + z[1,1,2] + z[2,1,2] + z[1,2,2] + z[2,2,2]"

            @objective(dot_m, Max, dot(x, ones(3)) - vecdot(y, ones(2,2)))
            for i in 1:3
                JuMP.setstartvalue(x[i], 1)
            end
            for i in 1:2, j in 1:2
                JuMP.setstartvalue(y[i,j], 1)
            end
            for i in 1:2, j in 1:2, k in 1:2
                JuMP.setstartvalue(z[i,j,k], 1)
            end
            @test dot(c, JuMP.startvalue.(x)) ≈ 6
            @test vecdot(A, JuMP.startvalue.(y)) ≈ 10
            @test vecdot(B, JuMP.startvalue.(z)) ≈ 8

            @testset "JuMP issue #656" begin
                issue656 = Model()
                @variable(issue656, x)
                floats = Float64[i for i in 1:2]
                anys   = Array{Any}(undef, 2)
                anys[1] = 10
                anys[2] = 20 + x
                @test dot(floats, anys) == 10 + 40 + 2x
            end

            @testset "JuMP PR #943" begin
                pull943 = Model()
                @variable(pull943, x[1 : 10^6]);
                JuMP.setstartvalue.(x, 1 : 10^6)
                @expression(pull943, testsum, sum(x[i] * i for i = 1 : 10^6))
                @expression(pull943, testdot1, dot(x, 1 : 10^6))
                @expression(pull943, testdot2, dot(1 : 10^6, x))
                @test JuMP.value(testsum, JuMP.startvalue) ≈ JuMP.value(testdot1, JuMP.startvalue)
                @test JuMP.value(testsum, JuMP.startvalue) ≈ JuMP.value(testdot2, JuMP.startvalue)
            end
        end
    end
end

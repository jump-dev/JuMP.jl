# For "expression^3 and unary*"
struct PowVariable <: JuMP.AbstractVariableRef
    pow::Int
end
Base.:^(x::PowVariable, i::Int) = PowVariable(x.pow*i)
Base.:*(x::PowVariable, y::PowVariable) = PowVariable(x.pow + y.pow)
Base.copy(x::PowVariable) = x

function expressions_test(ModelType::Type{<:JuMP.AbstractModel}, VariableRefType::Type{<:JuMP.AbstractVariableRef})
    AffExprType = JuMP.GenericAffExpr{Float64, VariableRefType}
    QuadExprType = JuMP.GenericQuadExpr{Float64, VariableRefType}

    @testset "isequal(::GenericAffExpr)" begin
        m = Model()
        @variable(m, x)
        @test isequal(x + 1, x + 1)
    end

    @testset "hash(::GenericAffExpr)" begin
        m = Model()
        @variable(m, x)
        @test hash(x + 1) == hash(x + 1)
    end

    @testset "isequal(::GenericQuadExpr)" begin
        m = Model()
        @variable(m, x)
        @test isequal(x^2 + 1, x^2 + 1)
    end

    @testset "hash(::GenericQuadExpr)" begin
        m = Model()
        @variable(m, x)
        @test hash(x^2 + 1) == hash(x^2 + 1)
    end

    @testset "value for GenericAffExpr" begin
        expr1 = JuMP.GenericAffExpr(3.0, 3 => -5.0, 2 => 4.0)
        @test @inferred(JuMP.value(expr1, -)) == 10.0
        expr2 = JuMP.GenericAffExpr{Int,Int}(2)
        @test typeof(@inferred(JuMP.value(expr2, i -> 1.0))) == Float64
        @test @inferred(JuMP.value(expr2, i -> 1.0)) == 2.0
    end

    @testset "add_to_expression!(::GenericAffExpr{C,V}, ::V)" begin
        aff = JuMP.GenericAffExpr(1.0, :a => 2.0)
        @test JuMP.isequal_canonical(JuMP.add_to_expression!(aff, :b),
                                     JuMP.GenericAffExpr(1.0, :a => 2.0, :b => 1.0))
    end

    @testset "add_to_expression!(::GenericAffExpr{C,V}, ::C)" begin
        aff = JuMP.GenericAffExpr(1.0, :a => 2.0)
        @test JuMP.isequal_canonical(JuMP.add_to_expression!(aff, 1.0),
                                     JuMP.GenericAffExpr(2.0, :a => 2.0))
    end

    @testset "linearterms(::AffExpr)" begin
        m = Model()
        @variable(m, x[1:10])

        aff = 1*x[1] + 2*x[2]
        k = 0
        @test length(linearterms(aff)) == 2
        for (coeff, var) in linearterms(aff)
            if k == 0
                @test coeff == 1
                @test var === x[1]
            elseif k == 1
                @test coeff == 2
                @test var === x[2]
            end
            k += 1
        end
        @test k == 2
    end

    @testset "linearterms(::AffExpr) for empty expression" begin
        k = 0
        aff = zero(AffExprType)
        @test length(linearterms(aff)) == 0
        for (coeff, var) in linearterms(aff)
            k += 1
        end
        @test k == 0
    end

    @testset "Copy AffExpr between models" begin
        m = Model()
        @variable(m, x)
        m2 = Model()
        aff = copy(2x + 1, m2)
        aff_expected = 2*copy(x, m2) + 1
        @test JuMP.isequal_canonical(aff, aff_expected)
    end

    @testset "destructive_add!(ex::Number, c::Number, x::GenericAffExpr)" begin
        aff = JuMP.destructive_add!(1.0, 2.0, JuMP.GenericAffExpr(1.0, :a => 1.0))
        @test JuMP.isequal_canonical(aff, JuMP.GenericAffExpr(3.0, :a => 2.0))
    end

    @testset "destructive_add!(ex::Number, c::Number, x::GenericQuadExpr) with c == 0" begin
        quad = JuMP.destructive_add!(2.0, 0.0, QuadExprType())
        @test JuMP.isequal_canonical(quad, convert(QuadExprType, 2.0))
    end

    @testset "destructive_add!(ex::Number, c::VariableRef, x::VariableRef)" begin
        m = Model()
        @variable(m, x)
        @variable(m, y)
        @test_expression_with_string JuMP.destructive_add!(5.0, x, y) "x*y + 5"
    end

    @testset "destructive_add!(ex::Number, c::T, x::T) where T<:GenericAffExpr" begin
        m = Model()
        @variable(m, x)
        @test_expression_with_string JuMP.destructive_add!(1.0, 2x, x+1) "2 x² + 2 x + 1"
    end

    @testset "destructive_add!(ex::Number, c::GenericAffExpr{C,V}, x::V) where {C,V}" begin
        m = Model()
        @variable(m, x)
        @test_expression_with_string JuMP.destructive_add!(1.0, 2x, x) "2 x² + 1"
    end

    @testset "destructive_add!(ex::Number, c::GenericQuadExpr, x::Number)" begin
        m = Model()
        @variable(m, x)
        @test_expression_with_string JuMP.destructive_add!(0.0, x^2, 1.0) "x²"
    end

    @testset "destructive_add!(ex::Number, c::GenericQuadExpr, x::Number) with c == 0" begin
        m = Model()
        @variable(m, x)
        @test_expression_with_string JuMP.destructive_add!(0.0, x^2, 0.0) "0"
    end

    @testset "destructive_add!(aff::AffExpr,c::VariableRef,x::AffExpr)" begin
        m = Model()
        @variable(m, x)
        @test_expression_with_string JuMP.destructive_add!(2x, x, x + 1) "x² + 3 x"
    end

    @testset "destructive_add!(aff::GenericAffExpr{C,V},c::GenericAffExpr{C,V},x::Number) where {C,V}" begin
        m = Model()
        @variable(m, x)
        @test_expression_with_string JuMP.destructive_add!(2x, x, 1) "3 x"
    end

    @testset "destructive_add!(aff::GenericAffExpr{C,V}, c::GenericQuadExpr{C,V}, x::Number) where {C,V}" begin
        m = Model()
        @variable(m, x)
        @test_expression_with_string JuMP.destructive_add!(2x, x^2, 1) "x² + 2 x"
    end

    @testset "destructive_add!(aff::GenericAffExpr{C,V}, c::GenericQuadExpr{C,V}, x::Number) where {C,V} with x == 0" begin
        m = Model()
        @variable(m, x)
        @test_expression_with_string JuMP.destructive_add!(2x, x^2, 0) "2 x"
    end

    @testset "destructive_add!(aff::GenericAffExpr{C,V}, c::Number, x::GenericQuadExpr{C,V}) where {C,V} with c == 0" begin
        m = Model()
        @variable(m, x)
        @test_expression_with_string JuMP.destructive_add!(2x, 0, x^2) "2 x"
    end

    @testset "destructive_add!(ex::GenericAffExpr{C,V}, c::GenericAffExpr{C,V}, x::GenericAffExpr{C,V}) where {C,V}" begin
        m = Model()
        @variable(m, x)
        @test_expression_with_string JuMP.destructive_add!(2x, x + 1, x + 0) "x² + 3 x"
    end

    @testset "destructive_add!(quad::GenericQuadExpr{C,V},c::GenericAffExpr{C,V},x::Number) where {C,V}" begin
        m = Model()
        @variable(m, x)
        @test_expression_with_string JuMP.destructive_add!(x^2, x + 1, 1) "x² + x + 1"
    end

    @testset "destructive_add!(quad::GenericQuadExpr{C,V},c::V,x::GenericAffExpr{C,V}) where {C,V}" begin
        m = Model()
        @variable(m, x)
        @test_expression_with_string JuMP.destructive_add!(x^2, x, x+1) "2 x² + x"
    end

    @testset "destructive_add!(quad::GenericQuadExpr{C,V},c::GenericQuadExpr{C,V},x::Number) where {C,V}" begin
        m = Model()
        @variable(m, x)
        @test_expression_with_string JuMP.destructive_add!(x^2 + x, x^2 + x, 2.0) "3 x² + 3 x"
    end

    @testset "destructive_add!(ex::GenericQuadExpr{C,V}, c::GenericAffExpr{C,V}, x::GenericAffExpr{C,V}) where {C,V}" begin
        m = Model()
        @variable(m, x)
        @test_expression_with_string JuMP.destructive_add!(x^2 + x, x + 0, x + 1) "2 x² + 2 x"
    end

    @testset "(+)(::AffExpr)" begin
        m = Model()
        @variable(m, x)
        @test_expression_with_string (+)(x + 1) "x + 1"
    end

    @testset "(+)(::QuadExpr)" begin
        m = Model()
        @variable(m, x)
        @test_expression_with_string (+)(x^2 + 1) "x² + 1"
    end

    @testset "sum(::Vector{VariableRef})" begin
        m = Model()
        @variable(m, x[1:2])
        @test_expression_with_string sum(x) "x[1] + x[2]"
    end

    @testset "expression^3 and unary*" begin
        m = Model()
        x = PowVariable(1)
        # Calls (*)((x*x)^6)
        y = @expression m (x*x)^3
        @test y.pow == 6
        z = @inferred (x*x)^3
        @test z.pow == 6
    end
end

@testset "Expressions for JuMP.Model" begin
    expressions_test(Model, VariableRef)
end

@testset "Expressions for JuMPExtension.MyModel" begin
    expressions_test(JuMPExtension.MyModel, JuMPExtension.MyVariableRef)
end

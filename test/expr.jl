#  Copyright 2017, Iain Dunning, Joey Huchette, Miles Lubin, and contributors
#  This Source Code Form is subject to the terms of the Mozilla Public
#  License, v. 2.0. If a copy of the MPL was not distributed with this
#  file, You can obtain one at http://mozilla.org/MPL/2.0/.
#############################################################################
# JuMP
# An algebraic modeling language for Julia
# See http://github.com/JuliaOpt/JuMP.jl
#############################################################################
# test/expr.jl
#############################################################################

module ExpressionTests

using JuMP
using Pukeko
import Test: @inferred

include("utilities.jl")
include("JuMPExtension.jl")
const BOTH_VAR = [VariableRef, JuMPExtension.MyVariableRef]

function test_isequal_hash_affexpr()
    m = Model()
    @variable(m, x)
    @test isequal(x + 1, x + 1)
    @test hash(x + 1) == hash(x + 1)
end

function test_isequal_hash_quadexpr()
    m = Model()
    @variable(m, x)
    @test isequal(x^2 + 1, x^2 + 1)
    @test hash(x^2 + 1) == hash(x^2 + 1)
end

function test_value_affexpr()
    expr1 = JuMP.GenericAffExpr(3.0, 3 => -5.0, 2 => 4.0)
    @test @inferred(JuMP.value(expr1, -)) == 10.0
    expr2 = JuMP.GenericAffExpr{Int,Int}(2)
    @test typeof(@inferred(JuMP.value(expr2, i -> 1.0))) == Float64
    @test @inferred(JuMP.value(expr2, i -> 1.0)) == 2.0
end

function test_value_quadexpr()
    # 1 + 2x(1) + 3x(2)
    affine_term = JuMP.GenericAffExpr(1.0, 1 => 2.0, 2 => 3.0)
    # 1 + 2x(1) + 3x(2) + 4x(1)^2 + 5x(1)*x(2) + 6x(2)^2
    expr = JuMP.GenericQuadExpr(affine_term,
        JuMP.UnorderedPair(1, 1) => 4.0,
        JuMP.UnorderedPair(1, 2) => 5.0,
        JuMP.UnorderedPair(2, 2) => 6.0)
    @test typeof(@inferred(JuMP.value(expr, i -> 1.0))) == Float64
    @test @inferred(JuMP.value(expr, i -> 1.0)) == 21
    @test @inferred(JuMP.value(expr, i -> 2.0)) == 71
end

function test_add_to_expr_by_var()
    # add_to_expression!(::GenericAffExpr{C,V}, ::V)
    aff = JuMP.GenericAffExpr(1.0, :a => 2.0)
    @test JuMP.isequal_canonical(JuMP.add_to_expression!(aff, :b),
                                    JuMP.GenericAffExpr(1.0, :a => 2.0, :b => 1.0))
end

function test_add_to_expr_coef()
    # add_to_expression!(::GenericAffExpr{C,V}, ::C)
    aff = JuMP.GenericAffExpr(1.0, :a => 2.0)
    @test JuMP.isequal_canonical(JuMP.add_to_expression!(aff, 1.0),
                                    JuMP.GenericAffExpr(2.0, :a => 2.0))
end

function test_linear_terms_affexpr()
    m = Model()
    @variable(m, x[1:10])
    aff = 1*x[1] + 2*x[2]
    k = 0
    @test length(linear_terms(aff)) == 2
    for (coeff, var) in linear_terms(aff)
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

function linear_terms_affexpr_empty(VariableRefType)
    AffExprType = JuMP.GenericAffExpr{Float64, VariableRefType}
    k = 0
    aff = zero(AffExprType)
    @test length(linear_terms(aff)) == 0
    for (coeff, var) in linear_terms(aff)
        k += 1
    end
    @test k == 0
end
@parametric linear_terms_affexpr_empty BOTH_VAR

function test_copy_affexpr_between_models()
    m = Model()
    @variable(m, x)
    m2 = Model()
    aff = copy(2x + 1, m2)
    aff_expected = 2*copy(x, m2) + 1
    @test JuMP.isequal_canonical(aff, aff_expected)
end

function test_destructive_add_affexpr()
    # destructive_add!(ex::Number, c::Number, x::GenericAffExpr)
    aff = JuMP.destructive_add!(1.0, 2.0, JuMP.GenericAffExpr(1.0, :a => 1.0))
    @test JuMP.isequal_canonical(aff, JuMP.GenericAffExpr(3.0, :a => 2.0))
end

function destructive_add_quadexpr(VariableRefType)
    QuadExprType = JuMP.GenericQuadExpr{Float64, VariableRefType}
    # destructive_add!(ex::Number, c::Number, x::GenericQuadExpr) with c == 0
    quad = JuMP.destructive_add!(2.0, 0.0, QuadExprType())
    @test JuMP.isequal_canonical(quad, convert(QuadExprType, 2.0))
end
@parametric destructive_add_quadexpr BOTH_VAR

function test_destructive_add_products()
    m = Model()
    @variable(m, x)
    @variable(m, y)
    # destructive_add!(ex::Number, c::VariableRef, x::VariableRef)
    @test_expression_with_string JuMP.destructive_add!(5.0, x, y) "x*y + 5"
    # destructive_add!(ex::Number, c::T, x::T) where T<:GenericAffExpr
    @test_expression_with_string JuMP.destructive_add!(1.0, 2x, x+1) "2 x² + 2 x + 1"
    # destructive_add!(ex::Number, c::GenericAffExpr{C,V}, x::V) where {C,V}
    @test_expression_with_string JuMP.destructive_add!(1.0, 2x, x) "2 x² + 1"
    # destructive_add!(ex::Number, c::GenericQuadExpr, x::Number)
    @test_expression_with_string JuMP.destructive_add!(0.0, x^2, 1.0) "x²"
    # destructive_add!(ex::Number, c::GenericQuadExpr, x::Number) with c == 0
    @test_expression_with_string JuMP.destructive_add!(0.0, x^2, 0.0) "0"
    # destructive_add!(aff::AffExpr,c::VariableRef,x::AffExpr)
    @test_expression_with_string JuMP.destructive_add!(2x, x, x + 1) "x² + 3 x"
    # destructive_add!(aff::GenericAffExpr{C,V},c::GenericAffExpr{C,V},x::Number) where {C,V}
    @test_expression_with_string JuMP.destructive_add!(2x, x, 1) "3 x"
    # destructive_add!(aff::GenericAffExpr{C,V}, c::GenericQuadExpr{C,V}, x::Number) where {C,V}
    @test_expression_with_string JuMP.destructive_add!(2x, x^2, 1) "x² + 2 x"
    # destructive_add!(aff::GenericAffExpr{C,V}, c::GenericQuadExpr{C,V}, x::Number) where {C,V} with x == 0
    @test_expression_with_string JuMP.destructive_add!(2x, x^2, 0) "2 x"
    # destructive_add!(aff::GenericAffExpr{C,V}, c::Number, x::GenericQuadExpr{C,V}) where {C,V} with c == 0
    @test_expression_with_string JuMP.destructive_add!(2x, 0, x^2) "2 x"
    # destructive_add!(ex::GenericAffExpr{C,V}, c::GenericAffExpr{C,V}, x::GenericAffExpr{C,V}) where {C,V}
    @test_expression_with_string JuMP.destructive_add!(2x, x + 1, x + 0) "x² + 3 x"
    # destructive_add!(quad::GenericQuadExpr{C,V},c::GenericAffExpr{C,V},x::Number) where {C,V}
    @test_expression_with_string JuMP.destructive_add!(x^2, x + 1, 1) "x² + x + 1"
    # destructive_add!(quad::GenericQuadExpr{C,V},c::V,x::GenericAffExpr{C,V}) where {C,V}
    @test_expression_with_string JuMP.destructive_add!(x^2, x, x+1) "2 x² + x"
    # destructive_add!(quad::GenericQuadExpr{C,V},c::GenericQuadExpr{C,V},x::Number) where {C,V}
    @test_expression_with_string JuMP.destructive_add!(x^2 + x, x^2 + x, 2.0) "3 x² + 3 x"
    # destructive_add!(ex::GenericQuadExpr{C,V}, c::GenericAffExpr{C,V}, x::GenericAffExpr{C,V}) where {C,V}
    @test_expression_with_string JuMP.destructive_add!(x^2 + x, x + 0, x + 1) "2 x² + 2 x"
end

function test_singleton_add_affexpr()
    m = Model()
    @variable(m, x)
    @test_expression_with_string (+)(x + 1) "x + 1"
end

function test_singleton_add_quadexpr()
    m = Model()
    @variable(m, x)
    @test_expression_with_string (+)(x^2 + 1) "x² + 1"
end

function test_sum_vector_variableref()
    m = Model()
    @variable(m, x[1:2])
    @test_expression_with_string sum(x) "x[1] + x[2]"  # sum(::Vector{VR})
end


# For "expression^3 and unary*"
struct PowVariable <: JuMP.AbstractVariableRef
    pow::Int
end
Base.:^(x::PowVariable, i::Int) = PowVariable(x.pow*i)
Base.:*(x::PowVariable, y::PowVariable) = PowVariable(x.pow + y.pow)
Base.copy(x::PowVariable) = x

function test_expression_cubed_and_unary_mult()
    m = Model()
    x = PowVariable(1)
    # Calls (*)((x*x)^6)
    y = @expression m (x*x)^3
    @test y.pow == 6
    z = @inferred (x*x)^3
    @test z.pow == 6
end

end  # module ExpressionTests

import Pukeko
Pukeko.run_tests(ExpressionTests)
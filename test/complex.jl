module TestComplexNumberSupport

using JuMP
using Test
import MutableArithmetics
const MA = MutableArithmetics

function runtests()
    for name in names(@__MODULE__; all = true)
        if startswith("$(name)", "test_")
            @testset "$(name)" begin
                getfield(@__MODULE__, name)()
            end
        end
    end
    return
end

function test_complex_aff_expr()
    model = Model()
    @variable(model, x)
    y = (1 + 2im) * x + 1
    @test typeof(y) == GenericAffExpr{Complex{Float64},VariableRef}
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

function test_complex_scalar_constraint()
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

function test_complex_print()
    model = Model()
    @variable(model, x)
    y = (1 + 2im) * x + 1
    @test sprint(show, y) == "(1.0 + 2.0im) x + (1.0 + 0.0im)"
    y = im * x
    @test sprint(show, y) == "(0.0 + 1.0im) x"
    return
end

function test_complex_conj()
    model = Model()
    @variable(model, x)
    @test conj(x) === x
    @test real(x) === x
    @test imag(x) === x
    real_aff = 2 * x + 3
    @test conj(real_aff) === real_aff
    @test real(real_aff) === real_aff
    @test imag(real_aff) === real_aff
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

end

TestComplexNumberSupport.runtests()

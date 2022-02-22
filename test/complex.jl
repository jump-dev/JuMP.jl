module TestComplexNumberSupport

using JuMP
using Test

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
    y = x + (1 + 2im)
    @test typeof(y) == GenericAffExpr{Complex{Float64},VariableRef}
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

function test_complex_constraint()
    model = Model()
    @variable(model, x)
    @constraint(model, [(1 + 2im) * x + 1] in MOI.Zeros(1))
    @test list_of_constraint_types(model) ==
          [(Vector{GenericAffExpr{Complex{Float64},VariableRef}}, MOI.Zeros)]
    return
end

function test_complex_print()
    model = Model()
    @variable(model, x)
    y = (1 + 2im) * x + 1
    @test sprint(show, y) == "(1.0 + 2.0im) x + (1.0 + 0.0im)"
    return
end

end

TestComplexNumberSupport.runtests()

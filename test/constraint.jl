function constraints_test(ModelType::Type{<:JuMP.AbstractModel})
    @testset "SingleVariable constraints" begin
        m = ModelType()
        @variable(m, x)

        # x <= 10.0 doesn't translate to a SingleVariable constraint because
        # the LHS is first subtracted to form x - 10.0 <= 0.
        @constraint(m, cref, x in MOI.LessThan(10.0))
        @test JuMP.name(cref) == "cref"
        c = JuMP.constraint_object(cref)
        @test c.func == x
        @test c.set == MOI.LessThan(10.0)

        @variable(m, y[1:2])
        @constraint(m, cref2[i=1:2], y[i] in MOI.LessThan(float(i)))
        @test JuMP.name(cref2[1]) == "cref2[1]"
        c = JuMP.constraint_object(cref2[1])
        @test c.func == y[1]
        @test c.set == MOI.LessThan(1.0)
    end

    @testset "VectorOfVariables constraints" begin
        m = ModelType()
        @variable(m, x[1:2])

        cref = @constraint(m, x in MOI.Zeros(2))
        c = JuMP.constraint_object(cref)
        @test c.func == x
        @test c.set == MOI.Zeros(2)

        cref = @constraint(m, [x[2],x[1]] in MOI.Zeros(2))
        c = JuMP.constraint_object(cref)
        @test c.func == [x[2],x[1]]
        @test c.set == MOI.Zeros(2)
    end

    @testset "AffExpr constraints" begin
        m = ModelType()
        @variable(m, x)

        cref = @constraint(m, 2x <= 10)
        @test JuMP.name(cref) == ""
        JuMP.set_name(cref, "c")
        @test JuMP.name(cref) == "c"

        c = JuMP.constraint_object(cref)
        @test JuMP.isequal_canonical(c.func, 2x)
        @test c.set == MOI.LessThan(10.0)

        cref = @constraint(m, 3x + 1 ≥ 10)
        c = JuMP.constraint_object(cref)
        @test JuMP.isequal_canonical(c.func, 3x)
        @test c.set == MOI.GreaterThan(9.0)

        cref = @constraint(m, 1 == -x)
        c = JuMP.constraint_object(cref)
        @test JuMP.isequal_canonical(c.func, 1.0x)
        @test c.set == MOI.EqualTo(-1.0)

        @test_throws ErrorException @constraint(m, [x, 2x] == [1-x, 3])
        @test_macro_throws ErrorException @constraint(m, [x == 1-x, 2x == 3])
        cref = @constraint(m, [x, 2x] .== [1-x, 3])
        c = JuMP.constraint_object.(cref)
        @test JuMP.isequal_canonical(c[1].func, 2.0x)
        @test c[1].set == MOI.EqualTo(1.0)
        @test JuMP.isequal_canonical(c[2].func, 2.0x)
        @test c[2].set == MOI.EqualTo(3.0)

        @test MOI.isvalid(m, cref[1])
        MOI.delete!(m, cref[1])
        @test !MOI.isvalid(m, cref[1])
    end

    @testset "Two-sided constraints" begin
        m = ModelType()
        @variable(m, x)
        @variable(m, y)

        @constraint(m, cref, 1.0 <= x + y + 1.0 <= 2.0)
        @test JuMP.name(cref) == "cref"

        c = JuMP.constraint_object(cref)
        @test JuMP.isequal_canonical(c.func, x + y)
        @test c.set == MOI.Interval(0.0, 1.0)
    end

    @testset "Broadcasted constraint (.==)" begin
        m = ModelType()
        @variable(m, x[1:2])

        A = [1.0 2.0; 3.0 4.0]
        b = [4.0, 5.0]

        cref = @constraint(m, A*x .== b)
        @test size(cref) == (2,)

        c1 = JuMP.constraint_object(cref[1])
        @test JuMP.isequal_canonical(c1.func, x[1] + 2x[2])
        @test c1.set == MOI.EqualTo(4.0)
        c2 = JuMP.constraint_object(cref[2])
        @test JuMP.isequal_canonical(c2.func, 3x[1] + 4x[2])
        @test c2.set == MOI.EqualTo(5.0)
    end

    @testset "Broadcasted constraint (.<=)" begin
        m = ModelType()
        @variable(m, x[1:2,1:2])

        UB = [1.0 2.0; 3.0 4.0]

        cref = @constraint(m, x + 1 .<= UB)
        @test size(cref) == (2,2)
        for i in 1:2
            for j in 1:2
                c = JuMP.constraint_object(cref[i,j])
                @test JuMP.isequal_canonical(c.func, x[i,j] + 0)
                @test c.set == MOI.LessThan(UB[i,j] - 1)
            end
        end
    end

    @testset "Broadcasted two-sided constraint" begin
        m = ModelType()
        @variable(m, x[1:2])
        @variable(m, y[1:2])
        l = [1.0, 2.0]
        u = [3.0, 4.0]

        cref = @constraint(m, l .<= x + y + 1 .<= u)
        @test size(cref) == (2,)

        for i in 1:2
            c = JuMP.constraint_object(cref[i])
            @test JuMP.isequal_canonical(c.func, x[i] + y[i])
            @test c.set == MOI.Interval(l[i]-1, u[i]-1)
        end
    end

    @testset "Broadcasted constraint with indices" begin
        m = ModelType()
        @variable m x[1:2]
        @constraint m cref1[i=2:4] x .== [i, i+1]
        ConstraintRefType = eltype(cref1[2])
        @test cref1 isa JuMP.JuMPArray{AbstractArray{ConstraintRefType}}
        @constraint m cref2[i=1:3, j=1:4] x .≤ [i+j, i-j]
        @test cref2 isa Matrix{AbstractArray{ConstraintRefType}}
        @variable m y[1:2, 1:2]
        @constraint m cref3[i=1:2] x[i,:] .== 1
        @test cref3 isa Vector{AbstractArray{ConstraintRefType}}
    end

    @testset "QuadExpr constraints" begin
        m = ModelType()
        @variable(m, x)
        @variable(m, y)

        cref = @constraint(m, x^2 + x <= 1)
        c = JuMP.constraint_object(cref)
        @test JuMP.isequal_canonical(c.func, x^2 + x)
        @test c.set == MOI.LessThan(1.0)

        cref = @constraint(m, y*x - 1.0 == 0.0)
        c = JuMP.constraint_object(cref)
        @test JuMP.isequal_canonical(c.func, x*y)
        @test c.set == MOI.EqualTo(1.0)

        # TODO: VectorQuadraticFunctions
        # cref = @constraint(m, [x^2 - 1] in MOI.SecondOrderCone(1))
        # c = JuMP.constraint_object(cref)
        # @test JuMP.isequal_canonical(c.func, -1 + x^2)
        # @test c.set == MOI.SecondOrderCone(1)
    end

    @testset "SDP constraint" begin
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

        @SDconstraint(m, cref, [x 1; 1 -y] ⪰ [1 x; x -2])
        @test JuMP.name(cref) == "cref"
        c = JuMP.constraint_object(cref)
        @test JuMP.isequal_canonical(c.func[1], x-1)
        @test JuMP.isequal_canonical(c.func[2], 1-x)
        @test JuMP.isequal_canonical(c.func[3], 2-y)
        @test c.set == MOI.PositiveSemidefiniteConeTriangle(2)
        @test c.shape isa JuMP.SymmetricMatrixShape

        @SDconstraint(m, iref[i=1:2], 0 ⪯ [x+i x+y; x+y -y])
        for i in 1:2
            @test JuMP.name(iref[i]) == "iref[$i]"
            c = JuMP.constraint_object(iref[i])
            @test JuMP.isequal_canonical(c.func[1], x+i)
            @test JuMP.isequal_canonical(c.func[2], x+y)
            @test JuMP.isequal_canonical(c.func[3], -y)
            @test c.set == MOI.PositiveSemidefiniteConeTriangle(2)
            @test c.shape isa JuMP.SymmetricMatrixShape
        end

        # Should throw "ERROR: function JuMP.add_constraint does not accept keyword arguments"
        # This tests that the keyword arguments are passed to add_constraint
        @test_macro_throws ErrorException @SDconstraint(m, [x 1; 1 -y] ⪰ [1 x; x -2], unknown_kw=1)
        # Invalid sense == in SDP constraint
        @test_macro_throws ErrorException @SDconstraint(m, [x 1; 1 -y] == [1 x; x -2])
    end

    @testset "Nonsensical SDPs" begin
        m = ModelType()
        @test_throws ErrorException @variable(m, unequal[1:5,1:6], PSD)
        # Some of these errors happen at compile time, so we can't use @test_throws
        @test_macro_throws ErrorException @variable(m, notone[1:5,2:6], PSD)
        @test_macro_throws ErrorException @variable(m, oneD[1:5], PSD)
        @test_macro_throws ErrorException @variable(m, threeD[1:5,1:5,1:5], PSD)
        @test_macro_throws ErrorException @variable(m, psd[2] <= rand(2,2), PSD)
        @test_macro_throws ErrorException @variable(m, -ones(3,4) <= foo[1:4,1:4] <= ones(4,4), PSD)
        @test_macro_throws ErrorException @variable(m, -ones(3,4) <= foo[1:4,1:4] <= ones(4,4), Symmetric)
        @test_macro_throws ErrorException @variable(m, -ones(4,4) <= foo[1:4,1:4] <= ones(4,5), Symmetric)
        @test_macro_throws ErrorException @variable(m, -rand(5,5) <= nonsymmetric[1:5,1:5] <= rand(5,5), Symmetric)
    end

    @testset "[macros] sum(generator)" begin
        m = ModelType()
        @variable(m, x[1:3,1:3])
        @variable(m, y)
        C = [1 2 3; 4 5 6; 7 8 9]

        @test_expression sum( C[i,j]*x[i,j] for i in 1:2, j = 2:3 )
        @test_expression sum( C[i,j]*x[i,j] for i = 1:3, j in 1:3 if i != j) - y
        @test JuMP.isequal_canonical(@expression(m, sum( C[i,j]*x[i,j] for i = 1:3, j = 1:i)),
                                                    sum( C[i,j]*x[i,j] for i = 1:3 for j = 1:i))
        @test_expression sum( C[i,j]*x[i,j] for i = 1:3 for j = 1:i)
        @test_expression sum( C[i,j]*x[i,j] for i = 1:3 if true for j = 1:i)
        @test_expression sum( C[i,j]*x[i,j] for i = 1:3 if true for j = 1:i if true)
        @test_expression sum( 0*x[i,1] for i=1:3)
        @test_expression sum( 0*x[i,1] + y for i=1:3)
        @test_expression sum( 0*x[i,1] + y for i=1:3 for j in 1:3)
    end
end

@testset "Constraints for JuMP.Model" begin
    constraints_test(Model)
end

@testset "Constraints for JuMPExtension.MyModel" begin
    constraints_test(JuMPExtension.MyModel)
end

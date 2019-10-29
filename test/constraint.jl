function test_constraint_name(constraint, name, F::Type, S::Type)
    @test name == @inferred JuMP.name(constraint)
    model = constraint.model
    @test constraint.index == JuMP.constraint_by_name(model, name).index
    if !(model isa JuMPExtension.MyModel)
        @test constraint.index == JuMP.constraint_by_name(model, name, F, S).index
    end
end

function constraints_test(ModelType::Type{<:JuMP.AbstractModel},
                          VariableRefType::Type{<:JuMP.AbstractVariableRef})
    AffExprType = JuMP.GenericAffExpr{Float64, VariableRefType}
    QuadExprType = JuMP.GenericQuadExpr{Float64, VariableRefType}

    @testset "SingleVariable constraints" begin
        m = ModelType()
        @variable(m, x)

        # x <= 10.0 doesn't translate to a SingleVariable constraint because
        # the LHS is first subtracted to form x - 10.0 <= 0.
        @constraint(m, cref, x in MOI.LessThan(10.0))
        test_constraint_name(cref, "cref", JuMP.VariableRef,
                             MOI.LessThan{Float64})
        c = JuMP.constraint_object(cref)
        @test c.func == x
        @test c.set == MOI.LessThan(10.0)

        @variable(m, y[1:2])
        @constraint(m, cref2[i=1:2], y[i] in MOI.LessThan(float(i)))
        test_constraint_name(cref2[1], "cref2[1]", JuMP.VariableRef,
                             MOI.LessThan{Float64})
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
        @testset "Scalar" begin
            model = ModelType()
            @variable(model, x)

            cref = @constraint(model, 2x <= 10)
            @test "" == @inferred JuMP.name(cref)
            JuMP.set_name(cref, "c")
            test_constraint_name(cref, "c", JuMP.AffExpr, MOI.LessThan{Float64})

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

        @testset "Vectorized" begin
            model = ModelType()
            @variable(model, x)

            @test_throws ErrorException @constraint(model, [x, 2x] == [1-x, 3])
            @test_macro_throws ErrorException begin
                @constraint(model, [x == 1-x, 2x == 3])
            end

            cref = @constraint(model, [x, 2x] .== [1-x, 3])
            c = JuMP.constraint_object.(cref)
            @test JuMP.isequal_canonical(c[1].func, 2.0x)
            @test c[1].set == MOI.EqualTo(1.0)
            @test JuMP.isequal_canonical(c[2].func, 2.0x)
            @test c[2].set == MOI.EqualTo(3.0)
        end

        @testset "Vector" begin
            model = ModelType()

            cref = @constraint(model, [1, 2] in MOI.Zeros(2))
            c = JuMP.constraint_object(cref)
            @test JuMP.isequal_canonical(c.func[1], zero(JuMP.AffExpr) + 1)
            @test JuMP.isequal_canonical(c.func[2], zero(JuMP.AffExpr) + 2)
            @test c.set == MOI.Zeros(2)
            @test c.shape isa JuMP.VectorShape
        end
    end
    @testset "delete / is_valid constraints" begin
        model = ModelType()
        @variable(model, x)
        constraint_ref = @constraint(model, 2x <= 1)
        @test JuMP.is_valid(model, constraint_ref)
        JuMP.delete(model, constraint_ref)
        @test !JuMP.is_valid(model, constraint_ref)
        second_model = ModelType()
        @test_throws Exception JuMP.delete(second_model, constraint_ref)
    end

    @testset "Two-sided constraints" begin
        m = ModelType()
        @variable(m, x)
        @variable(m, y)

        @constraint(m, cref, 1.0 <= x + y + 1.0 <= 2.0)
        test_constraint_name(cref, "cref", JuMP.AffExpr, MOI.Interval{Float64})

        c = JuMP.constraint_object(cref)
        @test JuMP.isequal_canonical(c.func, x + y)
        @test c.set == MOI.Interval(0.0, 1.0)

        cref = @constraint(m, 2x - y + 2.0 ∈ MOI.Interval(-1.0, 1.0))
        @test "" == @inferred JuMP.name(cref)

        c = JuMP.constraint_object(cref)
        @test JuMP.isequal_canonical(c.func, 2x - y)
        @test c.set == MOI.Interval(-3.0, -1.0)
    end

    @testset "Broadcasted constraint (.==)" begin
        m = ModelType()
        @variable(m, x[1:2])

        A = [1.0 2.0; 3.0 4.0]
        b = [4.0, 5.0]

        cref = @constraint(m, A*x .== b)
        @test (2,) == @inferred size(cref)

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
        @test (2,2) == @inferred size(cref)
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
        @test (2,) == @inferred size(cref)

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
        @test cref1 isa JuMP.Containers.DenseAxisArray{Vector{ConstraintRefType}}
        @constraint m cref2[i=1:3, j=1:4] x .≤ [i+j, i-j]
        ConstraintRefType = eltype(cref2[1])
        @test cref2 isa Matrix{Vector{ConstraintRefType}}
        @variable m y[1:2, 1:2]
        @constraint m cref3[i=1:2] y[i,:] .== 1
        ConstraintRefType = eltype(cref3[1])
        @test cref3 isa Vector{Vector{ConstraintRefType}}
    end

    @testset "QuadExpr constraints" begin
        model = ModelType()
        @variable(model, x)
        @variable(model, y)

        cref = @constraint(model, x^2 + x <= 1)
        c = JuMP.constraint_object(cref)
        @test JuMP.isequal_canonical(c.func, x^2 + x)
        @test c.set == MOI.LessThan(1.0)

        cref = @constraint(model, y*x - 1.0 == 0.0)
        c = JuMP.constraint_object(cref)
        @test JuMP.isequal_canonical(c.func, x*y)
        @test c.set == MOI.EqualTo(1.0)

        cref = @constraint(model, [ 2x - 4x*y + 3x^2 - 1,
                                   -3y + 2x*y - 2x^2 + 1] in SecondOrderCone())
        c = JuMP.constraint_object(cref)
        @test JuMP.isequal_canonical(c.func[1], -1 + 3x^2 - 4x*y + 2x)
        @test JuMP.isequal_canonical(c.func[2],  1 - 2x^2 + 2x*y - 3y)
        @test c.set == MOI.SecondOrderCone(2)
    end

    @testset "Syntax error" begin
        model = ModelType()
        @variable(model, x[1:2])
        err = ErrorException(
            "In `@constraint(model, [3, x] in SecondOrderCone())`: unable to " *
            "add the constraint because we don't recognize $([3, x]) as a " *
            "valid JuMP function."
        )
        @test_throws err @constraint(model, [3, x] in SecondOrderCone())
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

        @constraint(m, sym_ref, Symmetric([x 1; 1 -y] - [1 x; x -2]) in PSDCone())
        test_constraint_name(sym_ref, "sym_ref", Vector{AffExpr},
                             MOI.PositiveSemidefiniteConeTriangle)
        c = JuMP.constraint_object(sym_ref)
        @test JuMP.isequal_canonical(c.func[1], x-1)
        @test JuMP.isequal_canonical(c.func[2], 1-x)
        @test JuMP.isequal_canonical(c.func[3], 2-y)
        @test c.set == MOI.PositiveSemidefiniteConeTriangle(2)
        @test c.shape isa JuMP.SymmetricMatrixShape

        @SDconstraint(m, cref, [x 1; 1 -y] ⪰ [1 x; x -2])
        test_constraint_name(cref, "cref", Vector{AffExpr},
                             MOI.PositiveSemidefiniteConeSquare)
        c = JuMP.constraint_object(cref)
        @test JuMP.isequal_canonical(c.func[1], x-1)
        @test JuMP.isequal_canonical(c.func[2], 1-x)
        @test JuMP.isequal_canonical(c.func[3], 1-x)
        @test JuMP.isequal_canonical(c.func[4], 2-y)
        @test c.set == MOI.PositiveSemidefiniteConeSquare(2)
        @test c.shape isa JuMP.SquareMatrixShape

        @SDconstraint(m, iref[i=1:2], 0 ⪯ [x+i x+y; x+y -y])
        for i in 1:2
            test_constraint_name(iref[i], "iref[$i]", Vector{AffExpr},
                                 MOI.PositiveSemidefiniteConeSquare)
            c = JuMP.constraint_object(iref[i])
            @test JuMP.isequal_canonical(c.func[1], x+i)
            @test JuMP.isequal_canonical(c.func[2], x+y)
            @test JuMP.isequal_canonical(c.func[3], x+y)
            @test JuMP.isequal_canonical(c.func[4], -y)
            @test c.set == MOI.PositiveSemidefiniteConeSquare(2)
            @test c.shape isa JuMP.SquareMatrixShape
        end

        # Should throw "ERROR: function JuMP.add_constraint does not accept keyword arguments"
        # This tests that the keyword arguments are passed to add_constraint
        @test_macro_throws ErrorException @SDconstraint(m, [x 1; 1 -y] ⪰ [1 x; x -2], unknown_kw=1)
        # Invalid sense == in SDP constraint
        @test_macro_throws ErrorException @SDconstraint(m, [x 1; 1 -y] == [1 x; x -2])
    end

    @testset "Constraint name" begin
        model = ModelType()
        @variable(model, x)
        @constraint(model, con, x^2 == 1)
        test_constraint_name(con, "con", QuadExprType, MOI.EqualTo{Float64})
        JuMP.set_name(con, "kon")
        @test JuMP.constraint_by_name(model, "con") isa Nothing
        test_constraint_name(con, "kon", QuadExprType, MOI.EqualTo{Float64})
        y = @constraint(model, kon, [x^2, x] in SecondOrderCone())
        err(name) = ErrorException("Multiple constraints have the name $name.")
        @test_throws err("kon") JuMP.constraint_by_name(model, "kon")
        JuMP.set_name(kon, "con")
        test_constraint_name(con, "kon", QuadExprType, MOI.EqualTo{Float64})
        test_constraint_name(kon, "con", Vector{QuadExprType},
                             MOI.SecondOrderCone)
        JuMP.set_name(con, "con")
        @test_throws err("con") JuMP.constraint_by_name(model, "con")
        @test JuMP.constraint_by_name(model, "kon") isa Nothing
        JuMP.set_name(kon, "kon")
        test_constraint_name(con, "con", QuadExprType, MOI.EqualTo{Float64})
        test_constraint_name(kon, "kon", Vector{QuadExprType},
                             MOI.SecondOrderCone)
    end

    @testset "Useful PSD error message" begin
        model = ModelType()
        @variable(model, X[1:2, 1:2])
        err = ErrorException(
            "In `@constraint(model, X in MOI.PositiveSemidefiniteConeSquare(2))`:" *
            " instead of `MathOptInterface.PositiveSemidefiniteConeSquare(2)`," *
            " use `JuMP.PSDCone()`.")
        @test_throws err @constraint(model, X in MOI.PositiveSemidefiniteConeSquare(2))
        err = ErrorException(
            "In `@constraint(model, X in MOI.PositiveSemidefiniteConeTriangle(2))`:" *
            " instead of `MathOptInterface.PositiveSemidefiniteConeTriangle(2)`," *
            " use `JuMP.PSDCone()`.")
        @test_throws err @constraint(model, X in MOI.PositiveSemidefiniteConeTriangle(2))
    end

    @testset "Useful Matrix error message" begin
        model = ModelType()
        @variable(model, X[1:2, 1:2])
        err = ErrorException(
            "In `@constraint(model, X in MOI.SecondOrderCone(4))`: unexpected " *
            "matrix in vector constraint. Do you need to flatten the matrix " *
            "into a vector using `vec()`?")
        # Note: this should apply to any MOI.AbstractVectorSet. We just pick
        # SecondOrderCone for convenience.
        @test_throws err @constraint(model, X in MOI.SecondOrderCone(4))
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
        model = ModelType()
        @variable(model, x[1:3,1:3])
        @variable(model, y)
        C = [1 2 3; 4 5 6; 7 8 9]

        @test_expression sum( C[i,j]*x[i,j] for i in 1:2, j = 2:3 )
        @test_expression sum( C[i,j]*x[i,j] for i = 1:3, j in 1:3 if i != j) - y
        @test JuMP.isequal_canonical(@expression(model, sum( C[i,j]*x[i,j] for i = 1:3, j = 1:i)),
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
    constraints_test(Model, JuMP.VariableRef)
    @testset "all_constraints (scalar)" begin
        model = Model()
        @variable(model, x >= 0)
        @test 1 == @inferred num_constraints(model, VariableRef,
                                             MOI.GreaterThan{Float64})
        ref = @inferred all_constraints(
            model, VariableRef, MOI.GreaterThan{Float64})
        @test ref == [LowerBoundRef(x)]
        @test 0 == @inferred num_constraints(model, AffExpr,
                                             MOI.GreaterThan{Float64})
        aff_constraints = all_constraints(model, AffExpr,
                                          MOI.GreaterThan{Float64})
        @test isempty(aff_constraints)
        err = ErrorException("`MathOptInterface.GreaterThan` is not a " *
                             "concrete type. Did you miss a type parameter?")
        @test_throws err num_constraints(model, AffExpr,
                                         MOI.GreaterThan)
        @test_throws err all_constraints(model, AffExpr,
                                         MOI.GreaterThan)
        err = ErrorException("`GenericAffExpr` is not a concrete type. " *
                             "Did you miss a type parameter?")
        @test_throws err num_constraints(model, GenericAffExpr,
                                         MOI.ZeroOne)
        @test_throws err all_constraints(model, GenericAffExpr,
                                         MOI.ZeroOne)
    end
    @testset "all_constraints (vector)" begin
        model = Model()
        @variable(model, x[1:2, 1:2], Symmetric)
        csdp = @constraint(model, x in PSDCone())
        csoc = @constraint(model, [x[1], 1] in SecondOrderCone())
        csos = @constraint(model, [x[2]^2, 1] in MOI.SOS1([1.0, 2.0]))
        @test 1 == @inferred num_constraints(
            model, Vector{VariableRef}, MOI.PositiveSemidefiniteConeTriangle)
        ref = all_constraints(model, Vector{VariableRef},
                              MOI.PositiveSemidefiniteConeTriangle)
        @test ref == [csdp]
        @test 1 == @inferred num_constraints(
            model, Vector{AffExpr}, MOI.SecondOrderCone)
        ref = all_constraints(model, Vector{AffExpr}, MOI.SecondOrderCone)
        @test ref == [csoc]
        @test 1 == @inferred num_constraints(
             model, Vector{QuadExpr}, MOI.SOS1{Float64})
        ref = all_constraints(model, Vector{QuadExpr}, MOI.SOS1{Float64})
        @test ref == [csos]
        @test 0 == @inferred num_constraints(
            model, Vector{AffExpr}, MOI.PositiveSemidefiniteConeTriangle)
        aff_constraints = all_constraints(
            model, Vector{AffExpr}, MOI.PositiveSemidefiniteConeTriangle)
        @test isempty(aff_constraints)
        err = ErrorException("`GenericAffExpr{Float64,VarType} where VarType`" *
                             " is not a concrete type. Did you miss a type " *
                             "parameter?")
        @test_throws err num_constraints(
            model, Vector{GenericAffExpr{Float64}},
            MOI.PositiveSemidefiniteConeTriangle)
        @test_throws err all_constraints(
            model, Vector{GenericAffExpr{Float64}},
            MOI.SecondOrderCone)
        err = ErrorException("`MathOptInterface.SOS1` is not a " *
                             "concrete type. Did you miss a type parameter?")
        @test_throws err all_constraints(
            model, Vector{GenericQuadExpr{Float64,VariableRef}},
            MOI.SOS1)
    end
    @testset "list_of_constraint_types" begin
        model = Model()
        @variable(model, x >= 0, Bin)
        @constraint(model, 2x <= 1)
        @constraint(model, [x, x] in SecondOrderCone())
        @constraint(model, [2x  1; 1 x] in PSDCone())
        @constraint(model, [x^2, x] in RotatedSecondOrderCone())
        constraint_types = @inferred list_of_constraint_types(model)
        @test Set(constraint_types) == Set([(VariableRef, MOI.ZeroOne),
            (VariableRef, MOI.GreaterThan{Float64}),
            (AffExpr, MOI.LessThan{Float64}),
            (Vector{VariableRef}, MOI.SecondOrderCone),
            (Vector{AffExpr}, MOI.PositiveSemidefiniteConeSquare),
            (Vector{QuadExpr}, MOI.RotatedSecondOrderCone)])
    end
end

@testset "Constraints for JuMPExtension.MyModel" begin
    constraints_test(JuMPExtension.MyModel, JuMPExtension.MyVariableRef)
end

@testset "Modifications" begin
    @testset "Change coefficient" begin
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
            JuMP.constraint_object(quad_con).func, x^2 + 2x)
    end

    @testset "Change rhs" begin
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

    @testset "Add to function constant" begin
        model = JuMP.Model()
        x = @variable(model)
        @testset "Scalar" begin
            con_ref = @constraint(model, 2 <= 2 * x <= 3)
            con = constraint_object(con_ref)
            @test JuMP.isequal_canonical(JuMP.jump_function(con), 2x)
            @test JuMP.moi_set(con) == MOI.Interval(2.0, 3.0)
            JuMP.add_to_function_constant(con_ref, 1.0)
            con = constraint_object(con_ref)
            @test JuMP.isequal_canonical(JuMP.jump_function(con), 2x)
            @test JuMP.moi_set(con) == MOI.Interval(1.0, 2.0)
        end
        @testset "Vector" begin
            con_ref = @constraint(model, [x + 1, x - 1] in MOI.Nonnegatives(2))
            con = constraint_object(con_ref)
            @test JuMP.isequal_canonical(JuMP.jump_function(con), [x + 1, x - 1])
            @test JuMP.moi_set(con) == MOI.Nonnegatives(2)
            JuMP.add_to_function_constant(con_ref, [2, 3])
            con = constraint_object(con_ref)
            @test JuMP.isequal_canonical(JuMP.jump_function(con), [x + 3, x + 2])
            @test JuMP.moi_set(con) == MOI.Nonnegatives(2)
        end
    end

end

function test_shadow_price(model_string, constraint_dual, constraint_shadow)
    model = JuMP.Model()
    MOIU.loadfromstring!(JuMP.backend(model), model_string)
    set_optimizer(model, () -> MOIU.MockOptimizer(
                                MOIU.Model{Float64}(),
                                eval_objective_value=false,
                                eval_variable_constraint_dual=false))
    JuMP.optimize!(model)
    mock_optimizer = JuMP.backend(model).optimizer.model
    MOI.set(mock_optimizer, MOI.TerminationStatus(), MOI.OPTIMAL)
    MOI.set(mock_optimizer, MOI.DualStatus(), MOI.FEASIBLE_POINT)
    JuMP.optimize!(model)

    @testset "shadow price of $constraint_name" for constraint_name in keys(constraint_dual)
        ci = MOI.get(JuMP.backend(model), MOI.ConstraintIndex,
                     constraint_name)
        constraint_ref = JuMP.ConstraintRef(model, ci, JuMP.ScalarShape())
        MOI.set(mock_optimizer, MOI.ConstraintDual(),
                JuMP.optimizer_index(constraint_ref),
                constraint_dual[constraint_name])
        @test JuMP.dual(constraint_ref) ==
              constraint_dual[constraint_name]
        @test JuMP.shadow_price(constraint_ref) ==
              constraint_shadow[constraint_name]
    end
end

@testset "shadow_price" begin
    test_shadow_price("""
    variables: x, y
    minobjective: -1.0*x
    xub: x <= 2.0
    ylb: y >= 0.0
    c: x + y <= 1.0
    """,
    # Optimal duals
    Dict("xub" => 0.0, "ylb" => 1.0, "c" => -1.0),
    # Expected shadow prices
    Dict("xub" => 0.0, "ylb" => -1.0, "c" => -1.0))

    test_shadow_price("""
    variables: x, y
    maxobjective: 1.0*x
    xub: x <= 2.0
    ylb: y >= 0.0
    c: x + y <= 1.0
    """,
    # Optimal duals
    Dict("xub" => 0.0, "ylb" => 1.0, "c" => -1.0),
    # Expected shadow prices
    Dict("xub" => 0.0, "ylb" => 1.0, "c" => 1.0))

    test_shadow_price("""
    variables: x, y
    maxobjective: 1.0*x
    xub: x <= 2.0
    ylb: y >= 0.0
    """,
    # Optimal duals
    Dict("xub" => -1.0, "ylb" => 0.0),
    # Expected shadow prices
    Dict("xub" => 1.0, "ylb" => 0.0))

    test_shadow_price("""
    variables: x
    maxobjective: 1.0*x
    xeq: x == 2.0
    """,
    # Optimal duals
    Dict("xeq" => -1.0),
    # Expected shadow prices
    Dict("xeq" => 1.0))

    test_shadow_price("""
    variables: x
    minobjective: 1.0*x
    xeq: x == 2.0
    """,
    # Optimal duals
    Dict("xeq" => 1.0),
    # Expected shadow prices
    Dict("xeq" => -1.0))
end

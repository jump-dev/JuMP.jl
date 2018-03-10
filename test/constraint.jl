

@testset "Constraints" begin
    @testset "SingleVariable constraints" begin
        m = Model()
        @variable(m, x)

        # x <= 10.0 doesn't translate to a SingleVariable constraint because
        # the LHS is first subtracted to form x - 10.0 <= 0.
        @constraint(m, cref, x in MOI.LessThan(10.0))
        @test JuMP.name(cref) == "cref"
        c = JuMP.constraintobject(cref, Variable, MOI.LessThan)
        @test c.func == x
        @test c.set == MOI.LessThan(10.0)
        @test_throws TypeError JuMP.constraintobject(cref, QuadExpr, MOI.LessThan)
        @test_throws TypeError JuMP.constraintobject(cref, AffExpr, MOI.EqualTo)

        @variable(m, y[1:2])
        @constraint(m, cref2[i=1:2], y[i] in MOI.LessThan(float(i)))
        @test JuMP.name(cref2[1]) == "cref2[1]"
        c = JuMP.constraintobject(cref2[1], Variable, MOI.LessThan)
        @test c.func == y[1]
        @test c.set == MOI.LessThan(1.0)
    end

    @testset "VectorOfVariables constraints" begin
        m = Model()
        @variable(m, x[1:2])

        cref = @constraint(m, x in MOI.Zeros(2))
        c = JuMP.constraintobject(cref, Vector{Variable}, MOI.Zeros)
        @test c.func == x
        @test c.set == MOI.Zeros(2)
        @test_throws TypeError JuMP.constraintobject(cref, Vector{AffExpr}, MOI.Nonnegatives)
        @test_throws TypeError JuMP.constraintobject(cref, AffExpr, MOI.EqualTo)

        cref = @constraint(m, [x[2],x[1]] in MOI.Zeros(2))
        c = JuMP.constraintobject(cref, Vector{Variable}, MOI.Zeros)
        @test c.func == [x[2],x[1]]
        @test c.set == MOI.Zeros(2)
    end

    @testset "AffExpr constraints" begin
        m = Model()
        @variable(m, x)

        cref = @constraint(m, 2x <= 10)
        @test JuMP.name(cref) == ""
        JuMP.setname(cref, "c")
        @test JuMP.name(cref) == "c"

        c = JuMP.constraintobject(cref, AffExpr, MOI.LessThan)
        @test JuMP.isequal_canonical(c.func, 2x)
        @test c.set == MOI.LessThan(10.0)
        @test_throws TypeError JuMP.constraintobject(cref, QuadExpr, MOI.LessThan)
        @test_throws TypeError JuMP.constraintobject(cref, AffExpr, MOI.EqualTo)

        cref = @constraint(m, 3x + 1 ≥ 10)
        c = JuMP.constraintobject(cref, AffExpr, MOI.GreaterThan)
        @test JuMP.isequal_canonical(c.func, 3x)
        @test c.set == MOI.GreaterThan(9.0)

        cref = @constraint(m, 1 == -x)
        c = JuMP.constraintobject(cref, AffExpr, MOI.EqualTo)
        @test JuMP.isequal_canonical(c.func, 1.0x)
        @test c.set == MOI.EqualTo(-1.0)

        @test_throws ErrorException @constraint(m, [x, 2x] == [1-x, 3])
        cref = @constraint(m, [x, 2x] .== [1-x, 3])
        c = JuMP.constraintobject.(cref, AffExpr, MOI.EqualTo)
        @test JuMP.isequal_canonical(c[1].func, 2.0x)
        @test c[1].set == MOI.EqualTo(1.0)
        @test JuMP.isequal_canonical(c[2].func, 2.0x)
        @test c[2].set == MOI.EqualTo(3.0)

        @test MOI.isvalid(m, cref[1])
        MOI.delete!(m, cref[1])
        @test !MOI.isvalid(m, cref[1])
    end

    @testset "Two-sided constraints" begin
        m = Model()
        @variable(m, x)
        @variable(m, y)

        @constraint(m, cref, 1.0 <= x + y + 1.0 <= 2.0)
        @test JuMP.name(cref) == "cref"

        c = JuMP.constraintobject(cref, AffExpr, MOI.Interval)
        @test JuMP.isequal_canonical(c.func, x + y)
        @test c.set == MOI.Interval(0.0, 1.0)
    end

    @testset "Broadcasted constraint (.==)" begin
        m = Model()
        @variable(m, x[1:2])

        A = [1.0 2.0; 3.0 4.0]
        b = [4.0, 5.0]

        cref = @constraint(m, A*x .== b)
        @test size(cref) == (2,)

        c1 = JuMP.constraintobject(cref[1], AffExpr, MOI.EqualTo)
        @test JuMP.isequal_canonical(c1.func, x[1] + 2x[2])
        @test c1.set == MOI.EqualTo(4.0)
        c2 = JuMP.constraintobject(cref[2], AffExpr, MOI.EqualTo)
        @test JuMP.isequal_canonical(c2.func, 3x[1] + 4x[2])
        @test c2.set == MOI.EqualTo(5.0)
    end

    @testset "Broadcasted constraint (.<=)" begin
        m = Model()
        @variable(m, x[1:2,1:2])

        UB = [1.0 2.0; 3.0 4.0]

        cref = @constraint(m, x + 1 .<= UB)
        @test size(cref) == (2,2)
        for i in 1:2
            for j in 1:2
                c = JuMP.constraintobject(cref[i,j], AffExpr, MOI.LessThan)
                @test JuMP.isequal_canonical(c.func, AffExpr(x[i,j]))
                @test c.set == MOI.LessThan(UB[i,j] - 1)
            end
        end
    end

    @testset "Broadcasted two-sided constraint" begin
        m = Model()
        @variable(m, x[1:2])
        @variable(m, y[1:2])
        l = [1.0, 2.0]
        u = [3.0, 4.0]

        cref = @constraint(m, l .<= x + y + 1 .<= u)
        @test size(cref) == (2,)

        for i in 1:2
            c = JuMP.constraintobject(cref[i], AffExpr, MOI.Interval)
            @test JuMP.isequal_canonical(c.func, x[i] + y[i])
            @test c.set == MOI.Interval(l[i]-1, u[i]-1)
        end
    end

    @testset "QuadExpr constraints" begin
        m = Model()
        @variable(m, x)
        @variable(m, y)

        cref = @constraint(m, x^2 + x <= 1)
        c = JuMP.constraintobject(cref, QuadExpr, MOI.LessThan)
        @test JuMP.isequal_canonical(c.func, x^2 + x)
        @test c.set == MOI.LessThan(1.0)

        cref = @constraint(m, y*x - 1.0 == 0.0)
        c = JuMP.constraintobject(cref, QuadExpr, MOI.EqualTo)
        @test JuMP.isequal_canonical(c.func, x*y)
        @test c.set == MOI.EqualTo(1.0)
        @test_throws TypeError JuMP.constraintobject(cref, QuadExpr, MOI.LessThan)
        @test_throws TypeError JuMP.constraintobject(cref, AffExpr, MOI.EqualTo)

        # cref = @constraint(m, [x^2 - 1] in MOI.SecondOrderCone(1))
        # c = JuMP.constraintobject(cref, QuadExpr, MOI.SecondOrderCone)
        # @test JuMP.isequal_canonical(c.func, -1 + x^2)
        # @test c.set == MOI.SecondOrderCone(1)
    end

    @testset "SDP constraint" begin
        m = Model()
        @variable(m, x)
        @variable(m, y)
        @variable(m, z)
        @variable(m, w)

        cref = @constraint(m, [x y; z w] in PSDCone())
        c = JuMP.constraintobject(cref, Vector{Variable}, MOI.PositiveSemidefiniteConeSquare)
        @test c.func == [x, z, y, w]
        @test c.set == MOI.PositiveSemidefiniteConeSquare(2)

        cref = @SDconstraint(m, [x 1; 1 -y] ⪰ [1 x; x -2])
        c = JuMP.constraintobject(cref, Vector{AffExpr}, MOI.PositiveSemidefiniteConeTriangle)
        @test JuMP.isequal_canonical(c.func[1], x-1)
        @test JuMP.isequal_canonical(c.func[2], 1-x)
        @test JuMP.isequal_canonical(c.func[3], 2-y)
        @test c.set == MOI.PositiveSemidefiniteConeTriangle(2)
    end

    @testset "Nonsensical SDPs" begin
        m = Model()
        @test_throws ErrorException @variable(m, unequal[1:5,1:6], PSD)
        # Some of these errors happen at compile time, so we can't use @test_throws
        @test macroexpand(:(@variable(m, notone[1:5,2:6], PSD))).head == :error
        @test macroexpand(:(@variable(m, oneD[1:5], PSD))).head == :error
        @test macroexpand(:(@variable(m, threeD[1:5,1:5,1:5], PSD))).head == :error
        @test macroexpand(:(@variable(m, psd[2] <= rand(2,2), PSD))).head == :error
        @test macroexpand(:(@variable(m, -ones(3,4) <= foo[1:4,1:4] <= ones(4,4), PSD))).head == :error
        @test macroexpand(:(@variable(m, -ones(3,4) <= foo[1:4,1:4] <= ones(4,4), Symmetric))).head == :error
        @test macroexpand(:(@variable(m, -ones(4,4) <= foo[1:4,1:4] <= ones(4,5), Symmetric))).head == :error
        @test macroexpand(:(@variable(m, -rand(5,5) <= nonsymmetric[1:5,1:5] <= rand(5,5), Symmetric))).head == :error
    end

end

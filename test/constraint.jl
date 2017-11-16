

@testset "Constraints" begin
    @testset "AffExpr constraints" begin
        m = Model()
        @variable(m, x)

        cref = @constraint(m, 2x <= 10)
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

    @testset "Trivial symmetry constraints are removed (#766, #972)" begin
        q = 2
        m = 3
        angles1 = linspace(3*pi/4, pi, m)
        angles2 = linspace(0, -pi/2, m)
        V = [3.*cos.(angles1)' 1.5.*cos.(angles2)';
             3.*sin.(angles1)' 1.5.*sin.(angles2)']
        V[abs.(V) .< 1e-10] = 0.0
        p = 2*m
        n = 100

        mod = Model()
        @variable(mod, x[j=1:p] >= 1, Int)
        @variable(mod, u[i=1:q] >= 1)
        @objective(mod, Min, sum(u))
        @constraint(mod, sum(x) <= n)
        for i=1:q
            cref = @SDconstraint(mod, [V*diagm(x./n)*V' eye(q)[:,i] ; eye(q)[i:i,:] u[i]] >= 0)
            @test cref.instanceref.symref === nothing
            @test isempty(cref.instanceref.symidx)
        end
    end

end

@testset "Second-order Cone Programming" begin
    @testset "SOC1" begin
        m = Model()
        @variable(m, x)
        @variable(m, y)
        @variable(m, t >= 0)
        @objective(m, Min, t)
        @constraint(m, x + y >= 1)
        @constraint(m, [t,x,y] in MOI.SecondOrderCone(3))

        JuMP.attach(m, CSDPInstance(printlevel=0))
        JuMP.solve(m)

        @test JuMP.isattached(m)
        @test JuMP.hasvariableresult(m)
        @test JuMP.terminationstatus(m) == MOI.Success
        @test JuMP.primalstatus(m) == MOI.FeasiblePoint

        @test JuMP.objectivevalue(m) ≈ sqrt(1/2) atol=1e-6

        @test JuMP.resultvalue.([x,y,t]) ≈ [0.5,0.5,sqrt(1/2)] atol=1e-3
        @test JuMP.resultvalue(x + y) ≈ 1.0 atol=1e-3

    end

    @testset "RotatedSOC1" begin
        m = Model()

        @variable(m, x[1:5] >= 0)
        @variable(m, 0 <= u <= 5)
        @variable(m, v)
        @variable(m, t1 == 1)
        @variable(m, t2 == 1)

        @objective(m, Max, v)

        @constraint(m, [t1/sqrt(2),t2/sqrt(2),x...] in MOI.RotatedSecondOrderCone(7))
        @constraint(m, [x[1]/sqrt(2), u/sqrt(2), v] in MOI.RotatedSecondOrderCone(3))

        JuMP.attach(m, CSDPInstance(printlevel=0))
        JuMP.solve(m)

        @test JuMP.isattached(m)
        @test JuMP.hasvariableresult(m)
        @test JuMP.terminationstatus(m) == MOI.Success
        @test JuMP.primalstatus(m) == MOI.FeasiblePoint

        @test JuMP.resultvalue.(x) ≈ [1,0,0,0,0] atol=1e-2
        @test JuMP.resultvalue(u) ≈ 5 atol=1e-4
        @test JuMP.resultvalue(v) ≈ sqrt(5) atol=1e-6
    end
end

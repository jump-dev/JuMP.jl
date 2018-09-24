using JuMP, Compat.Test, Compat

!isdefined(@__MODULE__, :conic_solvers_with_duals) && include("solvers.jl")

const TOL = 1e-4

@testset "SOC duals" begin

    @testset "SOC dual with $solver" for solver in conic_solvers_with_duals
        m = Model(solver=solver)

        @variable(m, x[1:5])
        @objective(m, Max, x[1] + x[2] + x[3] + x[4] + 2*x[5])
        @constraint(m, constrcon1, norm(x[2:5]) <= x[1])
        @constraint(m, constrlin1, x[1] <= 5)
        @constraint(m, constrlin2, x[1]+2*x[2] - x[3] <= 10)
        @constraint(m, constrcon2, norm([2 3;1 1]*x[2:3]-[3;4]) <= x[5] - 2)

        solve(m)
        @test length(getdual(constrcon1)) == 5
        @test length(getdual(constrcon2)) == 3
        @test length(m.conicconstrDuals) == 10
        @test length(m.linconstrDuals) == 2

        @test dot(getdual(constrcon2),[getvalue(x[5]) - 2;[2 3;1 1]*[getvalue(x[2]);getvalue(x[3])]-[3;4]]) ≤ TOL

    end

    # Test for consistency between LP duals and SOC duals

    @testset "LP dual vs SOC dual (Max) with $lp_solver" for lp_solver in lp_solvers

        m1 = Model(solver=lp_solver)
        @variable(m1, x1[1:2] >= 0)
        @objective(m1, Max, x1[1] + 2x1[2])
        @constraint(m1, c11, 3x1[1] + x1[2] <= 4)
        @constraint(m1, c12, x1[1] + 2x1[2] >= 1)
        @constraint(m1, c13, -x1[1] + x1[2] == 0.5)

        solve(m1)

        @test isapprox(getdual(c11), 0.75, atol=TOL)
        @test isapprox(getdual(c12), 0.0, atol=TOL)
        @test isapprox(getdual(c13), 1.25, atol=TOL)
    end
    @testset "LP dual vs SOC dual (Max) with $conic_solver" for conic_solver in conic_solvers_with_duals

        m2 = Model(solver=conic_solver)
        @variable(m2, x2[1:2] >= 0)
        @objective(m2, Max, x2[1] + 2x2[2])
        @constraint(m2, c21, 3x2[1] + x2[2] <= 4)
        @constraint(m2, c22, x2[1] + 2x2[2] >= 1)
        @constraint(m2, c23, -x2[1] + x2[2] == 0.5)
        @constraint(m2, norm(x2[1]) <= x2[2])

        solve(m2)

        @test isapprox(getdual(c21), 0.75, atol=TOL)
        @test isapprox(getdual(c22), 0.0, atol=TOL)
        @test isapprox(getdual(c23), 1.25, atol=TOL)
    end

    @testset "LP vs SOC dual (Min) with $lp_solver" for lp_solver in lp_solvers

        m1 = Model(solver=lp_solver)
        @variable(m1, x1[1:2] >= 0)
        @objective(m1, Min, -x1[1] - 2x1[2])
        @constraint(m1, c11, 3x1[1] + x1[2] <= 4)
        @constraint(m1, c12, x1[1] + 2x1[2] >= 1)
        @constraint(m1, c13, -x1[1] + x1[2] == 0.5)

        solve(m1)

        @test isapprox(getdual(c11), -0.75, atol=TOL)
        @test isapprox(getdual(c12), -0.0, atol=TOL)
        @test isapprox(getdual(c13), -1.25, atol=TOL)

    end
    @testset "LP vs SOC dual (Min) with $conic_solver" for  conic_solver in conic_solvers_with_duals

        m2 = Model(solver=conic_solver)
        @variable(m2, x2[1:2] >= 0)
        @objective(m2, Min, -x2[1] - 2x2[2])
        @constraint(m2, c21, 3x2[1] + x2[2] <= 4)
        @constraint(m2, c22, x2[1] + 2x2[2] >= 1)
        @constraint(m2, c23, -x2[1] + x2[2] == 0.5)
        @constraint(m2, norm(x2[1]) <= x2[2])

        solve(m2)
        @test isapprox(getdual(c21), -0.75, atol=TOL)
        @test isapprox(getdual(c22), -0.0, atol=TOL)
        @test isapprox(getdual(c23), -1.25, atol=TOL)
    end

    @testset "LP vs SOC reduced costs with $lp_solver" for lp_solver in lp_solvers
        m1 = Model(solver=lp_solver)

        @variable(m1, x1 >= 0)
        @variable(m1, y1 >= 0)
        @constraint(m1, x1 + y1 == 1)
        @objective(m1, Max, y1)

        solve(m1)

        @test dot([getdual(y1), getdual(x1)],[getvalue(y1); getvalue(x1)]) ≤ TOL

        @test isapprox(getvalue(x1), 0.0, atol=TOL)
        @test isapprox(getvalue(y1), 1.0, atol=TOL)

        m1 = Model(solver=lp_solver)

        @variable(m1, x1 >= 0)
        @constraint(m1, x1 <= 1)
        @objective(m1, Max, x1)

        solve(m1)

        @test isapprox(getdual(x1), 0.0, atol=TOL)
        @test isapprox(getvalue(x1), 1.0, atol=TOL)

        RHS = 2
        X_LB = 0.5
        Y_UB = 5

        m1 = Model(solver=lp_solver)

        @variable(m1, x1 >= X_LB)
        @variable(m1, y1 <= Y_UB)
        @constraint(m1, c1,  x1 + y1 == RHS)
        @objective(m1, Max, 0.8*y1 + 0.1*x1)

        solve(m1)

        @test isapprox(getvalue(x1), 0.5, atol=TOL)
        @test isapprox(getvalue(y1), 1.5, atol=TOL)

        @test isapprox(getdual(x1), -0.7, atol=TOL)
        @test isapprox(getdual(y1), 0.0, atol=TOL)

        lp_dual_obj = getdual(c1)*RHS + getdual(x1) * X_LB + getdual(y1) * Y_UB
        @test isapprox(getobjectivevalue(m1), lp_dual_obj, atol=TOL)
    end
    @testset "LP vs SOC reduced costs with $conic_solver" for conic_solver in conic_solvers_with_duals
        m2 = Model(solver=conic_solver)

        @variable(m2, x2 >= 0)
        @variable(m2, y2 >= 0)
        @constraint(m2, x2 + y2 == 1)
        @objective(m2, Max, y2)
        @constraint(m2, norm(x2) <= y2)

        solve(m2)

        @test dot([getdual(y2), getdual(x2)],[getvalue(y2); getvalue(x2)]) ≤ TOL

        @test isapprox(getvalue(x2), 0.0, atol=TOL)
        @test isapprox(getvalue(y2), 1.0, atol=TOL)

        m2 = Model(solver=conic_solver)

        @variable(m2, x2 >= 0)
        @constraint(m2, x2 <= 1)
        @objective(m2, Max, x2)
        @constraint(m2, norm(x2) <= x2)

        solve(m2)

        @test isapprox(getdual(x2), 0.0, atol=TOL)
        @test isapprox(getvalue(x2), 1.0, atol=TOL)

        RHS = 2
        X_LB = 0.5
        Y_UB = 5

        m2 = Model(solver=conic_solver)

        @variable(m2, x2 >= X_LB)
        @variable(m2, y2 <= Y_UB)
        @constraint(m2, c2, x2 + y2 == RHS)
        @objective(m2, Max, 0.8*y2 + 0.1*x2)
        @constraint(m2, norm(x2) <= y2)

        solve(m2)

        @test isapprox(getvalue(x2), 0.5, atol=TOL)
        @test isapprox(getvalue(y2), 1.5, atol=TOL)

        @test isapprox(getdual(x2), -0.7, atol=TOL)
        @test isapprox(getdual(y2), 0.0, atol=TOL)

        conic_dual_obj = getdual(c2)*RHS + getdual(x2) * X_LB + getdual(y2) * Y_UB
        @test isapprox(getobjectivevalue(m2), conic_dual_obj, atol=TOL)

    end

    # || b-Ax || <= a - c^T x, x\in K
    # (a - c^Tx, b-Ax) \in SOC, x\in K
    # ((a,b) - (c^T,A)x) \in SOC, x\in K
    # if the problem is infeasible, we get a y satisfying -(a,b)^Ty > 0, (c,A^T)y \in K^*, y \in SOC^*

    @testset "SOC infeasibility ray with $conic_solver" for conic_solver in conic_solvers_with_duals
        m2 = Model(solver=conic_solver)

        @variable(m2, x2 >= 0)
        @constraint(m2, x2 <= 1)
        @objective(m2, Max, x2)
        @constraint(m2, c2, norm(x2) <= x2 - 1)

        status = solve(m2, suppress_warnings=true)

        inf_ray = getdual(c2)
        @test status == :Infeasible
        @test inf_ray[1] ≥ abs(inf_ray[2]) - TOL
        @test -(-inf_ray[1]) ≥ TOL

    end
end

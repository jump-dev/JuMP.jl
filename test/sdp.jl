@testset "Semidefinite Programming" begin

ispsd(x::Matrix) = minimum(eigvals(x)) ≥ -1e-3

    @testset "SDP1" begin
        # Problem SDP1 - sdo1 from MOSEK docs
        # From Mosek.jl/test/mathprogtestextra.jl, under license:
        #   Copyright (c) 2013 Ulf Worsoe, Mosek ApS
        #   Permission is hereby granted, free of charge, to any person obtaining a copy of this
        #   software and associated documentation files (the "Software"), to deal in the Software
        #   without restriction, including without limitation the rights to use, copy, modify, merge,
        #   publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons
        #   to whom the Software is furnished to do so, subject to the following conditions:
        #   The above copyright notice and this permission notice shall be included in all copies or
        #   substantial portions of the Software.
        #
        #     | 2 1 0 |
        # min | 1 2 1 | . X + x1
        #     | 0 1 2 |
        #
        #
        # s.t. | 1 0 0 |
        #      | 0 1 0 | . X + x1 = 1
        #      | 0 0 1 |
        #
        #      | 1 1 1 |
        #      | 1 1 1 | . X + x2 + x3 = 1/2
        #      | 1 1 1 |
        #
        #      (x1,x2,x3) in C^3_q
        #      X in C_sdp

        m = Model(solver=CSDPSolver(printlevel=0))

        @variable(m, x[1:3])
        @constraint(m, x in MOI.SecondOrderCone(3))
        @variable(m, X[1:3, 1:3], PSD)

        C = [2 1 0
             1 2 1
             0 1 2]
        A1 = [1 0 0
              0 1 0
              0 0 1]
        A2 = [1 1 1
              1 1 1
              1 1 1]

        @objective(m, Min, vecdot(C, X) + x[1])

        @constraint(m, vecdot(A1, X) + x[1] == 1)
        @constraint(m, vecdot(A2, X) + x[2] + x[3] == 1/2)

        JuMP.solve(m)

        @test JuMP.isattached(m)
        @test JuMP.hasvariableresult(m)
        @test JuMP.terminationstatus(m) == MOI.Success
        @test JuMP.primalstatus(m) == MOI.FeasiblePoint

        @test JuMP.objectivevalue(m) ≈ 0.705710509 atol=1e-6

        xv = JuMP.resultvalue.(x)
        Xv = JuMP.resultvalue.(X)
        @test vecdot(C, Xv) + xv[1] ≈ 0.705710509 atol=1e-6

        @test eigmin(Xv) > -1e-6
    end

#    @testset "Simple SDP" begin
#        #contains(string(typeof(solver)),"SCSSolver") && continue
#        solver = CSDPSolver(printlevel=0)
#        m = Model(solver=solver)
#        @variable(m, X[1:3,1:3], PSD)
#        @SDconstraint(m, X <= 1/2*eye(3,3))
#        @variable(m, Y[1:5,1:5], Symmetric)
#        @SDconstraint(m, -ones(5,5) <= Y)
#        @SDconstraint(m, Y <= 2*ones(5,5))
#        @variable(m, Z[1:4,1:4], Symmetric)
#        @SDconstraint(m, ones(4,4) >= Z)
#
#        @constraint(m, trace(X) == 1)
#        @constraint(m, trace(Y) == 3)
#        @constraint(m, trace(Z) == -1)
#
#        #@test JuMP.numsdconstr(m) == 4
#
#        @objective(m, Max, X[1,2] + Y[1,2] + Z[1,2])
#        JuMP.solve(m)
#        @test JuMP.terminationstatus(m) == MOI.Success
#        @test JuMP.primalstatus(m) == MOI.FeasiblePoint
#
#        XX, YY, ZZ, = JuMP.resultvalue.(X), JuMP.resultvalue.(Y), JuMP.resultvalue.(Z)
#        Xtrue = [0.25 0.25 0
#                 0.25 0.25 0
#                 0    0    0.5]
#        Ytrue = fill(0.6, 5, 5)
#        Ztrue = [-1.5  3.5 1 1
#                  3.5 -1.5 1 1
#                  1    1   1 1
#                  1    1   1 1]
#
#        @test XX ≈ Xtrue atol=1e-2
#        @test YY ≈ Ytrue atol=1e-2
#        @test ZZ ≈ Ztrue atol=1e-2
#        @test ispsd(XX)
#        @test ispsd(1/2*eye(3,3)-XX)
#        @test ispsd(YY+ones(5,5))
#        @test ispsd(2*ones(5,5)-YY)
#        @test ispsd(ones(4,4)-ZZ)
#        @test trace(XX) ≈ 1 atol=1e-4
#        @test trace(YY) ≈ 3 atol=1e-4
#        @test trace(ZZ) ≈ -1 atol=1e-4
#        @test JuMP.objectivevalue(m) ≈ 4.35 atol=1e-2
#
#        @objective(m, Min, X[1,2] + Y[1,2] + Z[1,2])
#        JuMP.solve(m)
#        @test JuMP.terminationstatus(m) == MOI.Success
#        @test JuMP.primalstatus(m) == MOI.FeasiblePoint
#
#        XX, YY, ZZ, = JuMP.resultvalue.(X), JuMP.resultvalue.(Y), JuMP.resultvalue.(Z)
#        Xtrue = [ 0.25 -0.25 0
#                 -0.25  0.25 0
#                  0     0    0.5]
#        Ytrue = fill(0.6, 5, 5)
#        Ztrue = [-1.5 -1.5 1 1
#                 -1.5 -1.5 1 1
#                  1    1   1 1
#                  1    1   1 1]
#
#        @test isapprox(XX, Xtrue, atol=1e-2)
#        @test isapprox(YY, Ytrue, atol=1e-2)
#        @test isapprox(ZZ, Ztrue, atol=1e-2)
#        @test ispsd(XX)
#        @test ispsd(1/2*eye(3,3)-XX)
#        @test ispsd(YY+ones(5,5))
#        @test ispsd(2*ones(5,5)-YY)
#        @test ispsd(ones(4,4)-ZZ)
#        @test isapprox(trace(XX), 1, atol=1e-4)
#        @test isapprox(trace(YY), 3, atol=1e-4)
#        @test isapprox(trace(ZZ), -1, atol=1e-4)
#
#        # Test SDP constraints
#        m = Model(solver=solver)
#        @variable(m, X[1:3,1:3], PSD)
#        @SDconstraint(m, ones(3,3) <= X)
#
#        @objective(m, Min, trace(ones(3,3)*X))
#        JuMP.solve(m)
#        @test JuMP.terminationstatus(m) == MOI.Success
#        @test JuMP.primalstatus(m) == MOI.FeasiblePoint
#
#        XX = JuMP.resultvalue.(X)
#        @test ispsd(XX)
#        @test ispsd(XX - ones(3,3))
#        @test isapprox(JuMP.objectivevalue(m), 9, atol=1e-4)
#
#        # Another test SDP
#        m = Model(solver=solver)
#        @variable(m, X[1:3,1:3], PSD)
#        @variable(m, Y[1:2,1:2], PSD)
#
#        C = eye(3,3)
#        A1 = zeros(3,3)
#        A1[1,1] = 1.0
#        A2 = zeros(3,3)
#        A2[2,2] = 1.0
#        A3 = zeros(3,3)
#        A3[3,3] = 1.0
#        D = eye(2,2)
#        B1 = ones(2,2)
#        B2 = zeros(2,2)
#        B2[1,1] = 1
#        B3 = zeros(2,2)
#        B3[1,2] = B3[2,1] = 2
#
#        @objective(m, Min, trace(C*X)+1+trace(D*Y))
#        @constraint(m, trace(A1*X-eye(3,3)/3) == 0)
#        @constraint(m, 2*trace(A2*X) == 1)
#        @constraint(m, trace(A3*X) >= 2)
#        @constraint(m, trace(B1*Y) == 1)
#        @constraint(m, trace(B2*Y) == 0)
#        @constraint(m, trace(B3*Y) <= 0)
#        @constraint(m, trace(A1*X)+trace(B1*Y) >= 1)
#        @constraint(m, Y[2,2] == 1)
#
#        JuMP.solve(m)
#        @test JuMP.terminationstatus(m) == MOI.Success
#        @test JuMP.primalstatus(m) == MOI.FeasiblePoint
#
#        XX, YY = JuMP.resultvalue.(X), JuMP.resultvalue.(Y)
#        @test isapprox(trace(A1*XX-eye(3,3)/3), 0, atol=1e-5)
#        @test isapprox(2*trace(A2*XX), 1, atol=1e-5)
#        @test trace(A3*XX) >= 2 - 1e-5
#        @test isapprox(trace(B1*YY), 1, atol=1e-5)
#        @test isapprox(trace(B2*YY), 0, atol=1e-5)
#        @test trace(B3*YY) <= 1e-3
#        @test trace(A1*XX)+trace(B1*YY) >= 1
#        @test isapprox(YY[2,2], 1, atol=1e-5)
#        @test isapprox(XX, diagm([1,.5,2]), atol=1e-3)
#        @test isapprox(YY, [0 0;0 1], atol=1e-3)
#    end

#    @testset "SDP with SOC" begin
#        m = Model(solver=CSDPSolver(printlevel=0))
#        @variable(m, X[1:2,1:2], PSD)
#        @variable(m, y[0:2])
#        @constraint(m, norm([y[1],y[2]]) <= y[0])
#        @SDconstraint(m, X <= eye(2))
#        @constraint(m, X[1,1] + X[1,2] == y[1] + y[2])
#        @objective(m, Max, trace(X) - y[0])
#
#        JuMP.solve(m)
#        @test JuMP.terminationstatus(m) == MOI.Success
#        @test JuMP.primalstatus(m) == MOI.FeasiblePoint
#
#        XX, yy = JuMP.resultvalue.(X), JuMP.resultvalue.(y)
#        @test ispsd(XX)
#        @test (yy[0] >= 0)
#        @test (yy[1]^2 + yy[2]^2 <= yy[0]^2 + 1e-4)
#        @test ispsd(eye(2)-XX)
#        @test isapprox(XX[1,1] + XX[1,2], yy[1] + yy[2], atol=1e-4)
#        @test isapprox(XX, eye(2), atol=1e-4)
#        @test isapprox(yy[:], [1/sqrt(2), 0.5, 0.5], atol=1e-4)
#        @test isapprox(JuMP.objectivevalue(m), 1.293, atol=1e-2)
#    end

    # min tr(Y)          max 4x_1 +3x2
    #     Y[2,1] <= 4        [ 0 y1 0    [1 0 0
    #     Y[2,2] >= 3         y1 y2 0  <= 0 1 0
    #     Y >= 0               0  0 0]    0 0 1]
    #                         y1 <= 0 y2 >= 0
    @testset "Test problem #2" begin
        m = Model(solver=CSDPSolver(printlevel=0))
        @variable(m, Y[1:3,1:3], PSD)
        c1 = @constraint(m, Y[2,1] <= 4)
        c2 = @constraint(m, Y[2,2] >= 3)
        @objective(m, Min, trace(Y))

        JuMP.solve(m)
        @test JuMP.terminationstatus(m) == MOI.Success
        @test JuMP.primalstatus(m) == MOI.FeasiblePoint

        @test isapprox(JuMP.objectivevalue(m), 3, atol=1e-5)
        #@test isapprox(getdual(c1), 0, atol=1e-5)
        #@test isapprox(getdual(c2), 1, atol=1e-5)
        #@test isapprox(getdual(Y), [1 0 0; 0 0 0; 0 0 1], atol=1e-5)
    end

    # min Y[1,2]          max y
    #     Y[2,1] = 1         [0   y/2 0     [ 0 .5 0
    #                         y/2 0   0  <=  .5  0 0
    #     Y >= 0              0   0   0]      0  0 0]
    #                         y free
    @testset "Test problem #3" begin
        m = Model(solver=CSDPSolver(printlevel=0))
        @variable(m, x >= 0)
        @variable(m, Y[1:3,1:3], PSD)
        c = @constraint(m, Y[2,1] == 1)
        @objective(m, Min, Y[1,2])

        JuMP.solve(m)
        @test JuMP.terminationstatus(m) == MOI.Success
        @test JuMP.primalstatus(m) == MOI.FeasiblePoint

        JuMP.solve(m)
        @test JuMP.terminationstatus(m) == MOI.Success
        @test JuMP.primalstatus(m) == MOI.FeasiblePoint

        @test isapprox(JuMP.objectivevalue(m), 1, atol=1e-4)
        Yval = JuMP.resultvalue.(Y)
        @test isapprox(Yval[1,2], 1, atol=1e-4)
        @test isapprox(Yval[2,1], 1, atol=1e-4)
        #@test isapprox(getdual(c), 1, atol=1e-5)
        #@test isapprox(getdual(Y), zeros(3,3), atol=1e-4)
    end

    # min x + Y[1,1]          max y + z
    #     Y[2,1] = 1         [0   y/2 0     [1 0 0
    #     x >= 1              y/2 0   0  <=  0 0 0
    #                         0   0   0]     0 0 0]
    #                         z <= 1
    #     Y >= 0              y free
    #     x >= 0              z <= 0
    @testset "Test problem #4" begin
        #solver = fixscs(solver, 2000000)
        m = Model(solver=CSDPSolver(printlevel=0))
        @variable(m, x >= 0)
        @variable(m, Y[1:3,1:3], PSD)
        c1 = @constraint(m, x >= 1)
        c2 = @constraint(m, Y[2,1] == 1)
        @objective(m, Min, x + Y[1,1])

        JuMP.solve(m)
        @test JuMP.terminationstatus(m) == MOI.Success
        @test JuMP.primalstatus(m) == MOI.FeasiblePoint

        @test isapprox(JuMP.objectivevalue(m), 1, atol=1e-3)
        @test isapprox(JuMP.resultvalue(x), 1, atol=1e-5)
        @test isapprox(JuMP.resultvalue.(Y)[1,1], 0, atol=1e-4)
        #@test isapprox(getdual(x), 0, atol=1e-5)
        #@test isapprox(getdual(Y), [1 0 0; 0 0 0; 0 0 0], atol=1e-4)
        #@test isapprox(getdual(c1), 1, atol=1e-5)
        #@test isapprox(getdual(c2), 0, atol=1e-4)
    end

#    function nuclear_norm(model, A)
#        m, n = size(A,1), size(A,2)
#        @variable(model, U[1:m,1:m])
#        @variable(model, V[1:n,1:n])
#        @SDconstraint(model, 0 ⪯ [U A; A' V])
#        return 0.5(trace(U) + trace(V'))
#    end
#
#    @testset "Test problem #5" begin
#        m = Model(solver=CSDPSolver(printlevel=0))
#        @variable(m, Y[1:3,1:3], PSD)
#        @constraint(m, Y[2,1] <= 4)
#        @constraint(m, Y[2,2] >= 3)
#        @constraint(m, Y[3,3] <= 2)
#        @objective(m, Min, nuclear_norm(m, Y))
#
#        JuMP.solve(m)
#        @test JuMP.terminationstatus(m) == MOI.Success
#        @test JuMP.primalstatus(m) == MOI.FeasiblePoint
#
#        @test isapprox(JuMP.objectivevalue(m), 3, atol=1e-5)
#    end

    function operator_norm(model, A)
        m, n = size(A,1), size(A,2)
        @variable(model, t >= 0)
        @SDconstraint(model, [t*eye(n) A; A' eye(n)*t] >= 0)
        return t
    end

    @testset "Test problem #6" begin
        m = Model(solver=CSDPSolver(printlevel=0))
        @variable(m, Y[1:3,1:3])
        @constraint(m, Y[2,1] <= 4)
        @constraint(m, Y[2,2] >= 3)
        @constraint(m, sum(Y) >= 12)
        @objective(m, Min, operator_norm(m, Y))

        JuMP.solve(m)
        @test JuMP.terminationstatus(m) == MOI.Success
        @test JuMP.primalstatus(m) == MOI.FeasiblePoint

        @test isapprox(JuMP.objectivevalue(m), 4, atol=1e-5)
    end

#    function lambda_max(model, A)
#        m, n = size(A,1), size(A,2)
#        @variable(model, t)
#        @SDconstraint(model, speye(n)*t - A ⪰ 0)
#        @SDconstraint(model, A >= 0)
#        return t
#    end
#
#    @testset "Test problem #7" begin
#        m = Model(solver=CSDPSolver(printlevel=0))
#        @variable(m, Y[1:3,1:3])
#        @constraint(m, Y[1,1] >= 4)
#        @objective(m, Min, lambda_max(m, Y))
#
#        JuMP.solve(m)
#        @test JuMP.terminationstatus(m) == MOI.Success
#        @test JuMP.primalstatus(m) == MOI.FeasiblePoint
#
#        @test isapprox(JuMP.objectivevalue(m), 4, atol=1e-5)
#    end

    function lambda_min(model, A)
        m, n = size(A,1), size(A,2)
        @variable(model, t)
        @SDconstraint(model, A - eye(n)*t >= 0)
        @SDconstraint(model, A >= 0)
        return t
    end

    @testset "Test problem #8" begin
        m = Model(solver=CSDPSolver(printlevel=0))
        @variable(m, Y[1:3,1:3], PSD)
        @constraint(m, trace(Y) <= 6)
        @objective(m, Max, lambda_min(m, Y))

        JuMP.solve(m)
        @test JuMP.terminationstatus(m) == MOI.Success
        @test JuMP.primalstatus(m) == MOI.FeasiblePoint

        @test isapprox(JuMP.objectivevalue(m), 2, atol=1e-5)
    end

    function matrix_frac(model, x, P)
        n = size(P,1)
        @variable(model, t)
        @variable(model, M[1:(n+1),1:(n+1)], PSD)
        @constraint(model, M[1:n,1:n] .== P)
        @constraint(model, M[1:n,n+1] .== x)
        @constraint(model, M[n+1,n+1] == t)
        return t
    end

#    @testset "Test problem #9" begin
#        m = Model(solver=CSDPSolver(printlevel=0))
#        n = 3
#        x = [1,2,3]
#        lb = 0.5eye(n)
#        ub = 2eye(n)
#        @variable(m, lb[i,j] <= P[i=1:n,j=1:n] <= ub[i,j])
#        @objective(m, Min, matrix_frac(m, x, P))
#
#        JuMP.solve(m)
#        @test JuMP.terminationstatus(m) == MOI.Success
#        @test JuMP.primalstatus(m) == MOI.FeasiblePoint
#
#        @test isapprox(JuMP.objectivevalue(m), 7, atol=1e-5)
#    end

    @testset "Correlation example" begin
        m = Model(solver=CSDPSolver(printlevel=0))

        @variable(m, X[1:3,1:3], PSD)

        # Diagonal is 1s
        @constraint(m, X[1,1] == 1)
        @constraint(m, X[2,2] == 1)
        @constraint(m, X[3,3] == 1)

        # Bounds on the known correlations
        @constraint(m, X[1,2] >= -0.2)
        @constraint(m, X[1,2] <= -0.1)
        @constraint(m, X[2,3] >=  0.4)
        @constraint(m, X[2,3] <=  0.5)

        # Find upper bound
        @objective(m, Max, X[1,3])
        JuMP.solve(m)
        @test JuMP.terminationstatus(m) == MOI.Success
        @test JuMP.primalstatus(m) == MOI.FeasiblePoint
        @test (+0.8719 <= JuMP.resultvalue.(X)[1,3] <= +0.8720)

#        # Find lower bound
#        @objective(m, Min, X[1,3])
#        JuMP.solve(m)
#        @test JuMP.terminationstatus(m) == MOI.Success
#        @test JuMP.primalstatus(m) == MOI.FeasiblePoint
#        @test (-0.9779 >= JuMP.resultvalue.(X)[1,3] >= -0.9799)
    end

    # min o                    max y + X11
    # Q11 - 1   = Q22        [y-X12-X21  0     [0 0
    #                             0     -y] <=  0 0]
    # [1   Q11
    #  Q11 o  ] >= 0          -X[2,2] = 1
    # Q >= 0                        y free
    # o free                        X <= 0
    @testset "Just another SDP" begin
        m = Model(solver=CSDPSolver(printlevel=0))
        @variable(m, Q[1:2, 1:2], PSD)
        c1 = @constraint(m, Q[1,1] - 1 == Q[2,2])
        @variable(m, objective)
        @SDconstraint(m, [1 Q[1,1]; Q[1,1] objective] ⪰ 0)
        @objective(m, Min, objective)

        JuMP.solve(m)
        @test JuMP.terminationstatus(m) == MOI.Success
        @test JuMP.primalstatus(m) == MOI.FeasiblePoint

        @test JuMP.resultvalue.(Q) ≈ [1 0; 0 0] atol=1e-3
        @test JuMP.objectivevalue(m) ≈ 1 atol=1e-4
        @test JuMP.resultvalue(objective) ≈ 1 atol=1e-4
        #@test getdual(objective) ≈ 0 atol=1e-5
        #@test getdual(Q) ≈ [0 0; 0 2] atol=1e-3
        #@test getdual(c1) ≈ 2 atol=1e-4 # y
        #@test getdual(c2) ≈ [-1 1; 1 -1] atol=1e-3 # X
    end

    # The four following tests are from Example 2.11, Example 2.13 and Example 2.27 of:
    # Blekherman, G., Parrilo, P. A., & Thomas, R. R. (Eds.).
    # Semidefinite optimization and convex algebraic geometry SIAM 2013

    # Example 2.11
    @testset "SDP variable and optimal objective not rational" begin
#       solver = fixscs(solver, 7000000)
        m = Model(solver=CSDPSolver(printlevel=0))
        @variable(m, X[1:2,1:2], PSD)
        c = @constraint(m, X[1,1]+X[2,2] == 1)
        @objective(m, Min, 2*X[1,1]+2*X[1,2])
#       @test all(isnan.(getdual(X)))

        JuMP.solve(m)

        @test JuMP.terminationstatus(m) == MOI.Success
        @test JuMP.primalstatus(m) == MOI.FeasiblePoint

        @test JuMP.objectivevalue(m) ≈ 1-sqrt(2) atol=1e-5
        @test JuMP.resultvalue.(X) ≈ [(2-sqrt(2))/4 -1/(2*sqrt(2)); -1/(2*sqrt(2)) (2+sqrt(2))/4] atol=1e-4
#       @test getdual(X) ≈ [1+sqrt(2) 1; 1 sqrt(2)-1] atol=1e-4
#       @test getdual(c) ≈ 1-sqrt(2) atol=1e-5
    end

    # Example 2.13
    @testset "SDP constraint and optimal objective not rational" begin
        #solver = fixscs(solver, 7000000)
        m = Model(solver=CSDPSolver(printlevel=0))
        @variable(m, y)
        c = @SDconstraint(m, [2-y 1; 1 -y] >= 0)
        @objective(m, Max, y)
        #@test all(isnan, getdual(c))
        JuMP.solve(m)
        @test JuMP.terminationstatus(m) == MOI.Success
        @test JuMP.primalstatus(m) == MOI.FeasiblePoint

        @test isapprox(JuMP.objectivevalue(m), 1-sqrt(2), atol=1e-5)
        @test isapprox(JuMP.resultvalue(y), 1-sqrt(2), atol=1e-5)

        #X = getdual(c)
        #@test isapprox(getdual(c), [(2-sqrt(2))/4 -1/(2*sqrt(2)); -1/(2*sqrt(2)) (2+sqrt(2))/4], atol=1e-4)
        #@test isapprox(getdual(y), 0, atol=1e-5)
    end

    # Example 2.27
    # min X[1,1]   max y
    # 2X[1,2] = 1  [0 y     [1 0
    # X ⪰ 0         y 0] ⪯   0 0]
    # The dual optimal solution is y=0 and there is a primal solution
    # [ eps  1/2
    #   1/2  1/eps]
    # for any eps > 0 however there is no primal solution with objective value 0.
    @testset "SDP with dual solution not attained" begin
        #solver = fixscs(solver, 7000000)
        m = Model(solver=CSDPSolver(printlevel=0))
        @variable(m, y)
        c = @SDconstraint(m, [0 y; y 0] <= [1 0; 0 0])
        @objective(m, Max, y)
        #@test all(isnan, getdual(c))
        JuMP.solve(m)
        @test JuMP.terminationstatus(m) == MOI.Success
        @test JuMP.primalstatus(m) == MOI.FeasiblePoint

        @test isapprox(JuMP.objectivevalue(m), 0, atol=1e-5)
        @test isapprox(JuMP.resultvalue(y), 0, atol=1e-5)

        #X = getdual(c)
        #@test isapprox(X[1,1], 0, atol=1e-5)
        #@test isapprox(X[1,2], 1/2, atol=1e-5)
        #@test isapprox(X[2,1], 1/2, atol=1e-5)
        #@test isapprox(getdual(y), 0, atol=1e-5)
    end

    @testset "SDP with primal solution not attained" begin
#       solver = fixscs(solver, 7000000)
        m = Model(solver=CSDPSolver(printlevel=0))
        @variable(m, X[1:2,1:2], PSD)
        c = @constraint(m, 2*X[1,2] == 1)
        @objective(m, Min, X[1,1])
#       @test all(isnan, getdual(X))
        JuMP.solve(m)
        @test JuMP.terminationstatus(m) == MOI.Success
        @test JuMP.primalstatus(m) == MOI.FeasiblePoint

        @test JuMP.terminationstatus(m) == MOI.Success
        @test JuMP.primalstatus(m) == MOI.FeasiblePoint

        @test JuMP.objectivevalue(m) ≈ 0 atol=1e-5
        Xval = JuMP.resultvalue.(X)
        @test Xval[1,1] ≈ 0 atol=1e-5
        @test Xval[1,2] ≈ 1/2 atol=1e-5
        @test Xval[2,1] ≈ 1/2 atol=1e-5

#       @test isapprox(getdual(X), [1 0; 0 0], atol=1e-4)
#       @test isapprox(getdual(c), 0, atol=1e-5)
    end

#    # min X[1,1]     max y/2+z/2
#    # X[1,2] = 1/2   [0 y     [1 0
#    # X[2,1] = 1/2    z 0] ⪯   0 0]
#    # X+Xᵀ ⪰ 0
#    @testset "SDP with dual solution not attained without symmetric A_i" begin
#        #solver = fixscs(solver, 10000000)
#        m = Model(solver=CSDPSolver(printlevel=0))
#        @variable(m, y)
#        @variable(m, z)
#        c = @SDconstraint(m, [0 y; z 0] <= [1 0; 0 0])
#        @objective(m, Max, y/2+z/2)
#        JuMP.solve(m)
#        @test JuMP.terminationstatus(m) == MOI.Success
#        @test JuMP.primalstatus(m) == MOI.FeasiblePoint
#
#        @test isapprox(JuMP.objectivevalue(m), 0, atol=1e-5)
#        @test isapprox(JuMP.resultvalue(y), 0, atol=1e-5)
#        @test isapprox(JuMP.resultvalue(z), 0, atol=1e-5)
#
#        #X = getdual(c)
#        #@test isapprox(X[1,1], 0, atol=1e-5)
#        #@test isapprox(X[1,2], 1/2, atol=1e-5)
#        #@test isapprox(X[2,1], 1/2, atol=1e-5)
#        #@test isapprox(getdual(y), 0, atol=1e-5)
#        #@test isapprox(getdual(z), 0, atol=1e-5)
#    end

#    # min X[1,1]     max y
#    # X[1,2] = 1     [0 y     [1 0
#    # X[2,1] = 0      z 0] ⪯   0 0]
#    # X+Xᵀ ⪰ 0
#    @testset "SDP with dual solution not attained without symmetry" begin
#        #solver = fixscs(solver, 10000000)
#        m = Model(solver=CSDPSolver(printlevel=0))
#        @variable(m, y)
#        @variable(m, z)
#        c = @SDconstraint(m, [0 y; z 0] <= [1 0; 0 0])
#        @objective(m, Max, y)
#
#        JuMP.solve(m)
#        @test JuMP.terminationstatus(m) == MOI.Success
#        @test JuMP.primalstatus(m) == MOI.FeasiblePoint
#
#        @test isapprox(JuMP.objectivevalue(m), 0, atol=1e-5)
#        @test isapprox(JuMP.resultvalue(y), 0, atol=1e-5)
#        @test isapprox(JuMP.resultvalue(z), 0, atol=1e-5)
#
#        #X = getdual(c)
#        #@test isapprox(X[1,1], 0, atol=1e-5)
#        #@test isapprox(X[1,2], 1, atol=1e-5) # X is not symmetric !
#        #@test isapprox(X[2,1], 0, atol=1e-5)
#        #@test isapprox(getdual(y), 0, atol=1e-5)
#        #@test isapprox(getdual(z), 0, atol=1e-5)
#    end

    @testset "Nonzero dual for a scalar variable with sdp solver" begin
        m = Model(solver=CSDPSolver(printlevel=0))
        @variable(m, x1 >= 0)
        @variable(m, x2 >= 0)
        @variable(m, x3 >= 0)
        # The following constraint could be written as 2 linear constrains
        # but and sdp constraint is used to make it a conic problem
        c = @SDconstraint(m, [2x1-x2-x3 0; 0 x1-x2+x3] >= [3 0; 0 2])
        @objective(m, Min, 2x1 - x2)

        JuMP.solve(m)
        @test JuMP.terminationstatus(m) == MOI.Success
        @test JuMP.primalstatus(m) == MOI.FeasiblePoint

        @test isapprox(JuMP.objectivevalue(m), 10/3, atol=1e-5)
        @test isapprox(JuMP.resultvalue(x1), 5/3, atol=1e-5)
        @test isapprox(JuMP.resultvalue(x2), 0, atol=1e-5)
        @test isapprox(JuMP.resultvalue(x3), 1/3, atol=1e-5)

        #@test isapprox(getdual(c), [-2/3 0; 0 -2/3], atol=1e-4)
        #@test isapprox(getdual(x1), 0, atol=1e-5)
        #@test isapprox(getdual(x2), 1/3, atol=1e-5)
        #@test isapprox(getdual(x3), 0, atol=1e-5)
    end


    @testset "No constraint" begin
        m = Model(solver=CSDPSolver(printlevel=0))
        @variable(m, X[1:3,1:3], PSD)
        @objective(m, Min, trace(X))

        JuMP.solve(m)
        @test JuMP.terminationstatus(m) == MOI.Success
        @test JuMP.primalstatus(m) == MOI.FeasiblePoint

        @test abs(JuMP.objectivevalue(m)) < 1e-5
        @test norm(JuMP.resultvalue.(X)) < 1e-5
        #@test isapprox(getdual(X), eye(3), atol=1e-5)
    end

end

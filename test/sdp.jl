@testset "Semidefinite Programming" begin
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

        m = Model(solver=CSDPSolver(verbose=false))

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

    # min o                    max y + X11
    # Q11 - 1   = Q22        [y-X12-X21  0     [0 0
    #                             0     -y] <=  0 0]
    # [1   Q11
    #  Q11 o  ] >= 0          -X[2,2] = 1
    # Q >= 0                        y free
    # o free                        X <= 0
#   @testset "Just another SDP" begin
#       model = Model(solver=CSDPSolver(verbose=false))
#       @variable(model, Q[1:2, 1:2], PSD)
#       c1 = @constraint(model, Q[1,1] - 1 == Q[2,2])
#       @variable(model, objective)
#       T = [1 Q[1,1]; Q[1,1] objective]
#       @test_throws ErrorException SDConstraint(T, 1)
#       c2 = JuMP.addconstraint(model, SDConstraint(T, 0))
#       @objective(model, Min, objective)
#
#       JuMP.solve(m)
#
#       @test JuMP.terminationstatus(m) == MOI.Success
#       @test JuMP.primalstatus(m) == MOI.FeasiblePoint
#
#       @test JuMP.resultvalue(Q) ≈ [1 0; 0 0] atol=1e-3
#       @test JuMP.objectivevalue(model) ≈ 1 atol=1e-4
#       @test JuMP.resultvalue(objective) ≈ 1 atol=1e-4
#       @test getdual(objective) ≈, 0, atol=1e-5
#       @test getdual(Q) ≈ [0 0; 0 2] atol=1e-3
#       @test getdual(c1) ≈ 2, atol=1e-4 # y
#       @test getdual(c2) ≈ [-1 1; 1 -1] atol=1e-3 # X
#   end

    # The four following tests are from Example 2.11, Example 2.13 and Example 2.27 of:
    # Blekherman, G., Parrilo, P. A., & Thomas, R. R. (Eds.).
    # Semidefinite optimization and convex algebraic geometry SIAM 2013

    # Example 2.11
    @testset "SDP variable and optimal objective not rational" begin
#       solver = fixscs(solver, 7000000)
        m = Model(solver=CSDPSolver(verbose=false))
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

#   # Example 2.13
#   @testset "SDP constraint and optimal objective not rational with $solver" for solver in sdp_solvers
#       solver = fixscs(solver, 7000000)
#       m = Model(solver=solver)
#       @variable(m, y)
#       c = @SDconstraint(m, [2-y 1; 1 -y] >= 0)
#       @objective(m, Max, y)
#       @test all(isnan, getdual(c))
#       status = solve(m)
#
#       @test status == :Optimal
#       @test isapprox(getobjectivevalue(m), 1-sqrt(2), atol=1e-5)
#       @test isapprox(getvalue(y), 1-sqrt(2), atol=1e-5)
#
#       X = getdual(c)
#       @test isapprox(getdual(c), [(2-sqrt(2))/4 -1/(2*sqrt(2)); -1/(2*sqrt(2)) (2+sqrt(2))/4], atol=1e-4)
#       @test isapprox(getdual(y), 0, atol=1e-5)
#   end
#
#   # Example 2.27
#   # min X[1,1]   max y
#   # 2X[1,2] = 1  [0 y     [1 0
#   # X ⪰ 0         y 0] ⪯   0 0]
#   # The dual optimal solution is y=0 and there is a primal solution
#   # [ eps  1/2
#   #   1/2  1/eps]
#   # for any eps > 0 however there is no primal solution with objective value 0.
#   @testset "SDP with dual solution not attained with $solver" for solver in sdp_solvers
#       solver = fixscs(solver, 7000000)
#       m = Model(solver=solver)
#       @variable(m, y)
#       c = @SDconstraint(m, [0 y; y 0] <= [1 0; 0 0])
#       @objective(m, Max, y)
#       @test all(isnan, getdual(c))
#       status = solve(m)
#
#       if contains(string(typeof(solver)),"MosekSolver")
#           # Mosek returns Stall on this instance
#           # Hack until we fix statuses in MPB
#           JuMP.fillConicDuals(m)
#       else
#           @test status == :Optimal
#       end
#       @test isapprox(getobjectivevalue(m), 0, atol=1e-5)
#       @test isapprox(getvalue(y), 0, atol=1e-5)
#
#       X = getdual(c)
#       @test isapprox(X[1,1], 0, atol=1e-5)
#       @test isapprox(X[1,2], 1/2, atol=1e-5)
#       @test isapprox(X[2,1], 1/2, atol=1e-5)
#       @test isapprox(getdual(y), 0, atol=1e-5)
#   end

    @testset "SDP with primal solution not attained" begin
#       solver = fixscs(solver, 7000000)
        m = Model(solver=CSDPSolver(verbose=false))
        @variable(m, X[1:2,1:2], PSD)
        c = @constraint(m, 2*X[1,2] == 1)
        @objective(m, Min, X[1,1])
#       @test all(isnan, getdual(X))
        status = solve(m)

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

    @testset "No constraint" begin
        m = Model(solver=CSDPSolver(verbose=false))
        @variable(m, X[1:3,1:3], PSD)
        @objective(m, Min, trace(X))

        JuMP.solve(m)

        status = solve(m)

        @test JuMP.terminationstatus(m) == MOI.Success
        @test JuMP.primalstatus(m) == MOI.FeasiblePoint

        @test abs(JuMP.objectivevalue(m)) < 1e-5
        @test norm(JuMP.resultvalue.(X)) < 1e-5
        #@test isapprox(getdual(X), eye(3), atol=1e-5)
    end

end

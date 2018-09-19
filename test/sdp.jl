using JuMP, Compat.Test

# If solvers not loaded, load them (i.e running just these tests)
if ((VERSION >= v"0.7-" && !(@isdefined sdp_solvers)) ||
    (VERSION < v"0.7-" && !isdefined(:sdp_solvers)))
    include("solvers.jl")
end

const TOL = 1e-4

@testset "Semidefinite" begin

ispsd(x::Matrix) = minimum(eigvals(x)) ≥ -1e-3
ispsd(x::JuMP.JuMPArray) = ispsd(x.innerArray)

    @testset "Simple SDP with $solver" for solver in sdp_solvers
        contains(string(typeof(solver)),"SCSSolver") && continue
        m = Model(solver=solver)
        @variable(m, X[1:3,1:3], SDP)
        @SDconstraint(m, X <= 1/2*eye(3,3))
        @variable(m, Y[1:5,1:5], Symmetric)
        @SDconstraint(m, -ones(5,5) <= Y)
        @SDconstraint(m, Y <= 2*ones(5,5))
        @variable(m, Z[1:4,1:4], Symmetric)
        @SDconstraint(m, ones(4,4) >= Z)

        @constraint(m, trace(X) == 1)
        @constraint(m, trace(Y) == 3)
        @constraint(m, trace(Z) == -1)

        @test JuMP.numsdconstr(m) == 4

        @objective(m, Max, X[1,2] + Y[1,2] + Z[1,2])
        solve(m)

        XX, YY, ZZ, = getvalue(X), getvalue(Y), getvalue(Z)
        Xtrue = [0.25 0.25 0
                 0.25 0.25 0
                 0    0    0.5]
        Ytrue = fill(0.6, 5, 5)
        Ztrue = [-1.5  3.5 1 1
                  3.5 -1.5 1 1
                  1    1   1 1
                  1    1   1 1]


        @test isapprox(XX, Xtrue, atol=1e-2)
        @test isapprox(YY, Ytrue, atol=1e-2)
        @test isapprox(ZZ, Ztrue, atol=1e-2)
        @test ispsd(XX)
        @test ispsd(1/2*eye(3,3)-XX)
        @test ispsd(YY+ones(5,5))
        @test ispsd(2*ones(5,5)-YY)
        @test ispsd(ones(4,4)-ZZ)
        @test isapprox(trace(XX), 1, atol=1e-4)
        @test isapprox(trace(YY), 3, atol=1e-4)
        @test isapprox(trace(ZZ), -1, atol=1e-4)
        @test isapprox(getobjectivevalue(m), 4.35, atol=1e-2)

        @objective(m, Min, X[1,2] + Y[1,2] + Z[1,2])
        solve(m)

        XX, YY, ZZ, = getvalue(X), getvalue(Y), getvalue(Z)
        Xtrue = [ 0.25 -0.25 0
                 -0.25  0.25 0
                  0     0    0.5]
        Ytrue = fill(0.6, 5, 5)
        Ztrue = [-1.5 -1.5 1 1
                 -1.5 -1.5 1 1
                  1    1   1 1
                  1    1   1 1]

        @test isapprox(XX, Xtrue, atol=1e-2)
        @test isapprox(YY, Ytrue, atol=1e-2)
        @test isapprox(ZZ, Ztrue, atol=1e-2)
        @test ispsd(XX)
        @test ispsd(1/2*eye(3,3)-XX)
        @test ispsd(YY+ones(5,5))
        @test ispsd(2*ones(5,5)-YY)
        @test ispsd(ones(4,4)-ZZ)
        @test isapprox(trace(XX), 1, atol=1e-4)
        @test isapprox(trace(YY), 3, atol=1e-4)
        @test isapprox(trace(ZZ), -1, atol=1e-4)

        # Test SDP constraints
        m = Model(solver=solver)
        @variable(m, X[1:3,1:3], SDP)

        @SDconstraint(m, ones(3,3) <= X)
        @objective(m, Min, trace(ones(3,3)*X))
        stat = solve(m)
        @test stat == :Optimal
        XX = getvalue(X)
        @test ispsd(XX)
        @test ispsd(XX - ones(3,3))
        @test isapprox(getobjectivevalue(m), 9, atol=1e-4)

        # Another test SDP
        m = Model(solver=solver)
        @variable(m, X[1:3,1:3], SDP)
        @variable(m, Y[1:2,1:2], SDP)

        C = eye(3,3)
        A1 = zeros(3,3)
        A1[1,1] = 1.0
        A2 = zeros(3,3)
        A2[2,2] = 1.0
        A3 = zeros(3,3)
        A3[3,3] = 1.0
        D = eye(2,2)
        B1 = ones(2,2)
        B2 = zeros(2,2)
        B2[1,1] = 1
        B3 = zeros(2,2)
        B3[1,2] = B3[2,1] = 2

        @objective(m, Min, trace(C*X)+1+trace(D*Y))
        @constraint(m, trace(A1*X-eye(3,3)/3) == 0)
        @constraint(m, 2*trace(A2*X) == 1)
        @constraint(m, trace(A3*X) >= 2)
        @constraint(m, trace(B1*Y) == 1)
        @constraint(m, trace(B2*Y) == 0)
        @constraint(m, trace(B3*Y) <= 0)
        @constraint(m, trace(A1*X)+trace(B1*Y) >= 1)
        @constraint(m, Y[2,2] == 1)

        stat = solve(m)
        @test stat == :Optimal
        XX, YY = getvalue(X), getvalue(Y)
        @test isapprox(trace(A1*XX-eye(3,3)/3), 0, atol=1e-5)
        @test isapprox(2*trace(A2*XX), 1, atol=1e-5)
        @test trace(A3*XX) >= 2 - 1e-5
        @test isapprox(trace(B1*YY), 1, atol=1e-5)
        @test isapprox(trace(B2*YY), 0, atol=1e-5)
        @test trace(B3*YY) <= 1e-3
        @test trace(A1*XX)+trace(B1*YY) >= 1
        @test isapprox(YY[2,2], 1, atol=1e-5)
        @test isapprox(XX, diagm([1,.5,2]), atol=1e-3)
        @test isapprox(YY, [0 0;0 1], atol=1e-3)
    end

    @testset "Nonsensical SDPs" begin
        m = Model()
        @test_throws ErrorException @variable(m, unequal[1:5,1:6], SDP)
        # Some of these errors happen at compile time, so we can't use @test_throws
        @test macroexpand(:(@variable(m, notone[1:5,2:6], SDP))).head == :error
        @test macroexpand(:(@variable(m, oneD[1:5], SDP))).head == :error
        @test macroexpand(:(@variable(m, threeD[1:5,1:5,1:5], SDP))).head == :error
        @test macroexpand(:(@variable(m, psd[2] <= rand(2,2), SDP))).head == :error
        @test macroexpand(:(@variable(m, -ones(3,4) <= foo[1:4,1:4] <= ones(4,4), SDP))).head == :error
        @test macroexpand(:(@variable(m, -ones(3,4) <= foo[1:4,1:4] <= ones(4,4), Symmetric))).head == :error
        @test macroexpand(:(@variable(m, -ones(4,4) <= foo[1:4,1:4] <= ones(4,5), Symmetric))).head == :error
        @test macroexpand(:(@variable(m, -rand(5,5) <= nonsymmetric[1:5,1:5] <= rand(5,5), Symmetric))).head == :error
    end

    @testset "SDP with SOC with $solver" for solver in sdp_solvers
        m = Model(solver=solver)
        @variable(m, X[1:2,1:2], SDP)
        @variable(m, y[0:2])
        @constraint(m, norm([y[1],y[2]]) <= y[0])
        @SDconstraint(m, X <= eye(2))
        @constraint(m, X[1,1] + X[1,2] == y[1] + y[2])
        @objective(m, Max, trace(X) - y[0])
        stat = solve(m)

        @test stat == :Optimal
        XX, yy = getvalue(X), getvalue(y)
        @test ispsd(XX)
        @test (yy[0] >= 0)
        @test (yy[1]^2 + yy[2]^2 <= yy[0]^2 + 1e-4)
        @test ispsd(eye(2)-XX)
        @test isapprox(XX[1,1] + XX[1,2], yy[1] + yy[2], atol=1e-4)
        @test isapprox(XX, eye(2), atol=1e-4)
        @test isapprox(yy[:], [1/sqrt(2), 0.5, 0.5], atol=1e-4)
        @test isapprox(getobjectivevalue(m), 1.293, atol=1e-2)
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
            @SDconstraint(mod, [V*diagm(x./n)*V' eye(q)[:,i] ; eye(q)[i:i,:] u[i]] >= 0)
        end
        f, A, b, var_cones, con_cones = JuMP.conicdata(mod)
        @test length(f) == 8
        @test size(A) == (21,8)
        @test minimum(abs.(nonzeros(A))) > 0.01
        @test length(b) == 21
        @test var_cones == [(:Free,[1,2,3,4,5,6,7,8])]
        @test con_cones == [(:NonNeg,[1]),(:NonPos,[2,3,4,5,6,7,8,9]),(:SDP,10:15),(:SDP,16:21)]
    end

    # Adapt SDP atom tests from Convex.jl:
    #   https://github.com/JuliaOpt/Convex.jl/blob/master/test/test_sdp.jl
    # facts("[sdp] Test problem #1") do
    # for solver in sdp_solvers
    # context("With solver $(typeof(solver))") do
    #     m = Model(solver=solver)
    #     @variable(m, X[1:2,1:2], SDP)
    #     @objective(m, Max, X[1,1])
    #     stat = solve(m)
    #     @test stat == :Unbounded

    #     setobjectivesense(m, :Min)
    #     stat = solve(m)
    #     @test stat == :Optimal
    #     @test getobjectivevalue(m) == roughly(0, 1e-4)

    #     @constraint(m, X[1,1] == 1)
    #     stat = solve(m)
    #     @test stat == :Optimal
    #     @test getobjectivevalue(m) == roughly(1, 1e-5)
    #     @test getvalue(X[1,1]) == roughly(1, 1e-5)

    #     @objective(m, Min, sum(diag(X)))
    #     stat = solve(m)
    #     @test stat == :Optimal
    #     @test getobjectivevalue(m) == roughly(1, 1e-5)
    # end; end; end

    # min tr(Y)          max 4x_1 +3x2
    #     Y[2,1] <= 4        [ 0 y1 0    [1 0 0
    #     Y[2,2] >= 3         y1 y2 0  <= 0 1 0
    #     Y >= 0               0  0 0]    0 0 1]
    #                         y1 <= 0 y2 >= 0
    @testset "Test problem #2 with $solver" for solver in sdp_solvers
        m = Model(solver=solver)
        @variable(m, Y[1:3,1:3], SDP)
        c1 = @constraint(m, Y[2,1] <= 4)
        c2 = @constraint(m, Y[2,2] >= 3)
        @objective(m, Min, trace(Y))
        stat = solve(m)
        @test stat == :Optimal
        @test isapprox(getobjectivevalue(m), 3, atol=1e-5)
        @test isapprox(getdual(c1), 0, atol=1e-5)
        @test isapprox(getdual(c2), 1, atol=1e-5)
        @test isapprox(getdual(Y), [1 0 0; 0 0 0; 0 0 1], atol=1e-5)
    end

    # min Y[1,2]          max y
    #     Y[2,1] = 1         [0   y/2 0     [ 0 .5 0
    #                         y/2 0   0  <=  .5  0 0
    #     Y >= 0              0   0   0]      0  0 0]
    #                         y free
    @testset "Test problem #3 with $solver" for solver in sdp_solvers
        m = Model(solver=solver)
        @variable(m, x >= 0)
        @variable(m, Y[1:3,1:3], SDP)
        c = @constraint(m, Y[2,1] == 1)
        @objective(m, Min, Y[1,2])
        stat = solve(m)
        @test stat == :Optimal
        @test isapprox(getobjectivevalue(m), 1, atol=1e-4)
        Yval = getvalue(Y)
        @test isapprox(Yval[1,2], 1, atol=1e-4)
        @test isapprox(Yval[2,1], 1, atol=1e-4)
        @test isapprox(getdual(c), 1, atol=1e-5)
        @test isapprox(getdual(Y), zeros(3,3), atol=1e-4)
    end

    # min x + Y[1,1]          max y + z
    #     Y[2,1] = 1         [0   y/2 0     [1 0 0
    #     x >= 1              y/2 0   0  <=  0 0 0
    #                         0   0   0]     0 0 0]
    #                         z <= 1
    #     Y >= 0              y free
    #     x >= 0              z <= 0
    @testset "Test problem #4 with $solver" for solver in sdp_solvers
        solver = fixscs(solver, 2000000)
        m = Model(solver=solver)
        @variable(m, x >= 0)
        @variable(m, Y[1:3,1:3], SDP)
        c1 = @constraint(m, x >= 1)
        c2 = @constraint(m, Y[2,1] == 1)
        @objective(m, Min, x + Y[1,1])
        stat = solve(m)
        @test stat == :Optimal
        @test isapprox(getobjectivevalue(m), 1, atol=1e-3)
        @test isapprox(getvalue(x), 1, atol=1e-5)
        @test isapprox(getvalue(Y)[1,1], 0, atol=1e-4)
        @test isapprox(getdual(x), 0, atol=1e-5)
        @test isapprox(getdual(Y), [1 0 0; 0 0 0; 0 0 0], atol=1e-4)
        @test isapprox(getdual(c1), 1, atol=1e-5)
        @test isapprox(getdual(c2), 0, atol=1e-4)
    end

    # SCS cannot solve it
    #facts("[sdp] Test problem #4.5") do
    #for solver in sdp_solvers
    #context("With solver $(typeof(solver))") do
    #    solver = fixscs(solver, 10000000)
    #    m = Model(solver=solver)
    #    @variable(m, x >= 1)
    #    @variable(m, Y[1:3,1:3], SDP)
    #    c = @constraint(m, Y[2,1] == 1)
    #    @objective(m, Min, x + Y[1,1])
    #    stat = solve(m)
    #    @test stat == :Optimal
    #    @test getobjectivevalue(m) == roughly(1, 1e-3)
    #    @show getdual(x)
    #    @show getdual(Y)
    #    @show getdual(c)
    #end; end; end

    function nuclear_norm(model, A)
        m, n = size(A,1), size(A,2)
        @variable(model, U[1:m,1:m])
        @variable(model, V[1:n,1:n])
        @SDconstraint(model, 0 ⪯ [U A; A' V])
        return 0.5(trace(U) + trace(V'))
    end

    @testset "Test problem #5 with $solver" for solver in sdp_solvers
        m = Model(solver=solver)
        @variable(m, Y[1:3,1:3], SDP)
        @constraint(m, Y[2,1] <= 4)
        @constraint(m, Y[2,2] >= 3)
        @constraint(m, Y[3,3] <= 2)
        @objective(m, Min, nuclear_norm(m, Y))
        stat = solve(m)
        @test stat == :Optimal
        @test isapprox(getobjectivevalue(m), 3, atol=1e-5)
    end

    function operator_norm(model, A)
        m, n = size(A,1), size(A,2)
        @variable(model, t >= 0)
        @SDconstraint(model, [t*eye(n) A; A' eye(n)*t] >= 0)
        return t
    end

    @testset "Test problem #6 with $solver" for solver in sdp_solvers
        m = Model(solver=solver)
        @variable(m, Y[1:3,1:3])
        @constraint(m, Y[2,1] <= 4)
        @constraint(m, Y[2,2] >= 3)
        @constraint(m, sum(Y) >= 12)
        @objective(m, Min, operator_norm(m, Y))
        stat = solve(m)
        @test stat == :Optimal
        @test isapprox(getobjectivevalue(m), 4, atol=1e-5)
    end

    function lambda_max(model, A)
        m, n = size(A,1), size(A,2)
        @variable(model, t)
        @SDconstraint(model, speye(n)*t - A ⪰ 0)
        @SDconstraint(model, A >= 0)
        return t
    end

    @testset "Test problem #7 with $solver" for solver in sdp_solvers
        m = Model(solver=solver)
        @variable(m, Y[1:3,1:3])
        @constraint(m, Y[1,1] >= 4)
        @objective(m, Min, lambda_max(m, Y))
        stat = solve(m)
        @test stat == :Optimal
        @test isapprox(getobjectivevalue(m), 4, atol=1e-5)
    end

    function lambda_min(model, A)
        m, n = size(A,1), size(A,2)
        @variable(model, t)
        @SDconstraint(model, A - eye(n)*t >= 0)
        @SDconstraint(model, A >= 0)
        return t
    end

    @testset "Test problem #8 with $solver" for solver in sdp_solvers
        m = Model(solver=solver)
        @variable(m, Y[1:3,1:3], SDP)
        @constraint(m, trace(Y) <= 6)
        @objective(m, Max, lambda_min(m, Y))
        stat = solve(m)
        @test stat == :Optimal
        @test isapprox(getobjectivevalue(m), 2, atol=1e-5)
    end

    function matrix_frac(model, x, P)
        n = size(P,1)
        @variable(model, t)
        @variable(model, M[1:(n+1),1:(n+1)], SDP)
        @constraint(model, M[1:n,1:n] .== P)
        @constraint(model, M[1:n,n+1] .== x)
        @constraint(model, M[n+1,n+1] == t)
        return t
    end

    @testset "Test problem #9 with $solver" for solver in sdp_solvers
        m = Model(solver=solver)
        n = 3
        x = [1,2,3]
        lb = 0.5eye(n)
        ub = 2eye(n)
        @variable(m, lb[i,j] <= P[i=1:n,j=1:n] <= ub[i,j])
        @objective(m, Min, matrix_frac(m, x, P))
        stat = solve(m)
        @test stat == :Optimal
        @test isapprox(getobjectivevalue(m), 7, atol=1e-5)
    end

    @testset "Correlation example with $solver" for solver in sdp_solvers
        m = Model(solver=solver)

        @variable(m, X[1:3,1:3], SDP)

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
        stat = solve(m)
        @test stat == :Optimal
        @test (+0.8719 <= getvalue(X)[1,3] <= +0.8720)

        # Find lower bound
        @objective(m, Min, X[1,3])
        stat = solve(m)
        @test stat == :Optimal
        @test (-0.9779 >= getvalue(X)[1,3] >= -0.9799)
    end

    @testset "Robust uncertainty example with $solver" for solver in sdp_solvers
        include(joinpath("data","robust_uncertainty.jl"))
        R = 1
        d = 3
        𝛿 = 0.05
        ɛ = 0.05
        N = ceil((2+2log(2/𝛿))^2) + 1

        Γ1(𝛿,N) = (R/sqrt(N))*(2+sqrt(2*log(1/𝛿)))
        Γ2(𝛿,N) = (2R^2/sqrt(N))*(2+sqrt(2*log(2/𝛿)))

        @testset "d = $d" for d in [3,5,8]

            μhat = μhats[d]
            M = Ms[d]
            Σhat = 1/(d-1)*(M-ones(d)*μhat')'*(M-ones(d)*μhat')

            m = Model(solver=solver)

            @variable(m, Σ[1:d,1:d], SDP)
            @variable(m, u[1:d])
            @variable(m, μ[1:d])

            @variable(m, L1[1:d])
            @constraint(m, L1 .== (μ-μhat))
            @constraint(m, norm(L1) <= Γ1(𝛿/2,N))

            @variable(m, L2[1:d,1:d])
            @constraint(m, L2 .== (Σ-Σhat))
            @constraint(m, vecnorm(L2) <= Γ2(𝛿/2,N))

            A = [(1-ɛ)/ɛ (u-μ)';
                 (u-μ)     Σ   ]
            @SDconstraint(m, A >= 0)

            c = cs[d]
            @objective(m, Max, dot(c,u))

            stat = solve(m)

            object = getobjectivevalue(m)
            exact = dot(μhat,c) + Γ1(𝛿/2,N)*norm(c) + sqrt((1-ɛ)/ɛ)*sqrt(dot(c,(Σhat+Γ2(𝛿/2,N)*eye(d,d))*c))
            @test stat == :Optimal
            @test isapprox(object, exact, atol=1e-5)
        end
    end

    @testset "Robust uncertainty example (with norms) with $solver" for solver in sdp_solvers
        include(joinpath("data","robust_uncertainty.jl"))
        R = 1
        d = 3
        𝛿 = 0.05
        ɛ = 0.05
        N = ceil((2+2log(2/𝛿))^2) + 1

        Γ1(𝛿,N) = (R/sqrt(N))*(2+sqrt(2*log(1/𝛿)))
        Γ2(𝛿,N) = (2R^2/sqrt(N))*(2+sqrt(2*log(2/𝛿)))

        @testset "d = $d" for d in [3,5,8]

            μhat = μhats[d]
            M = Ms[d]
            Σhat = 1/(d-1)*(M-ones(d)*μhat')'*(M-ones(d)*μhat')

            m = Model(solver=solver)

            @variable(m, Σ[1:d,1:d], SDP)
            @variable(m, u[1:d])
            @variable(m, μ[1:d])

            @constraint(m, norm(μ-μhat) <= Γ1(𝛿/2,N))
            @constraint(m, vecnorm(Σ-Σhat) <= Γ2(𝛿/2,N))

            A = [(1-ɛ)/ɛ (u-μ)';
                 (u-μ)     Σ   ]
            @SDconstraint(m, A >= 0)

            c = cs[d]
            @objective(m, Max, dot(c,u))

            stat = solve(m)

            object = getobjectivevalue(m)
            exact = dot(μhat,c) + Γ1(𝛿/2,N)*norm(c) + sqrt((1-ɛ)/ɛ)*sqrt(dot(c,(Σhat+Γ2(𝛿/2,N)*eye(d,d))*c))
            @test stat == :Optimal
            @test isapprox(object, exact, atol=1e-5)
        end
    end

    @testset "Can't mix SDP and QP (#665) with $solver" for solver in sdp_solvers
        model = Model(solver=solver)
        @variable(model, Q[1:2, 1:2], SDP)
        @constraint(model, Q[1,1] - 1 == Q[2,2])
        @objective(model, Min, Q[1,1] * Q[1,1])
        @test_throws ErrorException solve(model)
    end

    # min o                    max y + X11
    # Q11 - 1   = Q22        [y-X12-X21  0     [0 0
    #                             0     -y] <=  0 0]
    # [1   Q11
    #  Q11 o  ] >= 0          -X[2,2] = 1
    # Q >= 0                        y free
    # o free                        X <= 0
    @testset "Just another SDP with $solver" for solver in sdp_solvers
        model = Model(solver=solver)
        @variable(model, Q[1:2, 1:2], SDP)
        c1 = @constraint(model, Q[1,1] - 1 == Q[2,2])
        @variable(model, objective)
        T = [1 Q[1,1]; Q[1,1] objective]
        @test_throws ErrorException SDConstraint(T, 1)
        c2 = JuMP.addconstraint(model, SDConstraint(T, 0))
        @objective(model, Min, objective)
        @test solve(model) == :Optimal
        @test isapprox(getvalue(Q), [1 0;0 0], atol=1e-3)
        @test isapprox(getobjectivevalue(model), 1, atol=TOL)
        @test isapprox(getvalue(objective), 1, atol=TOL)
        @test isapprox(getdual(objective), 0, atol=1e-5)
        @test isapprox(getdual(Q), [0 0; 0 2], atol=1e-3)
        @test isapprox(getdual(c1), 2, atol=1e-4) # y
        @test isapprox(getdual(c2), [-1 1; 1 -1], atol=1e-3) # X
    end

    if length(lp_solvers) > 0
        @testset "Internal Model unloaded when SDP constraint added (#830)" begin
            model = Model(solver=lp_solvers[1])
            @variable(model, x)
            solve(model)
            T = [1 x; -x 1]
            c = @SDconstraint(model, T ⪰ 0)
            @test typeof(c) == JuMP.ConstraintRef{JuMP.Model,JuMP.SDConstraint}
            @test model.internalModelLoaded == false
        end
    end

    # The four following tests are from Example 2.11, Example 2.13 and Example 2.27 of:
    # Blekherman, G., Parrilo, P. A., & Thomas, R. R. (Eds.).
    # Semidefinite optimization and convex algebraic geometry SIAM 2013

    # Example 2.11
    @testset "SDP variable and optimal objective not rational with $solver" for solver in sdp_solvers
        solver = fixscs(solver, 7000000)
        m = Model(solver=solver)
        @variable(m, X[1:2,1:2], SDP)
        c = @constraint(m, X[1,1]+X[2,2] == 1)
        @objective(m, Min, 2*X[1,1]+2*X[1,2])
        @test all(isnan.(getdual(X)))
        status = solve(m)

        @test status == :Optimal
        @test isapprox(getobjectivevalue(m), 1-sqrt(2), atol=1e-5)
        @test isapprox(getvalue(X), [(2-sqrt(2))/4 -1/(2*sqrt(2)); -1/(2*sqrt(2)) (2+sqrt(2))/4], atol=1e-4)
        @test isapprox(getdual(X), [1+sqrt(2) 1; 1 sqrt(2)-1], atol=1e-4)
        @test isapprox(getdual(c), 1-sqrt(2), atol=1e-5)
    end

    # Example 2.13
    @testset "SDP constraint and optimal objective not rational with $solver" for solver in sdp_solvers
        solver = fixscs(solver, 7000000)
        m = Model(solver=solver)
        @variable(m, y)
        c = @SDconstraint(m, [2-y 1; 1 -y] >= 0)
        @objective(m, Max, y)
        @test all(isnan, getdual(c))
        status = solve(m)

        @test status == :Optimal
        @test isapprox(getobjectivevalue(m), 1-sqrt(2), atol=1e-5)
        @test isapprox(getvalue(y), 1-sqrt(2), atol=1e-5)

        X = getdual(c)
        @test isapprox(getdual(c), [(2-sqrt(2))/4 -1/(2*sqrt(2)); -1/(2*sqrt(2)) (2+sqrt(2))/4], atol=1e-4)
        @test isapprox(getdual(y), 0, atol=1e-5)
    end

    # Example 2.27
    # min X[1,1]   max y
    # 2X[1,2] = 1  [0 y     [1 0
    # X ⪰ 0         y 0] ⪯   0 0]
    # The dual optimal solution is y=0 and there is a primal solution
    # [ eps  1/2
    #   1/2  1/eps]
    # for any eps > 0 however there is no primal solution with objective value 0.
    @testset "SDP with dual solution not attained with $solver" for solver in sdp_solvers
        solver = fixscs(solver, 7000000)
        m = Model(solver=solver)
        @variable(m, y)
        c = @SDconstraint(m, [0 y; y 0] <= [1 0; 0 0])
        @objective(m, Max, y)
        @test all(isnan, getdual(c))
        status = solve(m)

        if contains(string(typeof(solver)),"MosekSolver")
            # Mosek returns Stall on this instance
            # Hack until we fix statuses in MPB
            JuMP.fillConicDuals(m)
        else
            @test status == :Optimal
        end
        @test isapprox(getobjectivevalue(m), 0, atol=1e-5)
        @test isapprox(getvalue(y), 0, atol=1e-5)

        X = getdual(c)
        @test isapprox(X[1,1], 0, atol=1e-5)
        @test isapprox(X[1,2], 1/2, atol=1e-5)
        @test isapprox(X[2,1], 1/2, atol=1e-5)
        @test isapprox(getdual(y), 0, atol=1e-5)
    end

    @testset "SDP with primal solution not attained with $solver" for solver in sdp_solvers
        solver = fixscs(solver, 7000000)
        m = Model(solver=solver)
        @variable(m, X[1:2,1:2], SDP)
        c = @constraint(m, 2*X[1,2] == 1)
        @objective(m, Min, X[1,1])
        @test all(isnan, getdual(X))
        status = solve(m)

        @test status == :Optimal
        @test isapprox(getobjectivevalue(m), 0, atol=1e-5)
        Xval = getvalue(X)
        @test isapprox(Xval[1,1], 0, atol=1e-5)
        @test isapprox(Xval[1,2], 1/2, atol=1e-5)
        @test isapprox(Xval[2,1], 1/2, atol=1e-5)

        @test isapprox(getdual(X), [1 0; 0 0], atol=1e-4)
        @test isapprox(getdual(c), 0, atol=1e-5)
    end

    # min X[1,1]     max y/2+z/2
    # X[1,2] = 1/2   [0 y     [1 0
    # X[2,1] = 1/2    z 0] ⪯   0 0]
    # X+Xᵀ ⪰ 0
    @testset "SDP with dual solution not attained without symmetric A_i with $solver" for solver in sdp_solvers
        solver = fixscs(solver, 10000000)
        m = Model(solver=solver)
        @variable(m, y)
        @variable(m, z)
        c = @SDconstraint(m, [0 y; z 0] <= [1 0; 0 0])
        @objective(m, Max, y/2+z/2)
        status = solve(m)

        if contains(string(typeof(solver)),"MosekSolver")
            # Mosek returns Stall on this instance
            # Hack until we fix statuses in MPB
            JuMP.fillConicDuals(m)
        else
            @test status == :Optimal
        end
        @test isapprox(getobjectivevalue(m), 0, atol=1e-5)
        @test isapprox(getvalue(y), 0, atol=1e-5)
        @test isapprox(getvalue(z), 0, atol=1e-5)

        X = getdual(c)
        @test isapprox(X[1,1], 0, atol=1e-5)
        @test isapprox(X[1,2], 1/2, atol=1e-5)
        @test isapprox(X[2,1], 1/2, atol=1e-5)
        @test isapprox(getdual(y), 0, atol=1e-5)
        @test isapprox(getdual(z), 0, atol=1e-5)
    end

    # min X[1,1]     max y
    # X[1,2] = 1     [0 y     [1 0
    # X[2,1] = 0      z 0] ⪯   0 0]
    # X+Xᵀ ⪰ 0
    @testset "SDP with dual solution not attained without symmetry with $solver" for solver in sdp_solvers
        solver = fixscs(solver, 10000000)
        m = Model(solver=solver)
        @variable(m, y)
        @variable(m, z)
        c = @SDconstraint(m, [0 y; z 0] <= [1 0; 0 0])
        @objective(m, Max, y)
        status = solve(m)

        if contains(string(typeof(solver)),"MosekSolver")
            # Mosek returns Stall on this instance
            # Hack until we fix statuses in MPB
            JuMP.fillConicDuals(m)
        else
            @test status == :Optimal
        end
        @test isapprox(getobjectivevalue(m), 0, atol=1e-5)
        @test isapprox(getvalue(y), 0, atol=1e-5)
        @test isapprox(getvalue(z), 0, atol=1e-5)

        X = getdual(c)
        @test isapprox(X[1,1], 0, atol=1e-5)
        @test isapprox(X[1,2], 1, atol=1e-5) # X is not symmetric !
        @test isapprox(X[2,1], 0, atol=1e-5)
        @test isapprox(getdual(y), 0, atol=1e-5)
        @test isapprox(getdual(z), 0, atol=1e-5)
    end

    @testset "Nonzero dual for a scalar variable with sdp solver with $solver" for solver in sdp_solvers
        m = Model(solver=solver)
        @variable(m, x1 >= 0)
        @variable(m, x2 >= 0)
        @variable(m, x3 >= 0)
        # The following constraint could be written as 2 linear constrains
        # but and sdp constraint is used to make it a conic problem
        c = @SDconstraint(m, [2x1-x2-x3 0; 0 x1-x2+x3] >= [3 0; 0 2])
        @objective(m, Min, 2x1 - x2)
        status = solve(m)

        @test status == :Optimal
        @test isapprox(getobjectivevalue(m), 10/3, atol=1e-5)
        @test isapprox(getvalue(x1), 5/3, atol=1e-5)
        @test isapprox(getvalue(x2), 0, atol=1e-5)
        @test isapprox(getvalue(x3), 1/3, atol=1e-5)

        @test isapprox(getdual(c), [-2/3 0; 0 -2/3], atol=1e-4)
        @test isapprox(getdual(x1), 0, atol=1e-5)
        @test isapprox(getdual(x2), 1/3, atol=1e-5)
        @test isapprox(getdual(x3), 0, atol=1e-5)
    end

    @testset "No constraint with $solver" for solver in sdp_solvers
        println(solver)
        m = Model(solver=solver)
        @variable(m, X[1:3,1:3], SDP)
        @objective(m, Min, trace(X))
        status = solve(m)

        @test status == :Optimal
        @test abs(getobjectivevalue(m)) < 1e-5
        @test norm(getvalue(X)) < 1e-5
        @test isapprox(getdual(X), eye(3), atol=1e-5)
    end
end

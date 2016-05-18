ispsd(x::Matrix) = minimum(eigvals(x)) ≥ -1e-3
ispsd(x::JuMP.JuMPArray) = ispsd(x.innerArray)

facts("[sdp] Test simple SDP") do
for solver in sdp_solvers
context("With solver $(typeof(solver))") do
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

    @fact norm(XX-Xtrue) --> roughly(0, 1e-2)
    @fact norm(YY-Ytrue) --> roughly(0, 1e-2)
    @fact norm(ZZ-Ztrue) --> roughly(0, 1e-2)
    @fact ispsd(XX) --> true
    @fact ispsd(1/2*eye(3,3)-XX) --> true
    @fact ispsd(YY+ones(5,5)) --> true
    @fact ispsd(2*ones(5,5)-YY) --> true
    @fact ispsd(ones(4,4)-ZZ) --> true
    @fact trace(XX) --> roughly( 1, 1e-4)
    @fact trace(YY) --> roughly( 3, 1e-4)
    @fact trace(ZZ) --> roughly(-1, 1e-4)
    @fact getobjectivevalue(m) --> roughly(4.35, 1e-2)

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

    @fact norm(XX-Xtrue) --> roughly(0, 1e-2)
    @fact norm(YY-Ytrue) --> roughly(0, 1e-2)
    @fact norm(ZZ-Ztrue) --> roughly(0, 1e-2)
    @fact ispsd(XX) --> true
    @fact ispsd(1/2*eye(3,3)-XX) --> true
    @fact ispsd(YY+ones(5,5)) --> true
    @fact ispsd(2*ones(5,5)-YY) --> true
    @fact ispsd(ones(4,4)-ZZ) --> true
    @fact trace(XX) --> roughly( 1, 1e-4)
    @fact trace(YY) --> roughly( 3, 1e-4)
    @fact trace(ZZ) --> roughly(-1, 1e-4)

    # Test SDP constraints
    m = Model(solver=solver)
    @variable(m, X[1:3,1:3], SDP)

    @SDconstraint(m, ones(3,3) <= X)
    @objective(m, Min, trace(ones(3,3)*X))
    stat = solve(m)
    @fact stat --> :Optimal
    XX = getvalue(X)
    @fact ispsd(XX) --> true
    @fact ispsd(XX - ones(3,3)) --> true
    @fact getobjectivevalue(m) --> roughly(9,1e-4)

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
    @fact stat --> :Optimal
    XX, YY = getvalue(X), getvalue(Y)
    @fact trace(A1*XX-eye(3,3)/3) --> roughly(0, 1e-5)
    @fact 2*trace(A2*XX) --> roughly(1, 1e-5)
    @fact (trace(A3*XX) >= 2 - 1e-5) --> true
    @fact trace(B1*YY) --> roughly(1, 1e-5)
    @fact trace(B2*YY) --> roughly(0, 1e-5)
    @fact (trace(B3*YY) <= 1e-3) --> true
    @fact (trace(A1*XX)+trace(B1*YY) >= 1) --> true
    @fact YY[2,2] --> roughly(1, 1e-5)
    @fact norm(XX - diagm([1,.5,2])) --> roughly(0, 1e-3)
    @fact norm(YY - [0 0;0 1]) --> roughly(0, 1e-3)
end; end; end

facts("[sdp] Nonsensical SDPs") do
    m = Model()
    @fact_throws @variable(m, unequal[1:5,1:6], SDP)
    # Some of these errors happen at compile time, so we can't use @fact_throws
    @fact macroexpand(:(@variable(m, notone[1:5,2:6], SDP))).head --> :error
    @fact macroexpand(:(@variable(m, oneD[1:5], SDP))).head --> :error
    @fact macroexpand(:(@variable(m, threeD[1:5,1:5,1:5], SDP))).head --> :error
    @fact macroexpand(:(@variable(m, psd[2] <= rand(2,2), SDP))).head --> :error
    @fact macroexpand(:(@variable(m, -ones(3,4) <= foo[1:4,1:4] <= ones(4,4), SDP))).head --> :error
    @fact_throws @variable(m, -ones(3,4) <= foo[1:4,1:4] <= ones(4,4), Symmetric)
    @fact_throws @variable(m, -ones(4,4) <= foo[1:4,1:4] <= ones(4,5), Symmetric)
    @fact_throws @variable(m, -rand(5,5) <= nonsymmetric[1:5,1:5] <= rand(5,5), Symmetric)
end

facts("[sdp] SDP with quadratics") do
for solver in sdp_solvers
context("With solver $(typeof(solver))") do
    m = Model(solver=solver)
    @variable(m, X[1:2,1:2], SDP)
    @variable(m, y[0:2])
    @constraint(m, y[0] >= 0)
    @constraint(m, y[1]^2 + y[2]^2 <= y[0]^2)
    @SDconstraint(m, X <= eye(2))
    @constraint(m, X[1,1] + X[1,2] == y[1] + y[2])
    @objective(m, Max, trace(X) - y[0])
    stat = solve(m)

    @fact stat --> :Optimal
    XX, yy = getvalue(X), getvalue(y)
    @fact ispsd(XX) --> true
    @fact (yy[0] >= 0) --> true
    @fact (yy[1]^2 + yy[2]^2 <= yy[0]^2 + 1e-4) --> true
    @fact ispsd(eye(2)-XX) --> true
    @fact ((XX[1,1] + XX[1,2]) - (yy[1] + yy[2])) --> roughly(0,1e-4)
    @fact norm(XX - eye(2)) --> roughly(0, 1e-4)
    @fact norm(yy[:] - [1/sqrt(2), 0.5, 0.5]) --> roughly(0, 1e-4)
    @fact getobjectivevalue(m) --> roughly(1.293, 1e-2)
end; end; end

facts("[sdp] Trivial symmetry constraints are removed (#766)") do
for solver in sdp_solvers
context("With solver $(typeof(solver))") do
    q = 2
    m = 3
    angles1 = linspace(3*pi/4, pi, m)
    angles2 = linspace(0, -pi/2, m)
    V = [3.*cos(angles1)' 1.5.*cos(angles2)';
         3.*sin(angles1)' 1.5.*sin(angles2)']
    p = 2*m
    n = 100

    mod = Model(solver=solver)
    @variable(mod, x[j=1:p] >= 0, Int)
    @variable(mod, u[i=1:q] >= 0)
    @objective(mod, Min, sum(u))
    @constraint(mod, sum(x) <= n)
    for i=1:q
        @SDconstraint(mod, [V*diagm(x./n)*V' eye(q)[:,i] ; eye(q)[i,:] u[i]] >= 0)
    end
    f, A, b, var_cones, con_cones = JuMP.conicdata(mod)
    @fact length(f) --> 8
    @fact size(A) --> (23,8)
    @fact length(b) --> 23
    @fact var_cones --> [(:NonNeg,[1,2,3,4,5,6,7,8])]
    @fact con_cones --> [(:NonNeg,[1]),(:NonPos,[2,3,4,5,6,7,8,9]),(:SDP,10:15),(:Zero,16:16),(:SDP,17:22),(:Zero,23:23)]
end; end; end

# Adapt SDP atom tests from Convex.jl:
#   https://github.com/JuliaOpt/Convex.jl/blob/master/test/test_sdp.jl
# facts("[sdp] Test problem #1") do
# for solver in sdp_solvers
# context("With solver $(typeof(solver))") do
#     m = Model(solver=solver)
#     @variable(m, X[1:2,1:2], SDP)
#     @objective(m, Max, X[1,1])
#     stat = solve(m)
#     @fact stat --> :Unbounded

#     setobjectivesense(m, :Min)
#     stat = solve(m)
#     @fact stat --> :Optimal
#     @fact getobjectivevalue(m) --> roughly(0, 1e-4)

#     @constraint(m, X[1,1] == 1)
#     stat = solve(m)
#     @fact stat --> :Optimal
#     @fact getobjectivevalue(m) --> roughly(1, 1e-5)
#     @fact getvalue(X[1,1]) --> roughly(1, 1e-5)

#     @objective(m, Min, sum(diag(X)))
#     stat = solve(m)
#     @fact stat --> :Optimal
#     @fact getobjectivevalue(m) --> roughly(1, 1e-5)
# end; end; end

facts("[sdp] Test problem #2") do
for solver in sdp_solvers
context("With solver $(typeof(solver))") do
    m = Model(solver=solver)
    @variable(m, Y[1:3,1:3], SDP)
    @constraint(m, Y[2,1] <= 4)
    @constraint(m, Y[2,2] >= 3)
    @objective(m, Min, trace(Y))
    stat = solve(m)
    @fact stat --> :Optimal
    @fact getobjectivevalue(m) --> roughly(3, 1e-5)
end; end; end

facts("[sdp] Test problem #3") do
for solver in sdp_solvers
context("With solver $(typeof(solver))") do
    m = Model(solver=solver)
    @variable(m, x >= 0)
    @variable(m, Y[1:3,1:3], SDP)
    @constraint(m, Y[2,1] == 1)
    @objective(m, Min, Y[1,2])
    stat = solve(m)
    # @fact stat --> :Optimal
    @fact getobjectivevalue(m) --> roughly(1, 1e-4)
end; end; end

facts("[sdp] Test problem #4") do
for solver in sdp_solvers
context("With solver $(typeof(solver))") do
    m = Model(solver=solver)
    @variable(m, x >= 0)
    @variable(m, Y[1:3,1:3], SDP)
    @constraint(m, x >= 1)
    @constraint(m, Y[2,1] == 1)
    @objective(m, Min, x + Y[1,1])
    stat = solve(m)
    # @fact stat --> :Optimal # TODO: remove this once SCS starts behaving
    @fact getobjectivevalue(m) --> roughly(1, 1e-3)
end; end; end

function nuclear_norm(model, A)
    m, n = size(A,1), size(A,2)
    @variable(model, U[1:m,1:m])
    @variable(model, V[1:n,1:n])
    @SDconstraint(model, 0 ⪯ [U A; A' V])
    return 0.5(trace(U) + trace(V'))
end

facts("[sdp] Test problem #5") do
for solver in sdp_solvers
context("With solver $(typeof(solver))") do
    m = Model(solver=solver)
    @variable(m, Y[1:3,1:3], SDP)
    @constraint(m, Y[2,1] <= 4)
    @constraint(m, Y[2,2] >= 3)
    @constraint(m, Y[3,3] <= 2)
    @objective(m, Min, nuclear_norm(m, Y))
    stat = solve(m)
    @fact stat --> :Optimal
    @fact getobjectivevalue(m) --> roughly(3, 1e-5)
end; end; end

function operator_norm(model, A)
    m, n = size(A,1), size(A,2)
    @variable(model, t >= 0)
    @SDconstraint(model, [t*eye(n) A; A' eye(n)*t] >= 0)
    return t
end

facts("[sdp] Test problem #6") do
for solver in sdp_solvers
context("With solver $(typeof(solver))") do
    m = Model(solver=solver)
    @variable(m, Y[1:3,1:3])
    @constraint(m, Y[2,1] <= 4)
    @constraint(m, Y[2,2] >= 3)
    @constraint(m, sum(Y) >= 12)
    @objective(m, Min, operator_norm(m, Y))
    stat = solve(m)
    @fact stat --> :Optimal
    @fact getobjectivevalue(m) --> roughly(4, 1e-5)
end; end; end

function lambda_max(model, A)
    m, n = size(A,1), size(A,2)
    @variable(model, t)
    @SDconstraint(model, speye(n)*t - A ⪰ 0)
    @SDconstraint(model, A >= 0)
    return t
end

facts("[sdp] Test problem #7") do
for solver in sdp_solvers
context("With solver $(typeof(solver))") do
    m = Model(solver=solver)
    @variable(m, Y[1:3,1:3])
    @constraint(m, Y[1,1] >= 4)
    @objective(m, Min, lambda_max(m, Y))
    stat = solve(m)
    @fact stat --> :Optimal
    @fact getobjectivevalue(m) --> roughly(4, 1e-5)
end; end; end

function lambda_min(model, A)
    m, n = size(A,1), size(A,2)
    @variable(model, t)
    @SDconstraint(model, A - eye(n)*t >= 0)
    @SDconstraint(model, A >= 0)
    return t
end

facts("[sdp] Test problem #8") do
for solver in sdp_solvers
context("With solver $(typeof(solver))") do
    m = Model(solver=solver)
    @variable(m, Y[1:3,1:3], SDP)
    @constraint(m, trace(Y) <= 6)
    @objective(m, Max, lambda_min(m, Y))
    stat = solve(m)
    @fact stat --> :Optimal
    @fact getobjectivevalue(m) --> roughly(2, 1e-5)
end; end; end

function matrix_frac(model, x, P)
    n = size(P,1)
    @variable(model, t)
    @variable(model, M[1:(n+1),1:(n+1)], SDP)
    @constraint(model, M[1:n,1:n] .== P)
    @constraint(model, M[1:n,n+1] .== x)
    @constraint(model, M[n+1,n+1] == t)
    return t
end

facts("[sdp] Test problem #9") do
for solver in sdp_solvers
context("With solver $(typeof(solver))") do
    m = Model(solver=solver)
    n = 3
    x = [1,2,3]
    lb = 0.5eye(n)
    ub = 2eye(n)
    @variable(m, lb[i,j] <= P[i=1:n,j=1:n] <= ub[i,j])
    @objective(m, Min, matrix_frac(m, x, P))
    stat = solve(m)
    @fact stat --> :Optimal
    @fact getobjectivevalue(m) --> roughly(7, 1e-5)
end; end; end

facts("[sdp] Correlation example") do
for solver in sdp_solvers
context("With solver $(typeof(solver))") do
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
    @fact stat --> :Optimal
    @fact (+0.8719 <= getvalue(X)[1,3] <= +0.8720) --> true

    # Find lower bound
    @objective(m, Min, X[1,3])
    stat = solve(m)
    @fact stat --> :Optimal
    @fact (-0.9779 >= getvalue(X)[1,3] >= -0.9799) --> true
end; end; end

facts("[sdp] Robust uncertainty example") do
for solver in sdp_solvers
context("With solver $(typeof(solver))") do
    include(joinpath("data","robust_uncertainty.jl"))
    R = 1
    d = 3
    𝛿 = 0.05
    ɛ = 0.05
    N = ceil((2+2log(2/𝛿))^2) + 1

    Γ1(𝛿,N) = (R/sqrt(N))*(2+sqrt(2*log(1/𝛿)))
    Γ2(𝛿,N) = (2R^2/sqrt(N))*(2+sqrt(2*log(2/𝛿)))

    for d in [3,5,8]; context("d = $d") do

        μhat = μhats[d]
        M = Ms[d]
        Σhat = 1/(d-1)*(M-ones(d)*μhat')'*(M-ones(d)*μhat')

        m = Model(solver=solver)

        @variable(m, Σ[1:d,1:d], SDP)
        @variable(m, u[1:d])
        @variable(m, μ[1:d])

        @variable(m, t1 >= 0)
        @variable(m, L1[1:d])
        @constraint(m, L1 .== (μ-μhat))
        @constraint(m, sum{L1[i]^2, i=1:d} <= t1^2)
        @constraint(m, t1 <= Γ1(𝛿/2,N))

        @variable(m, t2 >= 0)
        @variable(m, L2[1:d,1:d])
        @constraint(m, L2 .== (Σ-Σhat))
        @constraint(m, sum{L2[i,j]^2, i=1:d, j=1:d} <= t2^2)
        @constraint(m, t2 <= Γ2(𝛿/2,N))

        A = [(1-ɛ)/ɛ (u-μ)';
             (u-μ)     Σ   ]
        @SDconstraint(m, A >= 0)

        c = cs[d]
        @objective(m, Max, dot(c,u))

        stat = solve(m)

        object = getobjectivevalue(m)
        exact = dot(μhat,c) + Γ1(𝛿/2,N)*norm(c) + sqrt((1-ɛ)/ɛ)*sqrt(dot(c,(Σhat+Γ2(𝛿/2,N)*eye(d,d))*c))
        @fact stat --> :Optimal
        @fact abs(object - exact) --> roughly(0, 1e-5)
    end; end
end; end; end

facts("[sdp] Robust uncertainty example (with norms)") do
for solver in sdp_solvers
context("With solver $(typeof(solver))") do
    include(joinpath("data","robust_uncertainty.jl"))
    R = 1
    d = 3
    𝛿 = 0.05
    ɛ = 0.05
    N = ceil((2+2log(2/𝛿))^2) + 1

    Γ1(𝛿,N) = (R/sqrt(N))*(2+sqrt(2*log(1/𝛿)))
    Γ2(𝛿,N) = (2R^2/sqrt(N))*(2+sqrt(2*log(2/𝛿)))

    for d in [3,5,8]; context("d = $d") do

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
        @fact stat --> :Optimal
        @fact abs(object - exact) --> roughly(0, 1e-5)
    end; end
end; end; end

facts("[sdp] Can't mix SDP and QP (#665)") do
for solver in sdp_solvers
context("With solver $(typeof(solver))") do
    model = Model(solver=solver)
    @variable(model, Q[1:2, 1:2], SDP)
    @constraint(model, Q[1,1] - 1 == Q[2,2])
    @objective(model, Min, Q[1,1] * Q[1,1])
    @fact_throws solve(model)
end; end; end

facts("[sdp] Just another SDP") do
for solver in sdp_solvers
context("With solver $(typeof(solver))") do
    model = Model(solver=solver)
    @variable(model, Q[1:2, 1:2], SDP)
    @constraint(model, Q[1,1] - 1 == Q[2,2])
    @variable(model, objective)
    T = [1 Q[1,1]; Q[1,1] objective]
    @SDconstraint(model, T ⪰ 0)
    @objective(model, Min, objective)
    @fact solve(model) --> :Optimal
    @fact getvalue(Q) --> roughly([1 0;0 0], 1e-3)
    @fact getobjectivevalue(model) --> roughly(1, TOL)
end; end; end
